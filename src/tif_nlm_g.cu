/*
 * tif_nlm_g.cu - 3D Non-Local Means Filter for TIFF Image Stacks (CUDA)
 *
 * Compile:
 *   Windows: nvcc -O3 -I%CUDAINCL% -o tif_nlm_g.exe tif_nlm_g.cu -use_fast_math -Xcompiler "/wd 4819" libtiff.lib
 *   Linux:   nvcc -O3 -o tif_nlm_g tif_nlm_g.cu -use_fast_math -ltiff
 *
 * Usage:
 *   tif_nlm_g <input_dir> <output_dir> [patch_radius] [search_radius] [h]
 *
 *   patch_radius:  Half-size of comparison patch (default: 1 -> 3x3x3)
 *   search_radius: Half-size of search window   (default: 3 -> 7x7x7)
 *   h:             Filtering strength            (default: auto from noise sigma)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <cuda_runtime.h>
#include "cuda13_compat.h"

#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
    #include <direct.h>
    #include "tiffio.h"
    #define PATH_SEPARATOR "\\"
    #define snprintf _snprintf
    #ifndef isfinite
        #define isfinite(x) _finite(x)
    #endif
#else
    #include <unistd.h>
    #include <dirent.h>
    #include <sys/stat.h>
    #include <sys/sysinfo.h>
    #include <tiffio.h>
    #define PATH_SEPARATOR "/"
#endif

#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

#define MAX_PATH_LENGTH 1024
#define MAX_FILES 100000
#define BLOCK_SIZE_X 8
#define BLOCK_SIZE_Y 8
#define BLOCK_SIZE_Z 4
#define GPU_MEMORY_FRACTION 0.7f

typedef struct {
    int patch_radius;
    int search_radius;
    float h;
} FilterParams;

typedef struct {
    unsigned int width;
    unsigned int height;
    unsigned int depth;
    unsigned short bits_per_sample;
    unsigned short samples_per_pixel;
    unsigned short sample_format;
    size_t bytes_per_pixel;
    size_t bytes_per_slice;
} ImageInfo;

typedef struct {
    int start_z;
    int end_z;
    int chunk_depth;
    int valid_start;
    int valid_end;
    void *h_data;
    void *d_data;
    void *d_output;
} ChunkData;

/* ----------------------------------------------------------------
 * 3D NLM CUDA kernel
 *
 * For each voxel (x,y,z) in the valid region:
 *   1. Extract patch centered at (x,y,z) with radius patch_radius
 *   2. Search all voxels within search_radius
 *   3. Compute patch distance (normalized SSD)
 *   4. Weight = exp(-norm_ssd / h^2)
 *   5. Output = weighted average of center values
 *   6. Self-weight = max weight among all neighbors
 * ---------------------------------------------------------------- */
template<typename T>
__global__ void nlm_filter_3d_kernel(
    const T* __restrict__ input,
    T* __restrict__ output,
    int width, int height, int chunk_depth,
    int valid_start, int valid_end,
    int patch_radius, int search_radius,
    float h_sq_inv,
    float max_value)
{
    ENABLE_SMEM_SPILLING();
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int z = valid_start + blockIdx.z * blockDim.z + threadIdx.z;

    if (x >= width || y >= height || z > valid_end) return;

    const size_t slice_size = (size_t)width * height;
    const size_t center_idx = (size_t)z * slice_size + (size_t)y * width + x;

    float weighted_sum = 0.0f;
    float weight_sum = 0.0f;
    float max_weight = 0.0f;

    for (int sz = -search_radius; sz <= search_radius; sz++) {
        int nz = z + sz;
        if (nz < 0 || nz >= chunk_depth) continue;

        for (int sy = -search_radius; sy <= search_radius; sy++) {
            int ny = y + sy;
            if (ny < 0 || ny >= height) continue;

            for (int sx = -search_radius; sx <= search_radius; sx++) {
                int nx = x + sx;
                if (nx < 0 || nx >= width) continue;

                /* Skip self */
                if (sx == 0 && sy == 0 && sz == 0) continue;

                /* Compute patch distance (SSD) */
                float ssd = 0.0f;
                int valid_count = 0;

                for (int pz = -patch_radius; pz <= patch_radius; pz++) {
                    int cz = z + pz;
                    int qz = nz + pz;
                    if (cz < 0 || cz >= chunk_depth || qz < 0 || qz >= chunk_depth) continue;

                    for (int py = -patch_radius; py <= patch_radius; py++) {
                        int cy = y + py;
                        int qy = ny + py;
                        if (cy < 0 || cy >= height || qy < 0 || qy >= height) continue;

                        for (int px = -patch_radius; px <= patch_radius; px++) {
                            int cx = x + px;
                            int qx = nx + px;
                            if (cx < 0 || cx >= width || qx < 0 || qx >= width) continue;

                            size_t ci = (size_t)cz * slice_size + (size_t)cy * width + cx;
                            size_t qi = (size_t)qz * slice_size + (size_t)qy * width + qx;
                            float diff = (float)input[ci] - (float)input[qi];
                            ssd += diff * diff;
                            valid_count++;
                        }
                    }
                }

                float norm_ssd = (valid_count > 0) ? ssd / (float)valid_count : 0.0f;

                float weight = __expf(-norm_ssd * h_sq_inv);

                if (weight > max_weight) max_weight = weight;

                size_t neighbor_idx = (size_t)nz * slice_size + (size_t)ny * width + nx;
                weighted_sum += weight * (float)input[neighbor_idx];
                weight_sum += weight;
            }
        }
    }

    /* Self-weight: use max weight among all neighbors */
    float self_value = (float)input[center_idx];
    weighted_sum += max_weight * self_value;
    weight_sum += max_weight;

    float result;
    if (weight_sum > 0.0f) {
        result = weighted_sum / weight_sum;
        result = fminf(fmaxf(result, 0.0f), max_value);
    } else {
        result = self_value;
    }

    output[center_idx] = (T)result;
}

/* ----------------------------------------------------------------
 * Noise sigma estimation from image
 *
 * Uses adjacent pixel differences in the central region of
 * several sample slices.
 *
 * For i.i.d. noise with std sigma:
 *   E[|f(x+1) - f(x)|^2] = 2 * sigma^2  (+ signal gradient term)
 *
 * We sample from the central 50% of the image to avoid edges,
 * and use a sparse sampling to reduce the contribution of
 * structural gradients, giving a reasonable noise estimate.
 * ---------------------------------------------------------------- */
static float estimate_noise_sigma(const char *dir_path, char files[][MAX_PATH_LENGTH],
                                  int file_count, ImageInfo *info) {
    TIFF *tif;
    char full_path[MAX_PATH_LENGTH];
    tsize_t scanline_size;
    double sum_diff_sq = 0.0;
    long count = 0;

    int sample_count = file_count / 20;
    if (sample_count < 1) sample_count = 1;
    if (sample_count > 10) sample_count = 10;
    int sample_interval = file_count / sample_count;

    for (int i = 0; i < sample_count; i++) {
        int file_idx = i * sample_interval;
        snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, files[file_idx]);
        tif = TIFFOpen(full_path, "r");
        if (tif == NULL) continue;
        scanline_size = TIFFScanlineSize(tif);

        int start_y = info->height / 4;
        int end_y = 3 * info->height / 4;
        int y_step = (end_y - start_y) / 10;
        if (y_step < 1) y_step = 1;
        int start_x = info->width / 4;
        int end_x = 3 * info->width / 4;

        if (info->bits_per_sample == 8) {
            unsigned char *buf = (unsigned char*)_TIFFmalloc(scanline_size);
            if (buf == NULL) { TIFFClose(tif); continue; }
            for (int y = start_y; y < end_y; y += y_step) {
                if (TIFFReadScanline(tif, buf, y, 0) >= 0) {
                    for (int x = start_x; x < end_x - 1; x++) {
                        double diff = (double)buf[x+1] - (double)buf[x];
                        sum_diff_sq += diff * diff;
                        count++;
                    }
                }
            }
            _TIFFfree(buf);
        } else if (info->bits_per_sample == 16) {
            unsigned short *buf = (unsigned short*)_TIFFmalloc(scanline_size);
            if (buf == NULL) { TIFFClose(tif); continue; }
            for (int y = start_y; y < end_y; y += y_step) {
                if (TIFFReadScanline(tif, buf, y, 0) >= 0) {
                    for (int x = start_x; x < end_x - 1; x++) {
                        double diff = (double)buf[x+1] - (double)buf[x];
                        sum_diff_sq += diff * diff;
                        count++;
                    }
                }
            }
            _TIFFfree(buf);
        } else if (info->bits_per_sample == 32) {
            float *buf = (float*)_TIFFmalloc(scanline_size);
            if (buf == NULL) { TIFFClose(tif); continue; }
            for (int y = start_y; y < end_y; y += y_step) {
                if (TIFFReadScanline(tif, buf, y, 0) >= 0) {
                    for (int x = start_x; x < end_x - 1; x++) {
                        float val0 = buf[x], val1 = buf[x+1];
                        if (isfinite(val0) && isfinite(val1)) {
                            double diff = (double)val1 - (double)val0;
                            sum_diff_sq += diff * diff;
                            count++;
                        }
                    }
                }
            }
            _TIFFfree(buf);
        }
        TIFFClose(tif);
    }

    if (count < 2) return 1.0f;

    /* sigma = sqrt(mean_diff_sq / 2)
     * Adjacent pixel diffs include both noise and signal gradient.
     * For typical CT images the signal gradient contribution causes
     * this estimate to be somewhat larger than the true noise sigma,
     * which is desirable: it yields h ~ sigma_est * 1.2 in the
     * range that empirically produces good NLM results. */
    double mean_diff_sq = sum_diff_sq / (double)count;
    float sigma = (float)sqrt(mean_diff_sq / 2.0);

    return sigma;
}

/* ----------------------------------------------------------------
 * File I/O and utility functions
 * ---------------------------------------------------------------- */

static int create_directory(const char *path) {
#ifdef _WIN32
    if (GetFileAttributesA(path) == INVALID_FILE_ATTRIBUTES) {
        if (_mkdir(path) != 0) return -1;
    }
#else
    struct stat st;
    if (stat(path, &st) != 0) {
        if (mkdir(path, 0755) != 0) return -1;
    }
#endif
    return 0;
}

static int compare_strings(const void *a, const void *b) {
    return strcmp((const char *)a, (const char *)b);
}

static int get_tiff_files(const char *dir_path, char files[][MAX_PATH_LENGTH], int *file_count) {
#ifdef _WIN32
    WIN32_FIND_DATAA find_data;
    HANDLE find_handle;
    char search_path[MAX_PATH_LENGTH];
    snprintf(search_path, MAX_PATH_LENGTH, "%s\\*.tif*", dir_path);
    find_handle = FindFirstFileA(search_path, &find_data);
    if (find_handle == INVALID_HANDLE_VALUE) return -1;
    *file_count = 0;
    do {
        if (!(find_data.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
            strncpy(files[*file_count], find_data.cFileName, MAX_PATH_LENGTH - 1);
            files[*file_count][MAX_PATH_LENGTH - 1] = '\0';
            (*file_count)++;
            if (*file_count >= MAX_FILES) break;
        }
    } while (FindNextFileA(find_handle, &find_data));
    FindClose(find_handle);
#else
    DIR *dir;
    struct dirent *entry;
    char *ext;
    dir = opendir(dir_path);
    if (dir == NULL) return -1;
    *file_count = 0;
    while ((entry = readdir(dir)) != NULL) {
        ext = strrchr(entry->d_name, '.');
        if (ext != NULL && (strcmp(ext, ".tif") == 0 || strcmp(ext, ".tiff") == 0)) {
            strncpy(files[*file_count], entry->d_name, MAX_PATH_LENGTH - 1);
            files[*file_count][MAX_PATH_LENGTH - 1] = '\0';
            (*file_count)++;
            if (*file_count >= MAX_FILES) break;
        }
    }
    closedir(dir);
#endif
    qsort(files, *file_count, MAX_PATH_LENGTH, compare_strings);
    return 0;
}

static int get_image_info(const char *dir_path, const char *filename, ImageInfo *info) {
    TIFF *tif;
    char full_path[MAX_PATH_LENGTH];
    uint32_t width, height;
    uint16_t bits_per_sample, samples_per_pixel, sample_format;
    snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, filename);
    tif = TIFFOpen(full_path, "r");
    if (tif == NULL) return -1;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);
    sample_format = SAMPLEFORMAT_UINT;
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &sample_format);
    TIFFClose(tif);
    info->width = width; info->height = height;
    info->bits_per_sample = bits_per_sample; info->samples_per_pixel = samples_per_pixel;
    info->sample_format = sample_format;
    info->bytes_per_pixel = (bits_per_sample / 8) * samples_per_pixel;
    info->bytes_per_slice = (size_t)width * height * info->bytes_per_pixel;
    return 0;
}

static int calculate_overlap_size(FilterParams *params) {
    return params->search_radius + params->patch_radius + 1;
}

static int calculate_gpu_chunk_size(ImageInfo *info, int overlap_size) {
    size_t free_mem, total_mem;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
    size_t available_mb = (size_t)(free_mem * GPU_MEMORY_FRACTION) / (1024 * 1024);
    size_t slice_mb = info->bytes_per_slice / (1024 * 1024);
    if (slice_mb < 1) slice_mb = 1;
    int chunk_depth = (int)(available_mb / (2 * slice_mb));
    int min_chunk_depth = 4 * overlap_size;
    if (chunk_depth < min_chunk_depth) chunk_depth = min_chunk_depth;
    return chunk_depth;
}

static ChunkData* allocate_chunk(ImageInfo *info, int chunk_depth) {
    ChunkData *chunk = (ChunkData*)malloc(sizeof(ChunkData));
    if (chunk == NULL) return NULL;
    chunk->chunk_depth = chunk_depth;
    size_t chunk_size = (size_t)chunk_depth * info->bytes_per_slice;
    CUDA_CHECK(cudaMallocHost(&chunk->h_data, chunk_size));
    CUDA_CHECK(cudaMalloc(&chunk->d_data, chunk_size));
    CUDA_CHECK(cudaMalloc(&chunk->d_output, chunk_size));
    return chunk;
}

static void free_chunk(ChunkData *chunk) {
    if (chunk != NULL) {
        if (chunk->h_data != NULL) cudaFreeHost(chunk->h_data);
        if (chunk->d_data != NULL) cudaFree(chunk->d_data);
        if (chunk->d_output != NULL) cudaFree(chunk->d_output);
        free(chunk);
    }
}

static int load_chunk(const char *dir_path, char files[][MAX_PATH_LENGTH],
                      ImageInfo *info, ChunkData *chunk) {
    TIFF *tif;
    char full_path[MAX_PATH_LENGTH];
    int z, y;
    tsize_t scanline_size;
    unsigned char *buffer;
    void *slice_data;
    for (z = 0; z < chunk->chunk_depth && (chunk->start_z + z) < (int)info->depth; z++) {
        snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR,
                 files[chunk->start_z + z]);
        tif = TIFFOpen(full_path, "r");
        if (tif == NULL) return -1;
        scanline_size = TIFFScanlineSize(tif);
        buffer = (unsigned char*)_TIFFmalloc(scanline_size);
        if (buffer == NULL) { TIFFClose(tif); return -1; }
        slice_data = (char*)chunk->h_data + z * info->bytes_per_slice;
        for (y = 0; y < (int)info->height; y++) {
            if (TIFFReadScanline(tif, buffer, y, 0) < 0) {
                _TIFFfree(buffer); TIFFClose(tif); return -1;
            }
            memcpy((char*)slice_data + y * scanline_size, buffer, scanline_size);
        }
        _TIFFfree(buffer);
        TIFFClose(tif);
    }
    chunk->end_z = chunk->start_z + z - 1;
    return 0;
}

/* Copy metadata tags from input TIFF (in) to output TIFF (out).
 * Tags that affect pixel data layout (dimensions, bps, compression, photometric,
 * planar config, strip/tile geometry) are intentionally NOT copied here, because
 * the caller sets those explicitly to match the data actually being written.
 * Call this BEFORE closing the input TIFF (string pointers belong to the input). */
static void copy_tiff_metadata(TIFF *in, TIFF *out) {
    uint32_t u32;
    uint16_t u16, u16a, u16b, *u16ptr;
    float f;
    double dbl;
    char *str;
    uint32_t count;
    void *data;

    /* ASCII / string tags */
    if (TIFFGetField(in, TIFFTAG_IMAGEDESCRIPTION, &str)) TIFFSetField(out, TIFFTAG_IMAGEDESCRIPTION, str);
    if (TIFFGetField(in, TIFFTAG_MAKE,             &str)) TIFFSetField(out, TIFFTAG_MAKE,             str);
    if (TIFFGetField(in, TIFFTAG_MODEL,            &str)) TIFFSetField(out, TIFFTAG_MODEL,            str);
    if (TIFFGetField(in, TIFFTAG_SOFTWARE,         &str)) TIFFSetField(out, TIFFTAG_SOFTWARE,         str);
    if (TIFFGetField(in, TIFFTAG_DATETIME,         &str)) TIFFSetField(out, TIFFTAG_DATETIME,         str);
    if (TIFFGetField(in, TIFFTAG_ARTIST,           &str)) TIFFSetField(out, TIFFTAG_ARTIST,           str);
    if (TIFFGetField(in, TIFFTAG_HOSTCOMPUTER,     &str)) TIFFSetField(out, TIFFTAG_HOSTCOMPUTER,     str);
    if (TIFFGetField(in, TIFFTAG_COPYRIGHT,        &str)) TIFFSetField(out, TIFFTAG_COPYRIGHT,        str);
    if (TIFFGetField(in, TIFFTAG_DOCUMENTNAME,     &str)) TIFFSetField(out, TIFFTAG_DOCUMENTNAME,     str);
    if (TIFFGetField(in, TIFFTAG_PAGENAME,         &str)) TIFFSetField(out, TIFFTAG_PAGENAME,         str);
    if (TIFFGetField(in, TIFFTAG_TARGETPRINTER,    &str)) TIFFSetField(out, TIFFTAG_TARGETPRINTER,    str);

    /* uint32_t tags */
    if (TIFFGetField(in, TIFFTAG_SUBFILETYPE, &u32)) TIFFSetField(out, TIFFTAG_SUBFILETYPE, u32);

    /* uint16_t tags */
    if (TIFFGetField(in, TIFFTAG_ORIENTATION,      &u16)) TIFFSetField(out, TIFFTAG_ORIENTATION,      u16);
    if (TIFFGetField(in, TIFFTAG_RESOLUTIONUNIT,   &u16)) TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT,   u16);
    if (TIFFGetField(in, TIFFTAG_MINSAMPLEVALUE,   &u16)) TIFFSetField(out, TIFFTAG_MINSAMPLEVALUE,   u16);
    if (TIFFGetField(in, TIFFTAG_MAXSAMPLEVALUE,   &u16)) TIFFSetField(out, TIFFTAG_MAXSAMPLEVALUE,   u16);

    /* float tags */
    if (TIFFGetField(in, TIFFTAG_XRESOLUTION, &f)) TIFFSetField(out, TIFFTAG_XRESOLUTION, f);
    if (TIFFGetField(in, TIFFTAG_YRESOLUTION, &f)) TIFFSetField(out, TIFFTAG_YRESOLUTION, f);
    if (TIFFGetField(in, TIFFTAG_XPOSITION,   &f)) TIFFSetField(out, TIFFTAG_XPOSITION,   f);
    if (TIFFGetField(in, TIFFTAG_YPOSITION,   &f)) TIFFSetField(out, TIFFTAG_YPOSITION,   f);

    /* double tags (libtiff 4.x: SMin/SMaxSampleValue are double) */
    if (TIFFGetField(in, TIFFTAG_SMINSAMPLEVALUE, &dbl)) TIFFSetField(out, TIFFTAG_SMINSAMPLEVALUE, dbl);
    if (TIFFGetField(in, TIFFTAG_SMAXSAMPLEVALUE, &dbl)) TIFFSetField(out, TIFFTAG_SMAXSAMPLEVALUE, dbl);
    if (TIFFGetField(in, TIFFTAG_STONITS,         &dbl)) TIFFSetField(out, TIFFTAG_STONITS,         dbl);

    /* PageNumber: two uint16_t values */
    if (TIFFGetField(in, TIFFTAG_PAGENUMBER, &u16a, &u16b)) {
        TIFFSetField(out, TIFFTAG_PAGENUMBER, u16a, u16b);
    }

    /* ICC profile: uint32_t count + data pointer */
    if (TIFFGetField(in, TIFFTAG_ICCPROFILE, &count, &data)) {
        TIFFSetField(out, TIFFTAG_ICCPROFILE, count, data);
    }

    /* ExtraSamples: uint16_t count + uint16_t array */
    if (TIFFGetField(in, TIFFTAG_EXTRASAMPLES, &u16, &u16ptr)) {
        TIFFSetField(out, TIFFTAG_EXTRASAMPLES, u16, u16ptr);
    }
}

static int save_chunk_valid_region(const char *input_dir, const char *dir_path,
                                   char files[][MAX_PATH_LENGTH],
                                   ImageInfo *info, ChunkData *chunk) {
    TIFF *tif, *tif_in;
    char full_path[MAX_PATH_LENGTH];
    char in_path[MAX_PATH_LENGTH];
    int z, y;
    tsize_t scanline_size;
    void *slice_data;
    int global_z;
    for (z = chunk->valid_start; z <= chunk->valid_end; z++) {
        global_z = chunk->start_z + z;
        if (global_z >= (int)info->depth) break;
        snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, files[global_z]);
        snprintf(in_path,   MAX_PATH_LENGTH, "%s%s%s", input_dir, PATH_SEPARATOR, files[global_z]);
        tif = TIFFOpen(full_path, "w");
        if (tif == NULL) return -1;

        /* Copy metadata tags from corresponding input file first. */
        tif_in = TIFFOpen(in_path, "r");
        if (tif_in != NULL) {
            copy_tiff_metadata(tif_in, tif);
            TIFFClose(tif_in);
        }

        /* Set / override tags that must match the data being written. */
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, info->width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, info->height);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, info->bits_per_sample);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, info->samples_per_pixel);
        if (info->sample_format == SAMPLEFORMAT_UINT ||
            info->sample_format == SAMPLEFORMAT_INT ||
            info->sample_format == SAMPLEFORMAT_IEEEFP ||
            info->sample_format == SAMPLEFORMAT_VOID ||
            info->sample_format == SAMPLEFORMAT_COMPLEXINT ||
            info->sample_format == SAMPLEFORMAT_COMPLEXIEEEFP) {
            TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, info->sample_format);
        } else {
            if (info->bits_per_sample == 32) {
                TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
            } else {
                TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
            }
        }
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        slice_data = (char*)chunk->h_data + z * info->bytes_per_slice;
        scanline_size = info->width * info->bytes_per_pixel;
        for (y = 0; y < (int)info->height; y++) {
            if (TIFFWriteScanline(tif, (char*)slice_data + y * scanline_size, y, 0) < 0) {
                TIFFClose(tif); return -1;
            }
        }
        TIFFClose(tif);
    }
    return 0;
}

static void process_chunk_on_gpu(ChunkData *chunk, ImageInfo *info, FilterParams *params) {
    size_t chunk_bytes = (size_t)chunk->chunk_depth * info->bytes_per_slice;
    int valid_depth = chunk->valid_end - chunk->valid_start + 1;

    CUDA_CHECK(cudaMemcpy(chunk->d_data, chunk->h_data, chunk_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(chunk->d_output, chunk->d_data, chunk_bytes, cudaMemcpyDeviceToDevice));

    float h_sq_inv = 1.0f / (params->h * params->h);

    float max_value;
    if (info->bits_per_sample == 8) max_value = 255.0f;
    else if (info->bits_per_sample == 16) max_value = 65535.0f;
    else max_value = FLT_MAX;

    dim3 block(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
    dim3 grid(
        (info->width + block.x - 1) / block.x,
        (info->height + block.y - 1) / block.y,
        (valid_depth + block.z - 1) / block.z
    );

    if (info->bits_per_sample == 8) {
        nlm_filter_3d_kernel<unsigned char><<<grid, block>>>(
            (unsigned char*)chunk->d_data, (unsigned char*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            chunk->valid_start, chunk->valid_end,
            params->patch_radius, params->search_radius,
            h_sq_inv, max_value);
    } else if (info->bits_per_sample == 16) {
        nlm_filter_3d_kernel<unsigned short><<<grid, block>>>(
            (unsigned short*)chunk->d_data, (unsigned short*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            chunk->valid_start, chunk->valid_end,
            params->patch_radius, params->search_radius,
            h_sq_inv, max_value);
    } else if (info->bits_per_sample == 32 && info->sample_format == SAMPLEFORMAT_IEEEFP) {
        nlm_filter_3d_kernel<float><<<grid, block>>>(
            (float*)chunk->d_data, (float*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            chunk->valid_start, chunk->valid_end,
            params->patch_radius, params->search_radius,
            h_sq_inv, max_value);
    }

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(chunk->h_data, chunk->d_output, chunk_bytes, cudaMemcpyDeviceToHost));
}

/* ----------------------------------------------------------------
 * Main
 * ---------------------------------------------------------------- */
int main(int argc, char *argv[]) {
    char input_dir[MAX_PATH_LENGTH];
    char output_dir[MAX_PATH_LENGTH];
    FilterParams params;
    char (*files)[MAX_PATH_LENGTH];
    int file_count;
    ImageInfo info;
    ChunkData *chunk;
    int chunk_size, num_chunks, chunk_idx;
    int overlap_size, valid_chunk_size;

    /* GPU setup */
    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count == 0) { fprintf(stderr, "Error: No CUDA-capable devices found\n"); return 1; }
    int best_device = 0, max_sm = 0;
    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, i));
        if (prop.multiProcessorCount > max_sm) { max_sm = prop.multiProcessorCount; best_device = i; }
    }
    CUDA_CHECK(cudaSetDevice(best_device));

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input_dir> <output_dir> [patch_radius] [search_radius] [h]\n", argv[0]);
        fprintf(stderr, "  patch_radius:  Half-size of comparison patch (default: 1 -> 3x3x3)\n");
        fprintf(stderr, "  search_radius: Half-size of search window   (default: 3 -> 7x7x7)\n");
        fprintf(stderr, "  h:             Filtering strength            (default: auto)\n");
        fprintf(stderr, "\n");
        fprintf(stderr, "  Smaller h = sharper (more noise remains)\n");
        fprintf(stderr, "  Larger  h = smoother (more blurring)\n");
        fprintf(stderr, "  Auto h is estimated from image noise level.\n");
        return 1;
    }

    strncpy(input_dir, argv[1], MAX_PATH_LENGTH - 1); input_dir[MAX_PATH_LENGTH - 1] = '\0';
    strncpy(output_dir, argv[2], MAX_PATH_LENGTH - 1); output_dir[MAX_PATH_LENGTH - 1] = '\0';

    /* Defaults: patch_radius=1 (3x3x3), search_radius=3 (7x7x7), h=auto */
    params.patch_radius = 1;
    params.search_radius = 3;
    params.h = -1.0f;

    if (argc > 3) params.patch_radius = atoi(argv[3]);
    if (argc > 4) params.search_radius = atoi(argv[4]);
    if (argc > 5) params.h = (float)atof(argv[5]);

    if (params.patch_radius < 1 || params.patch_radius > 5) {
        fprintf(stderr, "Error: patch_radius must be between 1 and 5\n"); return 1;
    }
    if (params.search_radius < 1 || params.search_radius > 15) {
        fprintf(stderr, "Error: search_radius must be between 1 and 15\n"); return 1;
    }

    if (create_directory(output_dir) != 0) {
        fprintf(stderr, "Error: Cannot create output directory\n"); return 1;
    }

    files = (char (*)[MAX_PATH_LENGTH])malloc((size_t)MAX_FILES * MAX_PATH_LENGTH);
    if (files == NULL) { fprintf(stderr, "Error: Memory allocation failed\n"); return 1; }
    if (get_tiff_files(input_dir, files, &file_count) != 0 || file_count == 0) {
        fprintf(stderr, "Error: No TIFF files found\n"); free(files); return 1;
    }
    if (get_image_info(input_dir, files[0], &info) != 0) {
        fprintf(stderr, "Error: Cannot read image information\n"); free(files); return 1;
    }
    info.depth = file_count;

    /* Auto-estimate h from noise sigma if not specified */
    if (params.h <= 0.0f) {
        float sigma_noise = estimate_noise_sigma(input_dir, files, file_count, &info);
        /* h = 1.2 * sigma_noise
         *
         * The estimate_noise_sigma function measures the RMS of adjacent
         * pixel differences, which includes both noise and signal gradient.
         * For typical 16-bit CT images from BL47XU (sigma ~150-250),
         * h = 1.2 * sigma gives results in the range of ~180-300,
         * which empirically produces good denoising with edge preservation.
         *
         * Adjustment guide:
         *   h too small -> noisy result (insufficient denoising)
         *   h too large -> blurry result (over-smoothing)
         *   Good starting point: h ~ sigma_noise to 1.5*sigma_noise
         */
        params.h = 1.2f * sigma_noise;
        if (params.h <= 0.0f) params.h = 1.0f;
        fprintf(stderr, "  Estimated noise sigma: %.4f\n", sigma_noise);
    }

    fprintf(stderr, "NLM 3D Filter (GPU)\n");
    fprintf(stderr, "  Image: %u x %u x %u, %u-bit\n",
            info.width, info.height, info.depth, info.bits_per_sample);
    fprintf(stderr, "  Patch radius:  %d (patch size: %dx%dx%d)\n",
            params.patch_radius,
            2*params.patch_radius+1, 2*params.patch_radius+1, 2*params.patch_radius+1);
    fprintf(stderr, "  Search radius: %d (search size: %dx%dx%d)\n",
            params.search_radius,
            2*params.search_radius+1, 2*params.search_radius+1, 2*params.search_radius+1);
    fprintf(stderr, "  h: %.4f\n", params.h);

    overlap_size = calculate_overlap_size(&params);
    chunk_size = calculate_gpu_chunk_size(&info, overlap_size);
    valid_chunk_size = chunk_size - 2 * overlap_size;
    if (valid_chunk_size < 1) { valid_chunk_size = 1; chunk_size = 1 + 2 * overlap_size; }
    num_chunks = (file_count + valid_chunk_size - 1) / valid_chunk_size;

    fprintf(stderr, "  Overlap: %d slices, Chunk size: %d slices, Chunks: %d\n",
            overlap_size, chunk_size, num_chunks);

    chunk = allocate_chunk(&info, chunk_size);
    if (chunk == NULL) { fprintf(stderr, "Error: Memory allocation failed\n"); free(files); return 1; }

    for (chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {
        int vgs = chunk_idx * valid_chunk_size;
        int vge = vgs + valid_chunk_size - 1;
        if (vge >= file_count) vge = file_count - 1;
        int cgs = vgs - overlap_size; if (cgs < 0) cgs = 0;
        int cge = vge + overlap_size; if (cge >= file_count) cge = file_count - 1;

        chunk->start_z = cgs; chunk->end_z = cge;
        chunk->chunk_depth = cge - cgs + 1;
        chunk->valid_start = vgs - cgs; chunk->valid_end = vge - cgs;

        fprintf(stderr, "  Chunk %d/%d: z=%d..%d (valid %d..%d)\n",
                chunk_idx + 1, num_chunks, cgs, cge, vgs, vge);

        if (load_chunk(input_dir, files, &info, chunk) != 0) {
            fprintf(stderr, "Error: Failed to load chunk %d\n", chunk_idx);
            free_chunk(chunk); free(files); return 1;
        }
        process_chunk_on_gpu(chunk, &info, &params);
        if (save_chunk_valid_region(input_dir, output_dir, files, &info, chunk) != 0) {
            fprintf(stderr, "Error: Failed to save chunk %d\n", chunk_idx);
            free_chunk(chunk); free(files); return 1;
        }
    }

    free_chunk(chunk); free(files);

    /* append to log file */
    {
        FILE *f;
        int i;
        if ((f = fopen("cmd-hst.log", "a")) != NULL) {
            for (i = 0; i < argc; ++i) fprintf(f, "%s ", argv[i]);
            fprintf(f, "\n");
            fprintf(f, "   %% patch_radius %d  search_radius %d  h %g\n", params.patch_radius, params.search_radius, params.h);
            fclose(f);
        }
    }

    CUDA_CHECK(cudaDeviceReset());
    return 0;
}