/*
 * tif_adf_g.cu - 3D Anisotropic Diffusion Filter for TIFF Image Stacks (CUDA)
 *
 * Algorithm: Perona-Malik anisotropic diffusion (3D extension)
 *
 *   du/dt = div( g(|nabla u|) * nabla u )
 *   g(s) = 1 / (1 + (s/K)^2)        (Perona-Malik type 1)
 *   or
 *   g(s) = exp(-(s/K)^2)             (Perona-Malik type 2)
 *
 * Compile:
 *   Windows: nvcc -O3 -I%CUDAINCL% -o tif_adf_g.exe tif_adf_g.cu -use_fast_math -Xcompiler "/wd 4819" libtiff.lib
 *   Linux:   nvcc -O3 -o tif_adf_g tif_adf_g.cu -use_fast_math -ltiff
 *
 * Usage:
 *   tif_adf_g <input_dir> <output_dir> [iterations] [K] [dt] [mode]
 *
 *   iterations: Number of diffusion steps (default: 20)
 *   K:          Edge threshold (default: auto from noise sigma)
 *   dt:         Time step (default: 0.14, must be < 1/6 for 3D stability)
 *   mode:       1 = g(s)=1/(1+(s/K)^2),  2 = g(s)=exp(-(s/K)^2) (default: 1)
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
#define BLOCK_SIZE_Z 8
#define GPU_MEMORY_FRACTION 0.7f
#define AD_OVERLAP_MARGIN 2

typedef struct {
    int iterations;
    float K;
    float dt;
    int mode;  /* 1 or 2 */
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
 * Perona-Malik 3D diffusion kernel
 *
 * For each voxel, compute 6-connected finite differences,
 * apply edge-stopping function g(), and update.
 * ---------------------------------------------------------------- */
__global__ void aniso_diffusion_kernel(
    const float* __restrict__ u,
    float* __restrict__ u_new,
    int width, int height, int depth,
    float K_sq_inv, float dt, int mode)
{
    ENABLE_SMEM_SPILLING();
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x >= width || y >= height || z >= depth) return;

    size_t slice = (size_t)width * height;
    size_t idx = (size_t)z * slice + (size_t)y * width + x;
    float uc = u[idx];

    /* 6-connected differences */
    float dE = (x + 1 < width)  ? u[idx + 1]     - uc : 0.0f;
    float dW = (x - 1 >= 0)     ? u[idx - 1]     - uc : 0.0f;
    float dS = (y + 1 < height) ? u[idx + width]  - uc : 0.0f;
    float dN = (y - 1 >= 0)     ? u[idx - width]  - uc : 0.0f;
    float dU = (z + 1 < depth)  ? u[idx + slice]  - uc : 0.0f;
    float dD = (z - 1 >= 0)     ? u[idx - slice]  - uc : 0.0f;

    float gE, gW, gS, gN, gU, gD;

    if (mode == 1) {
        /* g(s) = 1 / (1 + (s/K)^2) */
        gE = 1.0f / (1.0f + dE * dE * K_sq_inv);
        gW = 1.0f / (1.0f + dW * dW * K_sq_inv);
        gS = 1.0f / (1.0f + dS * dS * K_sq_inv);
        gN = 1.0f / (1.0f + dN * dN * K_sq_inv);
        gU = 1.0f / (1.0f + dU * dU * K_sq_inv);
        gD = 1.0f / (1.0f + dD * dD * K_sq_inv);
    } else {
        /* g(s) = exp(-(s/K)^2) */
        gE = __expf(-dE * dE * K_sq_inv);
        gW = __expf(-dW * dW * K_sq_inv);
        gS = __expf(-dS * dS * K_sq_inv);
        gN = __expf(-dN * dN * K_sq_inv);
        gU = __expf(-dU * dU * K_sq_inv);
        gD = __expf(-dD * dD * K_sq_inv);
    }

    float update = dt * (gE * dE + gW * dW + gS * dS + gN * dN + gU * dU + gD * dD);
    u_new[idx] = uc + update;
}

template<typename T>
__global__ void to_float_kernel(const T* __restrict__ input, float* __restrict__ output,
                                int width, int height, int depth) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    if (x >= width || y >= height || z >= depth) return;
    size_t idx = (size_t)z * width * height + (size_t)y * width + x;
    output[idx] = (float)input[idx];
}

template<typename T>
__global__ void from_float_kernel(const float* __restrict__ input, T* __restrict__ output,
                                  int width, int height, int depth,
                                  float min_value, float max_value) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    int z = blockIdx.z * blockDim.z + threadIdx.z;
    if (x >= width || y >= height || z >= depth) return;
    size_t idx = (size_t)z * width * height + (size_t)y * width + x;
    float val = fminf(fmaxf(input[idx], min_value), max_value);
    output[idx] = (T)val;
}

/* ----------------------------------------------------------------
 * Noise sigma estimation (same as NLM version)
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
        int start_y = info->height / 4, end_y = 3 * info->height / 4;
        int y_step = (end_y - start_y) / 10;
        if (y_step < 1) y_step = 1;
        int start_x = info->width / 4, end_x = 3 * info->width / 4;

        if (info->bits_per_sample == 8) {
            unsigned char *buf = (unsigned char*)_TIFFmalloc(scanline_size);
            if (buf == NULL) { TIFFClose(tif); continue; }
            for (int y = start_y; y < end_y; y += y_step) {
                if (TIFFReadScanline(tif, buf, y, 0) >= 0) {
                    for (int x = start_x; x < end_x - 1; x++) {
                        double diff = (double)buf[x+1] - (double)buf[x];
                        sum_diff_sq += diff * diff; count++;
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
                        sum_diff_sq += diff * diff; count++;
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
                        if (isfinite(buf[x]) && isfinite(buf[x+1])) {
                            double diff = (double)buf[x+1] - (double)buf[x];
                            sum_diff_sq += diff * diff; count++;
                        }
                    }
                }
            }
            _TIFFfree(buf);
        }
        TIFFClose(tif);
    }
    if (count < 2) return 1.0f;
    return (float)sqrt(sum_diff_sq / (double)count / 2.0);
}

/* ----------------------------------------------------------------
 * File I/O (same pattern as other filters)
 * ---------------------------------------------------------------- */
static int create_directory(const char *path) {
#ifdef _WIN32
    if (GetFileAttributesA(path) == INVALID_FILE_ATTRIBUTES) { if (_mkdir(path) != 0) return -1; }
#else
    struct stat st;
    if (stat(path, &st) != 0) { if (mkdir(path, 0755) != 0) return -1; }
#endif
    return 0;
}

static int compare_strings(const void *a, const void *b) { return strcmp((const char *)a, (const char *)b); }

static int get_tiff_files(const char *dir_path, char files[][MAX_PATH_LENGTH], int *file_count) {
#ifdef _WIN32
    WIN32_FIND_DATAA fd; HANDLE fh;
    char sp[MAX_PATH_LENGTH];
    snprintf(sp, MAX_PATH_LENGTH, "%s\\*.tif*", dir_path);
    fh = FindFirstFileA(sp, &fd);
    if (fh == INVALID_HANDLE_VALUE) return -1;
    *file_count = 0;
    do { if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
        strncpy(files[*file_count], fd.cFileName, MAX_PATH_LENGTH - 1);
        files[*file_count][MAX_PATH_LENGTH - 1] = '\0'; (*file_count)++;
        if (*file_count >= MAX_FILES) break;
    }} while (FindNextFileA(fh, &fd));
    FindClose(fh);
#else
    DIR *dir; struct dirent *entry; char *ext;
    dir = opendir(dir_path); if (dir == NULL) return -1;
    *file_count = 0;
    while ((entry = readdir(dir)) != NULL) {
        ext = strrchr(entry->d_name, '.');
        if (ext != NULL && (strcmp(ext, ".tif") == 0 || strcmp(ext, ".tiff") == 0)) {
            strncpy(files[*file_count], entry->d_name, MAX_PATH_LENGTH - 1);
            files[*file_count][MAX_PATH_LENGTH - 1] = '\0'; (*file_count)++;
            if (*file_count >= MAX_FILES) break;
        }
    }
    closedir(dir);
#endif
    qsort(files, *file_count, MAX_PATH_LENGTH, compare_strings);
    return 0;
}

static int get_image_info(const char *dir_path, const char *filename, ImageInfo *info) {
    TIFF *tif; char full_path[MAX_PATH_LENGTH];
    uint32_t width, height; uint16_t bps, spp, sf;
    snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, filename);
    tif = TIFFOpen(full_path, "r"); if (tif == NULL) return -1;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    sf = SAMPLEFORMAT_UINT; TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &sf);
    TIFFClose(tif);
    info->width = width; info->height = height;
    info->bits_per_sample = bps; info->samples_per_pixel = spp;
    info->sample_format = sf;
    info->bytes_per_pixel = (bps / 8) * spp;
    info->bytes_per_slice = (size_t)width * height * info->bytes_per_pixel;
    return 0;
}

static int calculate_overlap_size(void) { return AD_OVERLAP_MARGIN + 1; }

static int calculate_gpu_chunk_size(ImageInfo *info, int overlap_size) {
    size_t free_mem, total_mem;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
    size_t available_mb = (size_t)(free_mem * GPU_MEMORY_FRACTION) / (1024 * 1024);
    size_t fslice_mb = ((size_t)info->width * info->height * sizeof(float)) / (1024 * 1024);
    if (fslice_mb < 1) fslice_mb = 1;
    /* Need: d_data + d_output + d_u + d_u_new = 4 float buffers */
    int chunk_depth = (int)(available_mb / (5 * fslice_mb));
    int min_cd = 4 * overlap_size;
    if (chunk_depth < min_cd) chunk_depth = min_cd;
    return chunk_depth;
}

static ChunkData* allocate_chunk(ImageInfo *info, int chunk_depth) {
    ChunkData *c = (ChunkData*)malloc(sizeof(ChunkData)); if (!c) return NULL;
    c->chunk_depth = chunk_depth;
    size_t cs = (size_t)chunk_depth * info->bytes_per_slice;
    CUDA_CHECK(cudaMallocHost(&c->h_data, cs));
    CUDA_CHECK(cudaMalloc(&c->d_data, cs));
    CUDA_CHECK(cudaMalloc(&c->d_output, cs));
    return c;
}

static void free_chunk(ChunkData *c) {
    if (c) { if (c->h_data) cudaFreeHost(c->h_data);
    if (c->d_data) cudaFree(c->d_data); if (c->d_output) cudaFree(c->d_output); free(c); }
}

static int load_chunk(const char *dir_path, char files[][MAX_PATH_LENGTH],
                      ImageInfo *info, ChunkData *chunk) {
    TIFF *tif; char fp[MAX_PATH_LENGTH]; int z, y; tsize_t ss; unsigned char *buf; void *sd;
    for (z = 0; z < chunk->chunk_depth && (chunk->start_z + z) < (int)info->depth; z++) {
        snprintf(fp, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, files[chunk->start_z + z]);
        tif = TIFFOpen(fp, "r"); if (!tif) return -1;
        ss = TIFFScanlineSize(tif); buf = (unsigned char*)_TIFFmalloc(ss);
        if (!buf) { TIFFClose(tif); return -1; }
        sd = (char*)chunk->h_data + z * info->bytes_per_slice;
        for (y = 0; y < (int)info->height; y++) {
            if (TIFFReadScanline(tif, buf, y, 0) < 0) { _TIFFfree(buf); TIFFClose(tif); return -1; }
            memcpy((char*)sd + y * ss, buf, ss);
        }
        _TIFFfree(buf); TIFFClose(tif);
    }
    chunk->end_z = chunk->start_z + z - 1; return 0;
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
    TIFF *tif, *tif_in; char fp[MAX_PATH_LENGTH]; char in_path[MAX_PATH_LENGTH];
    int z, y; tsize_t ss; void *sd; int gz;
    for (z = chunk->valid_start; z <= chunk->valid_end; z++) {
        gz = chunk->start_z + z; if (gz >= (int)info->depth) break;
        snprintf(fp, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, files[gz]);
        snprintf(in_path, MAX_PATH_LENGTH, "%s%s%s", input_dir, PATH_SEPARATOR, files[gz]);
        tif = TIFFOpen(fp, "w"); if (!tif) return -1;
        tif_in = TIFFOpen(in_path, "r");
        if (tif_in != NULL) { copy_tiff_metadata(tif_in, tif); TIFFClose(tif_in); }
        TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, info->width);
        TIFFSetField(tif, TIFFTAG_IMAGELENGTH, info->height);
        TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, info->bits_per_sample);
        TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, info->samples_per_pixel);
        if (info->sample_format == SAMPLEFORMAT_UINT || info->sample_format == SAMPLEFORMAT_INT ||
            info->sample_format == SAMPLEFORMAT_IEEEFP || info->sample_format == SAMPLEFORMAT_VOID ||
            info->sample_format == SAMPLEFORMAT_COMPLEXINT || info->sample_format == SAMPLEFORMAT_COMPLEXIEEEFP) {
            TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, info->sample_format);
        } else {
            TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, (info->bits_per_sample==32)?SAMPLEFORMAT_IEEEFP:SAMPLEFORMAT_UINT);
        }
        TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
        TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
        sd = (char*)chunk->h_data + z * info->bytes_per_slice;
        ss = info->width * info->bytes_per_pixel;
        for (y = 0; y < (int)info->height; y++) {
            if (TIFFWriteScanline(tif, (char*)sd + y * ss, y, 0) < 0) { TIFFClose(tif); return -1; }
        }
        TIFFClose(tif);
    }
    return 0;
}

static void process_chunk_on_gpu(ChunkData *chunk, ImageInfo *info, FilterParams *params) {
    size_t chunk_bytes = (size_t)chunk->chunk_depth * info->bytes_per_slice;
    int w = info->width, h = info->height, d = chunk->chunk_depth;
    size_t nv = (size_t)w * h * d;
    size_t fb = nv * sizeof(float);

    CUDA_CHECK(cudaMemcpy(chunk->d_data, chunk->h_data, chunk_bytes, cudaMemcpyHostToDevice));

    float *d_u, *d_u_new;
    CUDA_CHECK(cudaMalloc(&d_u,     fb));
    CUDA_CHECK(cudaMalloc(&d_u_new, fb));

    dim3 block(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
    dim3 grid((w+block.x-1)/block.x, (h+block.y-1)/block.y, (d+block.z-1)/block.z);

    if (info->bits_per_sample == 8) {
        to_float_kernel<unsigned char><<<grid, block>>>((unsigned char*)chunk->d_data, d_u, w, h, d);
    } else if (info->bits_per_sample == 16) {
        to_float_kernel<unsigned short><<<grid, block>>>((unsigned short*)chunk->d_data, d_u, w, h, d);
    } else {
        to_float_kernel<float><<<grid, block>>>((float*)chunk->d_data, d_u, w, h, d);
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    float K_sq_inv = 1.0f / (params->K * params->K);

    for (int iter = 0; iter < params->iterations; iter++) {
        aniso_diffusion_kernel<<<grid, block>>>(d_u, d_u_new, w, h, d, K_sq_inv, params->dt, params->mode);
        CUDA_CHECK(cudaDeviceSynchronize());
        /* Swap: u <- u_new */
        float *tmp = d_u; d_u = d_u_new; d_u_new = tmp;
    }

    float min_value, max_value;
    if (info->bits_per_sample == 8) {
        min_value = 0.0f;     max_value = 255.0f;
    } else if (info->bits_per_sample == 16) {
        min_value = 0.0f;     max_value = 65535.0f;
    } else {
        /* 32-bit float: allow negative values. */
        min_value = -FLT_MAX; max_value =  FLT_MAX;
    }

    if (info->bits_per_sample == 8) {
        from_float_kernel<unsigned char><<<grid, block>>>(d_u, (unsigned char*)chunk->d_output, w, h, d, min_value, max_value);
    } else if (info->bits_per_sample == 16) {
        from_float_kernel<unsigned short><<<grid, block>>>(d_u, (unsigned short*)chunk->d_output, w, h, d, min_value, max_value);
    } else {
        from_float_kernel<float><<<grid, block>>>(d_u, (float*)chunk->d_output, w, h, d, min_value, max_value);
    }
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(chunk->h_data, chunk->d_output, chunk_bytes, cudaMemcpyDeviceToHost));

    cudaFree(d_u); cudaFree(d_u_new);
}

int main(int argc, char *argv[]) {
    char input_dir[MAX_PATH_LENGTH], output_dir[MAX_PATH_LENGTH];
    FilterParams params;
    char (*files)[MAX_PATH_LENGTH]; int file_count;
    ImageInfo info; ChunkData *chunk;
    int chunk_size, num_chunks, chunk_idx, overlap_size, valid_chunk_size;

    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count == 0) { fprintf(stderr, "Error: No CUDA devices\n"); return 1; }
    int best_device = 0, max_sm = 0;
    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp prop; CUDA_CHECK(cudaGetDeviceProperties(&prop, i));
        if (prop.multiProcessorCount > max_sm) { max_sm = prop.multiProcessorCount; best_device = i; }
    }
    CUDA_CHECK(cudaSetDevice(best_device));

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input_dir> <output_dir> [iterations] [K] [dt] [mode]\n", argv[0]);
        fprintf(stderr, "  iterations: Diffusion steps (default: 20)\n");
        fprintf(stderr, "  K:          Edge threshold  (default: auto from noise)\n");
        fprintf(stderr, "  dt:         Time step       (default: 0.14, <1/6 for stability)\n");
        fprintf(stderr, "  mode:       1=1/(1+(s/K)^2)  2=exp(-(s/K)^2) (default: 1)\n");
        return 1;
    }

    strncpy(input_dir, argv[1], MAX_PATH_LENGTH-1); input_dir[MAX_PATH_LENGTH-1]='\0';
    strncpy(output_dir, argv[2], MAX_PATH_LENGTH-1); output_dir[MAX_PATH_LENGTH-1]='\0';

    params.iterations = 20; params.K = -1.0f; params.dt = 0.14f; params.mode = 1;
    if (argc > 3) params.iterations = atoi(argv[3]);
    if (argc > 4) params.K = (float)atof(argv[4]);
    if (argc > 5) params.dt = (float)atof(argv[5]);
    if (argc > 6) params.mode = atoi(argv[6]);

    if (params.iterations < 1 || params.iterations > 10000) {
        fprintf(stderr, "Error: iterations must be 1-10000\n"); return 1;
    }
    if (params.dt <= 0.0f || params.dt >= 1.0f/6.0f) {
        fprintf(stderr, "Warning: dt=%.4f. For 3D stability, dt < 1/6 ≈ 0.1667 is required.\n", params.dt);
    }
    if (params.mode != 1 && params.mode != 2) {
        fprintf(stderr, "Error: mode must be 1 or 2\n"); return 1;
    }

    if (create_directory(output_dir) != 0) { fprintf(stderr, "Error: Cannot create output dir\n"); return 1; }

    files = (char (*)[MAX_PATH_LENGTH])malloc((size_t)MAX_FILES * MAX_PATH_LENGTH);
    if (!files) { fprintf(stderr, "Error: Memory allocation failed\n"); return 1; }
    if (get_tiff_files(input_dir, files, &file_count) != 0 || file_count == 0) {
        fprintf(stderr, "Error: No TIFF files found\n"); free(files); return 1;
    }
    if (get_image_info(input_dir, files[0], &info) != 0) {
        fprintf(stderr, "Error: Cannot read image info\n"); free(files); return 1;
    }
    info.depth = file_count;

    if (params.K <= 0.0f) {
        float sigma = estimate_noise_sigma(input_dir, files, file_count, &info);
        /* K ~ 2*sigma is a common heuristic for Perona-Malik */
        params.K = 2.0f * sigma;
        if (params.K <= 0.0f) params.K = 1.0f;
        fprintf(stderr, "  Estimated noise sigma: %.4f\n", sigma);
    }

    fprintf(stderr, "Anisotropic Diffusion 3D (Perona-Malik, GPU)\n");
    fprintf(stderr, "  Image: %u x %u x %u, %u-bit\n", info.width, info.height, info.depth, info.bits_per_sample);
    fprintf(stderr, "  Iterations: %d, K: %.4f, dt: %.4f, mode: %d\n",
            params.iterations, params.K, params.dt, params.mode);

    overlap_size = calculate_overlap_size();
    chunk_size = calculate_gpu_chunk_size(&info, overlap_size);
    valid_chunk_size = chunk_size - 2 * overlap_size;
    if (valid_chunk_size < 1) { valid_chunk_size = 1; chunk_size = 1 + 2 * overlap_size; }
    num_chunks = (file_count + valid_chunk_size - 1) / valid_chunk_size;

    fprintf(stderr, "  Overlap: %d, Chunk: %d, Chunks: %d\n", overlap_size, chunk_size, num_chunks);

    chunk = allocate_chunk(&info, chunk_size);
    if (!chunk) { fprintf(stderr, "Error: Memory allocation failed\n"); free(files); return 1; }

    for (chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {
        int vgs = chunk_idx * valid_chunk_size;
        int vge = vgs + valid_chunk_size - 1; if (vge >= file_count) vge = file_count - 1;
        int cgs = vgs - overlap_size; if (cgs < 0) cgs = 0;
        int cge = vge + overlap_size; if (cge >= file_count) cge = file_count - 1;
        chunk->start_z = cgs; chunk->end_z = cge;
        chunk->chunk_depth = cge - cgs + 1;
        chunk->valid_start = vgs - cgs; chunk->valid_end = vge - cgs;
        fprintf(stderr, "  Chunk %d/%d: z=%d..%d\n", chunk_idx+1, num_chunks, cgs, cge);
        if (load_chunk(input_dir, files, &info, chunk) != 0) {
            fprintf(stderr, "Error: Load failed\n"); free_chunk(chunk); free(files); return 1;
        }
        process_chunk_on_gpu(chunk, &info, &params);
        if (save_chunk_valid_region(input_dir, output_dir, files, &info, chunk) != 0) {
            fprintf(stderr, "Error: Save failed\n"); free_chunk(chunk); free(files); return 1;
        }
    }

    free_chunk(chunk); free(files);

    { FILE *f; int i;
      if ((f = fopen("cmd-hst.log", "a")) != NULL) {
          for (i = 0; i < argc; ++i) fprintf(f, "%s ", argv[i]);
          fprintf(f, "\n");
          fprintf(f, "   %% iterations %d  K %g  dt %g  mode %d\n", params.iterations, params.K, params.dt, params.mode);
          fclose(f);
    }}

    CUDA_CHECK(cudaDeviceReset());
    return 0;
}