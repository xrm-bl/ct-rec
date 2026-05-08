/*
 * tif_gsf_g.cu - 3D Gaussian Filter for TIFF Image Stacks (CUDA version)
 *
 * Compile:
 *   nvcc -O3 -o tif_gsf_g tif_gsf_g.cu -ltiff
 *
 * Usage:
 *   tif_gsf_g <input_dir> <output_dir> [sigma]
 *   Default sigma: 2.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
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
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#define BLOCK_SIZE_Z 4
#define GPU_MEMORY_FRACTION 0.7f
#define MAX_KERNEL_WEIGHTS 4913  /* 17^3 */

__constant__ float d_weights[MAX_KERNEL_WEIGHTS];

typedef struct {
    float sigma;
    int kernel_size;
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

template<typename T>
__global__ void gaussian_filter_3d_kernel(
    const T* __restrict__ input,
    T* __restrict__ output,
    int width, int height, int chunk_depth,
    int valid_start, int valid_end,
    int kernel_size,
    float max_value)
{
    ENABLE_SMEM_SPILLING();
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int z = valid_start + blockIdx.z * blockDim.z + threadIdx.z;
    if (x >= width || y >= height || z > valid_end) return;
    const int half_kernel = kernel_size / 2;
    const size_t slice_size = (size_t)width * height;
    float weighted_sum = 0.0f;
    int widx = 0;
    for (int kz = -half_kernel; kz <= half_kernel; kz++) {
        int nz = z + kz;
        for (int ky = -half_kernel; ky <= half_kernel; ky++) {
            int ny = y + ky;
            for (int kx = -half_kernel; kx <= half_kernel; kx++) {
                int nx = x + kx;
                if (nz >= 0 && nz < chunk_depth &&
                    ny >= 0 && ny < height &&
                    nx >= 0 && nx < width) {
                    size_t neighbor_idx = (size_t)nz * slice_size + (size_t)ny * width + nx;
                    weighted_sum += (float)input[neighbor_idx] * d_weights[widx];
                }
                widx++;
            }
        }
    }
    size_t center_idx = (size_t)z * slice_size + (size_t)y * width + x;
    float result = fminf(fmaxf(weighted_sum, 0.0f), max_value);
    output[center_idx] = (T)result;
}

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

static int calculate_kernel_size(float sigma) {
    int half = (int)ceilf(3.0f * sigma);
    return 2 * half + 1;
}

static void upload_gaussian_kernel(int kernel_size, float sigma) {
    int half = kernel_size / 2;
    int total = kernel_size * kernel_size * kernel_size;
    float *weights = (float*)malloc(total * sizeof(float));
    float sigma_sq2 = 2.0f * sigma * sigma;
    float sum = 0.0f;
    int idx = 0;
    for (int kz = -half; kz <= half; kz++) {
        for (int ky = -half; ky <= half; ky++) {
            for (int kx = -half; kx <= half; kx++) {
                float dist_sq = (float)(kx*kx + ky*ky + kz*kz);
                weights[idx] = expf(-dist_sq / sigma_sq2);
                sum += weights[idx];
                idx++;
            }
        }
    }
    for (idx = 0; idx < total; idx++) weights[idx] /= sum;
    CUDA_CHECK(cudaMemcpyToSymbol(d_weights, weights, total * sizeof(float)));
    free(weights);
}

static int calculate_overlap_size(FilterParams *params) {
    int half_kernel = params->kernel_size / 2;
    int sigma_range = (int)ceilf(3.0f * params->sigma);
    int overlap = (half_kernel > sigma_range) ? half_kernel : sigma_range;
    return overlap + 1;
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

static int save_chunk_valid_region(const char *dir_path, char files[][MAX_PATH_LENGTH],
                                   ImageInfo *info, ChunkData *chunk) {
    TIFF *tif;
    char full_path[MAX_PATH_LENGTH];
    int z, y;
    tsize_t scanline_size;
    void *slice_data;
    int global_z;
    for (z = chunk->valid_start; z <= chunk->valid_end; z++) {
        global_z = chunk->start_z + z;
        if (global_z >= (int)info->depth) break;
        snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, files[global_z]);
        tif = TIFFOpen(full_path, "w");
        if (tif == NULL) return -1;
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
        gaussian_filter_3d_kernel<unsigned char><<<grid, block>>>(
            (unsigned char*)chunk->d_data, (unsigned char*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            chunk->valid_start, chunk->valid_end, params->kernel_size, max_value);
    } else if (info->bits_per_sample == 16) {
        gaussian_filter_3d_kernel<unsigned short><<<grid, block>>>(
            (unsigned short*)chunk->d_data, (unsigned short*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            chunk->valid_start, chunk->valid_end, params->kernel_size, max_value);
    } else if (info->bits_per_sample == 32 && info->sample_format == SAMPLEFORMAT_IEEEFP) {
        gaussian_filter_3d_kernel<float><<<grid, block>>>(
            (float*)chunk->d_data, (float*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            chunk->valid_start, chunk->valid_end, params->kernel_size, max_value);
    }
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(chunk->h_data, chunk->d_output, chunk_bytes, cudaMemcpyDeviceToHost));
}

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
        fprintf(stderr, "Usage: %s <input_dir> <output_dir> [sigma]\n", argv[0]); return 1;
    }
    strncpy(input_dir, argv[1], MAX_PATH_LENGTH - 1); input_dir[MAX_PATH_LENGTH - 1] = '\0';
    strncpy(output_dir, argv[2], MAX_PATH_LENGTH - 1); output_dir[MAX_PATH_LENGTH - 1] = '\0';

    params.sigma = 2.0f;
    if (argc > 3) params.sigma = (float)atof(argv[3]);
    if (params.sigma <= 0) { fprintf(stderr, "Error: Sigma must be positive\n"); return 1; }

    params.kernel_size = calculate_kernel_size(params.sigma);
    int total_weights = params.kernel_size * params.kernel_size * params.kernel_size;
    if (total_weights > MAX_KERNEL_WEIGHTS) {
        fprintf(stderr, "Error: Kernel too large (sigma=%.1f, kernel=%d, weights=%d > %d)\n",
                params.sigma, params.kernel_size, total_weights, MAX_KERNEL_WEIGHTS);
        return 1;
    }
    upload_gaussian_kernel(params.kernel_size, params.sigma);

    if (create_directory(output_dir) != 0) { fprintf(stderr, "Error: Cannot create output directory\n"); return 1; }

    files = (char (*)[MAX_PATH_LENGTH])malloc((size_t)MAX_FILES * MAX_PATH_LENGTH);
    if (files == NULL) { fprintf(stderr, "Error: Memory allocation failed\n"); return 1; }
    if (get_tiff_files(input_dir, files, &file_count) != 0 || file_count == 0) {
        fprintf(stderr, "Error: No TIFF files found\n"); free(files); return 1;
    }
    if (get_image_info(input_dir, files[0], &info) != 0) {
        fprintf(stderr, "Error: Cannot read image information\n"); free(files); return 1;
    }
    info.depth = file_count;

    overlap_size = calculate_overlap_size(&params);
    chunk_size = calculate_gpu_chunk_size(&info, overlap_size);
    valid_chunk_size = chunk_size - 2 * overlap_size;
    if (valid_chunk_size < 1) { valid_chunk_size = 1; chunk_size = 1 + 2 * overlap_size; }
    num_chunks = (file_count + valid_chunk_size - 1) / valid_chunk_size;

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
        if (load_chunk(input_dir, files, &info, chunk) != 0) {
            fprintf(stderr, "Error: Failed to load chunk %d\n", chunk_idx);
            free_chunk(chunk); free(files); return 1;
        }
        process_chunk_on_gpu(chunk, &info, &params);
        if (save_chunk_valid_region(output_dir, files, &info, chunk) != 0) {
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
            fprintf(f, "   %% sigma %g  kernel_size %d\n", params.sigma, params.kernel_size);
            fclose(f);
        }
    }

    CUDA_CHECK(cudaDeviceReset());
    return 0;
}