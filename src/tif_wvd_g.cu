/*
 * tif_wvd_g.cu - 3D Wavelet Denoising for TIFF Image Stacks (CUDA)
 *
 * Algorithm: 3D Haar Discrete Wavelet Transform + BayesShrink soft thresholding
 *
 * Compile:
 *   Windows: nvcc -O3 -I%CUDAINCL% -o tif_wvd_g.exe tif_wvd_g.cu -use_fast_math -Xcompiler "/wd 4819" libtiff.lib
 *   Linux:   nvcc -O3 -o tif_wvd_g tif_wvd_g.cu -use_fast_math -ltiff
 *
 * Usage:
 *   tif_wvd_g <input_dir> <output_dir> [levels] [threshold_scale]
 *
 *   levels:          DWT decomposition levels (default: 3, range: 1-5)
 *   threshold_scale: Multiplier for BayesShrink threshold (default: 1.0)
 *                    >1.0 = stronger denoising, <1.0 = less denoising
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
    #include <stdint.h>
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
#define BLOCK_SIZE 256
#define BLOCK_SIZE_3D_X 8
#define BLOCK_SIZE_3D_Y 8
#define BLOCK_SIZE_3D_Z 8
#define GPU_MEMORY_FRACTION 0.7f
#define MAX_DWT_LEVELS 5

typedef struct {
    int levels;
    float threshold_scale;
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
 * 3D Haar DWT CUDA Kernels
 *
 * Forward DWT (analysis):
 *   For each axis, split signal into low-pass (average) and
 *   high-pass (detail) coefficients:
 *     low[i]  = (s[2i] + s[2i+1]) * 0.707107
 *     high[i] = (s[2i] - s[2i+1]) * 0.707107
 *
 * Inverse DWT (synthesis):
 *     s[2i]   = (low[i] + high[i]) * 0.707107
 *     s[2i+1] = (low[i] - high[i]) * 0.707107
 *
 * 3D DWT is performed as separable 1D transforms along X, Y, Z.
 * ---------------------------------------------------------------- */

#define SQRT2_INV 0.70710678118f  /* 1/sqrt(2) */

/* Forward Haar DWT along X axis */
__global__ void haar_forward_x(float *data, float *temp, int w, int h, int d, int len) {
    int y = blockIdx.x * blockDim.x + threadIdx.x;
    int z = blockIdx.y * blockDim.y + threadIdx.y;
    if (y >= h || z >= d) return;

    size_t slice = (size_t)w * h;
    size_t row_base = (size_t)z * slice + (size_t)y * w;
    int half = len / 2;

    for (int i = 0; i < half; i++) {
        float a = data[row_base + 2 * i];
        float b = (2 * i + 1 < len) ? data[row_base + 2 * i + 1] : a;
        temp[row_base + i]        = (a + b) * SQRT2_INV;
        temp[row_base + half + i] = (a - b) * SQRT2_INV;
    }

    for (int i = 0; i < len; i++) {
        data[row_base + i] = temp[row_base + i];
    }
}

/* Inverse Haar DWT along X axis */
__global__ void haar_inverse_x(float *data, float *temp, int w, int h, int d, int len) {
    int y = blockIdx.x * blockDim.x + threadIdx.x;
    int z = blockIdx.y * blockDim.y + threadIdx.y;
    if (y >= h || z >= d) return;

    size_t slice = (size_t)w * h;
    size_t row_base = (size_t)z * slice + (size_t)y * w;
    int half = len / 2;

    for (int i = 0; i < half; i++) {
        float lo = data[row_base + i];
        float hi = data[row_base + half + i];
        temp[row_base + 2 * i]     = (lo + hi) * SQRT2_INV;
        if (2 * i + 1 < len)
            temp[row_base + 2 * i + 1] = (lo - hi) * SQRT2_INV;
    }

    for (int i = 0; i < len; i++) {
        data[row_base + i] = temp[row_base + i];
    }
}

/* Forward Haar DWT along Y axis */
__global__ void haar_forward_y(float *data, float *temp, int w, int h, int d, int len) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int z = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= w || z >= d) return;

    size_t slice = (size_t)w * h;
    size_t base = (size_t)z * slice + x;
    int half = len / 2;

    for (int i = 0; i < half; i++) {
        float a = data[base + (size_t)(2 * i) * w];
        float b = (2 * i + 1 < len) ? data[base + (size_t)(2 * i + 1) * w] : a;
        temp[base + (size_t)i * w]          = (a + b) * SQRT2_INV;
        temp[base + (size_t)(half + i) * w] = (a - b) * SQRT2_INV;
    }

    for (int i = 0; i < len; i++) {
        data[base + (size_t)i * w] = temp[base + (size_t)i * w];
    }
}

/* Inverse Haar DWT along Y axis */
__global__ void haar_inverse_y(float *data, float *temp, int w, int h, int d, int len) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int z = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= w || z >= d) return;

    size_t slice = (size_t)w * h;
    size_t base = (size_t)z * slice + x;
    int half = len / 2;

    for (int i = 0; i < half; i++) {
        float lo = data[base + (size_t)i * w];
        float hi = data[base + (size_t)(half + i) * w];
        temp[base + (size_t)(2 * i) * w]     = (lo + hi) * SQRT2_INV;
        if (2 * i + 1 < len)
            temp[base + (size_t)(2 * i + 1) * w] = (lo - hi) * SQRT2_INV;
    }

    for (int i = 0; i < len; i++) {
        data[base + (size_t)i * w] = temp[base + (size_t)i * w];
    }
}

/* Forward Haar DWT along Z axis */
__global__ void haar_forward_z(float *data, float *temp, int w, int h, int d, int len) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= w || y >= h) return;

    size_t slice = (size_t)w * h;
    size_t base = (size_t)y * w + x;
    int half = len / 2;

    for (int i = 0; i < half; i++) {
        float a = data[base + (size_t)(2 * i) * slice];
        float b = (2 * i + 1 < len) ? data[base + (size_t)(2 * i + 1) * slice] : a;
        temp[base + (size_t)i * slice]          = (a + b) * SQRT2_INV;
        temp[base + (size_t)(half + i) * slice] = (a - b) * SQRT2_INV;
    }

    for (int i = 0; i < len; i++) {
        data[base + (size_t)i * slice] = temp[base + (size_t)i * slice];
    }
}

/* Inverse Haar DWT along Z axis */
__global__ void haar_inverse_z(float *data, float *temp, int w, int h, int d, int len) {
    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= w || y >= h) return;

    size_t slice = (size_t)w * h;
    size_t base = (size_t)y * w + x;
    int half = len / 2;

    for (int i = 0; i < half; i++) {
        float lo = data[base + (size_t)i * slice];
        float hi = data[base + (size_t)(half + i) * slice];
        temp[base + (size_t)(2 * i) * slice]     = (lo + hi) * SQRT2_INV;
        if (2 * i + 1 < len)
            temp[base + (size_t)(2 * i + 1) * slice] = (lo - hi) * SQRT2_INV;
    }

    for (int i = 0; i < len; i++) {
        data[base + (size_t)i * slice] = temp[base + (size_t)i * slice];
    }
}

/* Soft thresholding kernel applied to detail coefficients */
__global__ void soft_threshold_kernel(float *data, size_t n, float threshold) {
    size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;

    float val = data[idx];
    if (val > threshold) {
        data[idx] = val - threshold;
    } else if (val < -threshold) {
        data[idx] = val + threshold;
    } else {
        data[idx] = 0.0f;
    }
}

/* Compute sum of absolute values (for MAD estimation) - reduction kernel */
__global__ void abs_sum_kernel(const float *data, double *partial_sums,
                               size_t n) {
    __shared__ double sdata[BLOCK_SIZE];
    size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t tid = threadIdx.x;

    sdata[tid] = (idx < n) ? fabs((double)data[idx]) : 0.0;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }
    if (tid == 0) partial_sums[blockIdx.x] = sdata[0];
}

/* Compute sum of squares - reduction kernel */
__global__ void sq_sum_kernel(const float *data, double *partial_sums,
                              size_t n) {
    __shared__ double sdata[BLOCK_SIZE];
    size_t idx = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t tid = threadIdx.x;

    double val = (idx < n) ? (double)data[idx] : 0.0;
    sdata[tid] = val * val;
    __syncthreads();

    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        __syncthreads();
    }
    if (tid == 0) partial_sums[blockIdx.x] = sdata[0];
}

/* Convert input data to float */
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

/* Convert float back to original type */
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
 * GPU reduction helper: sum an array of doubles
 * ---------------------------------------------------------------- */
static double gpu_reduce_sum(double *d_partial, int num_blocks) {
    double *h_partial = (double*)malloc(num_blocks * sizeof(double));
    CUDA_CHECK(cudaMemcpy(h_partial, d_partial, num_blocks * sizeof(double), cudaMemcpyDeviceToHost));
    double total = 0.0;
    for (int i = 0; i < num_blocks; i++) total += h_partial[i];
    free(h_partial);
    return total;
}

/* ----------------------------------------------------------------
 * Estimate noise sigma from finest-level detail coefficients
 * using MAD (Median Absolute Deviation) approximation:
 *   sigma_noise = MAD / 0.6745 ≈ mean(|coeffs|) * sqrt(pi/2) / 0.6745
 *
 * For large sample sizes, mean(|X|) ≈ sigma * sqrt(2/pi) for
 * Gaussian X, so sigma ≈ mean(|X|) * sqrt(pi/2).
 * This avoids sorting for true median on GPU.
 * ---------------------------------------------------------------- */
static float estimate_noise_from_coeffs(float *d_coeffs, size_t n) {
    int num_blocks = (int)((n + BLOCK_SIZE - 1) / BLOCK_SIZE);
    double *d_partial;
    CUDA_CHECK(cudaMalloc(&d_partial, num_blocks * sizeof(double)));

    abs_sum_kernel<<<num_blocks, BLOCK_SIZE>>>(d_coeffs, d_partial, n);
    CUDA_CHECK(cudaDeviceSynchronize());
    double abs_sum = gpu_reduce_sum(d_partial, num_blocks);

    cudaFree(d_partial);

    double mean_abs = abs_sum / (double)n;
    /* sigma = mean_abs * sqrt(pi/2) */
    float sigma = (float)(mean_abs * 1.2533141373);
    return sigma;
}

/* ----------------------------------------------------------------
 * BayesShrink threshold for a subband:
 *   sigma_x = sqrt(max(sigma_y^2 - sigma_n^2, 0))
 *   threshold = sigma_n^2 / sigma_x
 *
 * sigma_n = noise sigma (estimated from finest detail)
 * sigma_y = std of the subband coefficients
 * ---------------------------------------------------------------- */
static float compute_bayes_threshold(float *d_coeffs, size_t n, float sigma_noise) {
    int num_blocks = (int)((n + BLOCK_SIZE - 1) / BLOCK_SIZE);
    double *d_partial;
    CUDA_CHECK(cudaMalloc(&d_partial, num_blocks * sizeof(double)));

    sq_sum_kernel<<<num_blocks, BLOCK_SIZE>>>(d_coeffs, d_partial, n);
    CUDA_CHECK(cudaDeviceSynchronize());
    double sq_sum = gpu_reduce_sum(d_partial, num_blocks);

    cudaFree(d_partial);

    double sigma_y_sq = sq_sum / (double)n;
    double sigma_n_sq = (double)sigma_noise * sigma_noise;
    double sigma_x_sq = sigma_y_sq - sigma_n_sq;

    if (sigma_x_sq <= 0.0) {
        /* Subband is pure noise; threshold everything */
        return (float)(sqrtf(sigma_y_sq) * 10.0f);
    }

    double sigma_x = sqrt(sigma_x_sq);
    float threshold = (float)(sigma_n_sq / sigma_x);
    return threshold;
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

/* DWT needs chunk dimensions to be multiples of 2^levels.
 * Overlap must cover boundary effects. */
static int calculate_overlap_size(FilterParams *params) {
    int margin = 1 << params->levels;  /* 2^levels */
    return margin + 1;
}

static int calculate_gpu_chunk_size(ImageInfo *info, int overlap_size, int levels) {
    size_t free_mem, total_mem;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
    /* Need: d_data + d_output + d_float + d_temp = 4 float-sized buffers + original */
    size_t available_mb = (size_t)(free_mem * GPU_MEMORY_FRACTION) / (1024 * 1024);
    size_t float_slice_mb = ((size_t)info->width * info->height * sizeof(float)) / (1024 * 1024);
    if (float_slice_mb < 1) float_slice_mb = 1;
    int chunk_depth = (int)(available_mb / (5 * float_slice_mb));
    int min_chunk_depth = 4 * overlap_size;
    if (chunk_depth < min_chunk_depth) chunk_depth = min_chunk_depth;
    /* Round down to multiple of 2^levels for clean DWT */
    int block = 1 << levels;
    chunk_depth = (chunk_depth / block) * block;
    if (chunk_depth < block) chunk_depth = block;
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

/* ----------------------------------------------------------------
 * Main wavelet denoising pipeline for one chunk:
 *   1. Convert to float
 *   2. Forward 3D Haar DWT (multi-level)
 *   3. Estimate noise sigma from finest detail coefficients
 *   4. BayesShrink soft thresholding on all detail subbands
 *   5. Inverse 3D Haar DWT
 *   6. Convert back to original type
 * ---------------------------------------------------------------- */
static void process_chunk_on_gpu(ChunkData *chunk, ImageInfo *info, FilterParams *params) {
    size_t chunk_bytes = (size_t)chunk->chunk_depth * info->bytes_per_slice;
    int w = info->width, h = info->height, d = chunk->chunk_depth;
    size_t num_voxels = (size_t)w * h * d;
    size_t float_bytes = num_voxels * sizeof(float);

    CUDA_CHECK(cudaMemcpy(chunk->d_data, chunk->h_data, chunk_bytes, cudaMemcpyHostToDevice));

    /* Allocate float buffers */
    float *d_float, *d_temp;
    CUDA_CHECK(cudaMalloc(&d_float, float_bytes));
    CUDA_CHECK(cudaMalloc(&d_temp,  float_bytes));

    dim3 block3d(BLOCK_SIZE_3D_X, BLOCK_SIZE_3D_Y, BLOCK_SIZE_3D_Z);
    dim3 grid3d(
        (w + block3d.x - 1) / block3d.x,
        (h + block3d.y - 1) / block3d.y,
        (d + block3d.z - 1) / block3d.z
    );

    /* Step 1: Convert to float */
    if (info->bits_per_sample == 8) {
        to_float_kernel<unsigned char><<<grid3d, block3d>>>(
            (unsigned char*)chunk->d_data, d_float, w, h, d);
    } else if (info->bits_per_sample == 16) {
        to_float_kernel<unsigned short><<<grid3d, block3d>>>(
            (unsigned short*)chunk->d_data, d_float, w, h, d);
    } else {
        to_float_kernel<float><<<grid3d, block3d>>>(
            (float*)chunk->d_data, d_float, w, h, d);
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    /* Step 2: Forward 3D Haar DWT (multi-level) */
    int lw = w, lh = h, ld = d;
    int levels = params->levels;

    /* Adjust levels if dimensions are too small */
    while (levels > 0 && (lw < 2 || lh < 2 || ld < 2)) levels--;
    for (int lv = 0; lv < levels; lv++) {
        if (lw < 2 || lh < 2 || ld < 2) break;

        /* Forward X */
        dim3 gx((lh + 15) / 16, (ld + 15) / 16);
        dim3 bx(16, 16);
        haar_forward_x<<<gx, bx>>>(d_float, d_temp, w, h, d, lw);
        CUDA_CHECK(cudaDeviceSynchronize());

        /* Forward Y */
        dim3 gy((lw + 15) / 16, (ld + 15) / 16);
        dim3 by(16, 16);
        haar_forward_y<<<gy, by>>>(d_float, d_temp, w, h, d, lh);
        CUDA_CHECK(cudaDeviceSynchronize());

        /* Forward Z */
        dim3 gz((lw + 15) / 16, (lh + 15) / 16);
        dim3 bz(16, 16);
        haar_forward_z<<<gz, bz>>>(d_float, d_temp, w, h, d, ld);
        CUDA_CHECK(cudaDeviceSynchronize());

        lw /= 2; lh /= 2; ld /= 2;
    }
    int actual_levels = 0;
    {
        int tw = w, th = h, td = d;
        for (int lv = 0; lv < levels; lv++) {
            if (tw < 2 || th < 2 || td < 2) break;
            tw /= 2; th /= 2; td /= 2;
            actual_levels++;
        }
    }

    /* Step 3: Estimate noise sigma from finest detail subband (HHH at level 1)
     * After one level of 3D DWT, the HHH subband occupies
     * the upper-right-back octant: x=[w/2..w), y=[h/2..h), z=[d/2..d)
     * We estimate sigma from this subband. */
    float sigma_noise;
    {
        int hw = w / 2, hh = h / 2, hd = d / 2;
        size_t detail_count = (size_t)(w - hw) * (h - hh) * (d - hd);

        /* Extract HHH subband into a contiguous buffer for reduction */
        float *d_hhh;
        CUDA_CHECK(cudaMalloc(&d_hhh, detail_count * sizeof(float)));

        /* Copy HHH subband to contiguous memory */
        size_t slice_full = (size_t)w * h;
        float *h_hhh_temp = (float*)malloc(detail_count * sizeof(float));
        float *h_full = (float*)malloc(float_bytes);
        CUDA_CHECK(cudaMemcpy(h_full, d_float, float_bytes, cudaMemcpyDeviceToHost));

        size_t cnt = 0;
        for (int z = hd; z < d; z++) {
            for (int y = hh; y < h; y++) {
                for (int x = hw; x < w; x++) {
                    h_hhh_temp[cnt++] = h_full[z * slice_full + y * w + x];
                }
            }
        }
        CUDA_CHECK(cudaMemcpy(d_hhh, h_hhh_temp, detail_count * sizeof(float), cudaMemcpyHostToDevice));
        free(h_hhh_temp);
        free(h_full);

        sigma_noise = estimate_noise_from_coeffs(d_hhh, detail_count);
        cudaFree(d_hhh);
    }

    fprintf(stderr, "    Noise sigma (from DWT): %.4f\n", sigma_noise);

    /* Step 4: BayesShrink soft thresholding on all detail subbands
     * At each level, there are 7 detail subbands (all octants except LLL).
     * We apply BayesShrink to each separately. */
    {
        float *h_data_full = (float*)malloc(float_bytes);
        CUDA_CHECK(cudaMemcpy(h_data_full, d_float, float_bytes, cudaMemcpyDeviceToHost));

        size_t slice_full = (size_t)w * h;
        int cw = w, ch = h, cd = d;

        for (int lv = 0; lv < actual_levels; lv++) {
            int half_w = cw / 2, half_h = ch / 2, half_d = cd / 2;

            /* 7 detail subbands: iterate over all octants except (0,0,0) = LLL */
            for (int oz = 0; oz <= 1; oz++) {
                for (int oy = 0; oy <= 1; oy++) {
                    for (int ox = 0; ox <= 1; ox++) {
                        if (ox == 0 && oy == 0 && oz == 0) continue; /* skip LLL */

                        int x0 = ox * half_w, x1 = (ox == 0) ? half_w : cw;
                        int y0 = oy * half_h, y1 = (oy == 0) ? half_h : ch;
                        int z0 = oz * half_d, z1 = (oz == 0) ? half_d : cd;
                        size_t sub_count = (size_t)(x1 - x0) * (y1 - y0) * (z1 - z0);

                        /* Extract subband */
                        float *d_sub;
                        CUDA_CHECK(cudaMalloc(&d_sub, sub_count * sizeof(float)));
                        float *h_sub = (float*)malloc(sub_count * sizeof(float));
                        size_t si = 0;
                        for (int z = z0; z < z1; z++) {
                            for (int y = y0; y < y1; y++) {
                                for (int x = x0; x < x1; x++) {
                                    h_sub[si++] = h_data_full[z * slice_full + y * w + x];
                                }
                            }
                        }
                        CUDA_CHECK(cudaMemcpy(d_sub, h_sub, sub_count * sizeof(float), cudaMemcpyHostToDevice));

                        /* Compute BayesShrink threshold */
                        float threshold = compute_bayes_threshold(d_sub, sub_count, sigma_noise);
                        threshold *= params->threshold_scale;

                        /* Apply soft thresholding */
                        int nblocks = (int)((sub_count + BLOCK_SIZE - 1) / BLOCK_SIZE);
                        soft_threshold_kernel<<<nblocks, BLOCK_SIZE>>>(d_sub, sub_count, threshold);
                        CUDA_CHECK(cudaDeviceSynchronize());

                        /* Write back */
                        CUDA_CHECK(cudaMemcpy(h_sub, d_sub, sub_count * sizeof(float), cudaMemcpyDeviceToHost));
                        si = 0;
                        for (int z = z0; z < z1; z++) {
                            for (int y = y0; y < y1; y++) {
                                for (int x = x0; x < x1; x++) {
                                    h_data_full[z * slice_full + y * w + x] = h_sub[si++];
                                }
                            }
                        }

                        free(h_sub);
                        cudaFree(d_sub);
                    }
                }
            }

            cw = half_w; ch = half_h; cd = half_d;
        }

        CUDA_CHECK(cudaMemcpy(d_float, h_data_full, float_bytes, cudaMemcpyHostToDevice));
        free(h_data_full);
    }

    /* Step 5: Inverse 3D Haar DWT (multi-level) */
    {
        /* Reconstruct sizes for each level */
        int sizes_w[MAX_DWT_LEVELS + 1], sizes_h[MAX_DWT_LEVELS + 1], sizes_d[MAX_DWT_LEVELS + 1];
        sizes_w[0] = w; sizes_h[0] = h; sizes_d[0] = d;
        for (int lv = 1; lv <= actual_levels; lv++) {
            sizes_w[lv] = sizes_w[lv-1] / 2;
            sizes_h[lv] = sizes_h[lv-1] / 2;
            sizes_d[lv] = sizes_d[lv-1] / 2;
        }

        for (int lv = actual_levels - 1; lv >= 0; lv--) {
            int rw = sizes_w[lv], rh = sizes_h[lv], rd = sizes_d[lv];

            /* Inverse Z */
            dim3 gz((rw + 15) / 16, (rh + 15) / 16);
            dim3 bz(16, 16);
            haar_inverse_z<<<gz, bz>>>(d_float, d_temp, w, h, d, rd);
            CUDA_CHECK(cudaDeviceSynchronize());

            /* Inverse Y */
            dim3 gy((rw + 15) / 16, (rd + 15) / 16);
            dim3 by(16, 16);
            haar_inverse_y<<<gy, by>>>(d_float, d_temp, w, h, d, rh);
            CUDA_CHECK(cudaDeviceSynchronize());

            /* Inverse X */
            dim3 gx((rh + 15) / 16, (rd + 15) / 16);
            dim3 bx(16, 16);
            haar_inverse_x<<<gx, bx>>>(d_float, d_temp, w, h, d, rw);
            CUDA_CHECK(cudaDeviceSynchronize());
        }
    }

    /* Step 6: Convert back to original type */
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
        from_float_kernel<unsigned char><<<grid3d, block3d>>>(
            d_float, (unsigned char*)chunk->d_output, w, h, d, min_value, max_value);
    } else if (info->bits_per_sample == 16) {
        from_float_kernel<unsigned short><<<grid3d, block3d>>>(
            d_float, (unsigned short*)chunk->d_output, w, h, d, min_value, max_value);
    } else {
        from_float_kernel<float><<<grid3d, block3d>>>(
            d_float, (float*)chunk->d_output, w, h, d, min_value, max_value);
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    CUDA_CHECK(cudaMemcpy(chunk->h_data, chunk->d_output, chunk_bytes, cudaMemcpyDeviceToHost));

    cudaFree(d_float);
    cudaFree(d_temp);
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
        fprintf(stderr, "Usage: %s <input_dir> <output_dir> [levels] [threshold_scale]\n", argv[0]);
        fprintf(stderr, "  levels:          DWT decomposition levels (default: 3, range: 1-5)\n");
        fprintf(stderr, "  threshold_scale: BayesShrink threshold multiplier (default: 1.0)\n");
        fprintf(stderr, "                   >1.0 = stronger denoising\n");
        fprintf(stderr, "                   <1.0 = less denoising\n");
        return 1;
    }

    strncpy(input_dir, argv[1], MAX_PATH_LENGTH - 1); input_dir[MAX_PATH_LENGTH - 1] = '\0';
    strncpy(output_dir, argv[2], MAX_PATH_LENGTH - 1); output_dir[MAX_PATH_LENGTH - 1] = '\0';

    params.levels = 3;
    params.threshold_scale = 1.0f;

    if (argc > 3) params.levels = atoi(argv[3]);
    if (argc > 4) params.threshold_scale = (float)atof(argv[4]);

    if (params.levels < 1 || params.levels > MAX_DWT_LEVELS) {
        fprintf(stderr, "Error: levels must be between 1 and %d\n", MAX_DWT_LEVELS); return 1;
    }
    if (params.threshold_scale <= 0.0f) {
        fprintf(stderr, "Error: threshold_scale must be positive\n"); return 1;
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

    fprintf(stderr, "Wavelet Denoising 3D (Haar + BayesShrink, GPU)\n");
    fprintf(stderr, "  Image: %u x %u x %u, %u-bit\n",
            info.width, info.height, info.depth, info.bits_per_sample);
    fprintf(stderr, "  Levels: %d\n", params.levels);
    fprintf(stderr, "  Threshold scale: %.2f\n", params.threshold_scale);

    overlap_size = calculate_overlap_size(&params);
    chunk_size = calculate_gpu_chunk_size(&info, overlap_size, params.levels);
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
            fprintf(f, "   %% levels %d  threshold_scale %g\n", params.levels, params.threshold_scale);
            fclose(f);
        }
    }

    CUDA_CHECK(cudaDeviceReset());
    return 0;
}