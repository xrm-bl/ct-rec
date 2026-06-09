/*
 * tif_tvd_g.cu - 3D Total Variation Denoising for TIFF Image Stacks (CUDA)
 *
 * Algorithm: Chambolle-Pock (Primal-Dual) method for isotropic TV
 *
 *   min_u  (lambda/2) * ||u - f||^2  +  ||nabla u||_1
 *
 * Compile:
 *   Windows: nvcc -O3 -I%CUDAINCL% -o tif_tvd_g.exe tif_tvd_g.cu -use_fast_math -Xcompiler "/wd 4819" libtiff.lib
 *   Linux:   nvcc -O3 -o tif_tvd_g tif_tvd_g.cu -use_fast_math -ltiff
 *
 * Usage:
 *   tif_tvd_g <input_dir> <output_dir> [lambda] [iterations]
 *
 *   lambda:     Regularization strength (default: auto by bit depth)
 *               Larger  = more faithful to original (less denoising)
 *               Smaller = stronger smoothing
 *   iterations: Number of Chambolle-Pock iterations (default: 100)
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
#define TV_OVERLAP_MARGIN 2

typedef struct {
    float lambda;
    int iterations;
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
 * Chambolle-Pock CUDA Kernels
 *
 * Variables:
 *   f      - input (noisy) image (float, on GPU)
 *   u      - primal variable (denoised image)
 *   u_bar  - extrapolated primal variable
 *   px, py, pz - dual variables (gradient components)
 *
 * Step sizes:
 *   tau   = primal step size
 *   sigma = dual step size
 *   theta = extrapolation parameter (=1 for standard CP)
 *
 * Each iteration:
 *   1. Dual update:   p = proj( p + sigma * grad(u_bar) )
 *   2. Primal update: u_new = (u - tau * (-div(p)) + tau*lambda*f) / (1 + tau*lambda)
 *   3. Extrapolation: u_bar = u_new + theta * (u_new - u_old)
 * ---------------------------------------------------------------- */

/* Kernel: Dual update - compute gradient and update dual variables with projection */
__global__ void dual_update_kernel(
    const float* __restrict__ u_bar,
    float* __restrict__ px,
    float* __restrict__ py,
    float* __restrict__ pz,
    int width, int height, int depth,
    float sigma)
{
    ENABLE_SMEM_SPILLING();
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x >= width || y >= height || z >= depth) return;

    const size_t slice_size = (size_t)width * height;
    const size_t idx = (size_t)z * slice_size + (size_t)y * width + x;
    float ub = u_bar[idx];

    /* Forward differences for gradient */
    float gx = (x + 1 < width)  ? u_bar[idx + 1]          - ub : 0.0f;
    float gy = (y + 1 < height) ? u_bar[idx + width]       - ub : 0.0f;
    float gz = (z + 1 < depth)  ? u_bar[idx + slice_size]  - ub : 0.0f;

    /* Update dual: p = p + sigma * grad(u_bar) */
    float new_px = px[idx] + sigma * gx;
    float new_py = py[idx] + sigma * gy;
    float new_pz = pz[idx] + sigma * gz;

    /* Project onto unit ball (isotropic TV): p / max(1, |p|) */
    float norm = sqrtf(new_px * new_px + new_py * new_py + new_pz * new_pz);
    if (norm > 1.0f) {
        float inv_norm = 1.0f / norm;
        new_px *= inv_norm;
        new_py *= inv_norm;
        new_pz *= inv_norm;
    }

    px[idx] = new_px;
    py[idx] = new_py;
    pz[idx] = new_pz;
}

/* Kernel: Primal update - compute divergence and update u */
__global__ void primal_update_kernel(
    const float* __restrict__ f,
    const float* __restrict__ px,
    const float* __restrict__ py,
    const float* __restrict__ pz,
    float* __restrict__ u,
    float* __restrict__ u_bar,
    int width, int height, int depth,
    float tau, float lambda, float theta)
{
    ENABLE_SMEM_SPILLING();
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int z = blockIdx.z * blockDim.z + threadIdx.z;

    if (x >= width || y >= height || z >= depth) return;

    const size_t slice_size = (size_t)width * height;
    const size_t idx = (size_t)z * slice_size + (size_t)y * width + x;

    /* Backward differences for divergence: div(p) = dpx/dx + dpy/dy + dpz/dz */
    float div_px, div_py, div_pz;

    /* dpx/dx */
    if (x == 0) {
        div_px = px[idx];
    } else if (x == width - 1) {
        div_px = -px[idx - 1];
    } else {
        div_px = px[idx] - px[idx - 1];
    }

    /* dpy/dy */
    if (y == 0) {
        div_py = py[idx];
    } else if (y == height - 1) {
        div_py = -py[idx - width];
    } else {
        div_py = py[idx] - py[idx - width];
    }

    /* dpz/dz */
    if (z == 0) {
        div_pz = pz[idx];
    } else if (z == depth - 1) {
        div_pz = -pz[idx - slice_size];
    } else {
        div_pz = pz[idx] - pz[idx - slice_size];
    }

    float div_p = div_px + div_py + div_pz;

    /* Primal update: u_new = (u + tau * div(p) + tau*lambda*f) / (1 + tau*lambda) */
    float u_old = u[idx];
    float u_new = (u_old + tau * div_p + tau * lambda * f[idx]) / (1.0f + tau * lambda);

    /* Extrapolation: u_bar = u_new + theta * (u_new - u_old) */
    u[idx] = u_new;
    u_bar[idx] = u_new + theta * (u_new - u_old);
}

/* Kernel: Convert input data to float */
template<typename T>
__global__ void to_float_kernel(const T* __restrict__ input, float* __restrict__ output,
                                int width, int height, int depth) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int z = blockIdx.z * blockDim.z + threadIdx.z;
    if (x >= width || y >= height || z >= depth) return;
    const size_t idx = (size_t)z * (size_t)width * height + (size_t)y * width + x;
    output[idx] = (float)input[idx];
}

/* Kernel: Convert float result back to original type */
template<typename T>
__global__ void from_float_kernel(const float* __restrict__ input, T* __restrict__ output,
                                  int width, int height, int depth, float max_value) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int z = blockIdx.z * blockDim.z + threadIdx.z;
    if (x >= width || y >= height || z >= depth) return;
    const size_t idx = (size_t)z * (size_t)width * height + (size_t)y * width + x;
    float val = fminf(fmaxf(input[idx], 0.0f), max_value);
    output[idx] = (T)val;
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

static int calculate_overlap_size(void) {
    return TV_OVERLAP_MARGIN + 1;
}

static int calculate_gpu_chunk_size(ImageInfo *info, int overlap_size) {
    size_t free_mem, total_mem;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
    /* TV needs: d_data + d_output + d_f + d_u + d_u_bar + d_px + d_py + d_pz = 8 float buffers
     * Plus the original data buffer. Conservative: divide by 10 */
    size_t available_mb = (size_t)(free_mem * GPU_MEMORY_FRACTION) / (1024 * 1024);
    size_t float_slice_mb = ((size_t)info->width * info->height * sizeof(float)) / (1024 * 1024);
    if (float_slice_mb < 1) float_slice_mb = 1;
    int chunk_depth = (int)(available_mb / (10 * float_slice_mb));
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

/* ----------------------------------------------------------------
 * TV Denoising on GPU: convert to float, run CP iterations,
 * copy valid region back, convert to original type
 * ---------------------------------------------------------------- */
static void process_chunk_on_gpu(ChunkData *chunk, ImageInfo *info, FilterParams *params) {
    size_t chunk_bytes = (size_t)chunk->chunk_depth * info->bytes_per_slice;
    int w = info->width, h = info->height, d = chunk->chunk_depth;
    size_t num_voxels = (size_t)w * h * d;
    size_t float_bytes = num_voxels * sizeof(float);

    /* Upload raw data */
    CUDA_CHECK(cudaMemcpy(chunk->d_data, chunk->h_data, chunk_bytes, cudaMemcpyHostToDevice));

    /* Allocate float working buffers on GPU */
    float *d_f, *d_u, *d_u_bar, *d_px, *d_py, *d_pz;
    CUDA_CHECK(cudaMalloc(&d_f,     float_bytes));
    CUDA_CHECK(cudaMalloc(&d_u,     float_bytes));
    CUDA_CHECK(cudaMalloc(&d_u_bar, float_bytes));
    CUDA_CHECK(cudaMalloc(&d_px,    float_bytes));
    CUDA_CHECK(cudaMalloc(&d_py,    float_bytes));
    CUDA_CHECK(cudaMalloc(&d_pz,    float_bytes));

    /* Zero dual variables */
    CUDA_CHECK(cudaMemset(d_px, 0, float_bytes));
    CUDA_CHECK(cudaMemset(d_py, 0, float_bytes));
    CUDA_CHECK(cudaMemset(d_pz, 0, float_bytes));

    dim3 block(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
    dim3 grid(
        (w + block.x - 1) / block.x,
        (h + block.y - 1) / block.y,
        (d + block.z - 1) / block.z
    );

    /* Convert input to float and initialize u, u_bar */
    if (info->bits_per_sample == 8) {
        to_float_kernel<unsigned char><<<grid, block>>>(
            (unsigned char*)chunk->d_data, d_f, w, h, d);
    } else if (info->bits_per_sample == 16) {
        to_float_kernel<unsigned short><<<grid, block>>>(
            (unsigned short*)chunk->d_data, d_f, w, h, d);
    } else {
        to_float_kernel<float><<<grid, block>>>(
            (float*)chunk->d_data, d_f, w, h, d);
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    /* Initialize u = f, u_bar = f */
    CUDA_CHECK(cudaMemcpy(d_u,     d_f, float_bytes, cudaMemcpyDeviceToDevice));
    CUDA_CHECK(cudaMemcpy(d_u_bar, d_f, float_bytes, cudaMemcpyDeviceToDevice));

    /* Chambolle-Pock step sizes
     * For TV: L = sqrt(12) (operator norm of 3D gradient)
     * sigma * tau * L^2 < 1
     * Standard choice: tau = sigma = 1/L */
    float L = sqrtf(12.0f);
    float tau   = 1.0f / L;
    float sigma = 1.0f / L;
    float theta = 1.0f;

    /* Main iteration loop */
    for (int iter = 0; iter < params->iterations; iter++) {
        /* 1. Dual update */
        dual_update_kernel<<<grid, block>>>(
            d_u_bar, d_px, d_py, d_pz, w, h, d, sigma);

        /* 2. Primal update + extrapolation */
        primal_update_kernel<<<grid, block>>>(
            d_f, d_px, d_py, d_pz, d_u, d_u_bar,
            w, h, d, tau, params->lambda, theta);
    }
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    /* Convert back to original type */
    float max_value;
    if (info->bits_per_sample == 8) max_value = 255.0f;
    else if (info->bits_per_sample == 16) max_value = 65535.0f;
    else max_value = FLT_MAX;

    if (info->bits_per_sample == 8) {
        from_float_kernel<unsigned char><<<grid, block>>>(
            d_u, (unsigned char*)chunk->d_output, w, h, d, max_value);
    } else if (info->bits_per_sample == 16) {
        from_float_kernel<unsigned short><<<grid, block>>>(
            d_u, (unsigned short*)chunk->d_output, w, h, d, max_value);
    } else {
        from_float_kernel<float><<<grid, block>>>(
            d_u, (float*)chunk->d_output, w, h, d, max_value);
    }
    CUDA_CHECK(cudaDeviceSynchronize());

    /* Copy result back to host */
    CUDA_CHECK(cudaMemcpy(chunk->h_data, chunk->d_output, chunk_bytes, cudaMemcpyDeviceToHost));

    /* Free working buffers */
    cudaFree(d_f);
    cudaFree(d_u);
    cudaFree(d_u_bar);
    cudaFree(d_px);
    cudaFree(d_py);
    cudaFree(d_pz);
}

static float get_float_dynamic_range(const char *dir_path, char files[][MAX_PATH_LENGTH],
                                    int file_count, ImageInfo *info) {
    TIFF *tif;
    char full_path[MAX_PATH_LENGTH];
    tsize_t scanline_size;
    float *buffer;
    int sample_count = file_count / 20;
    if (sample_count < 1) sample_count = 1;
    if (sample_count > 10) sample_count = 10;
    int sample_interval = file_count / sample_count;
    float min_val = FLT_MAX, max_val = -FLT_MAX;
    for (int i = 0; i < sample_count; i++) {
        int file_idx = i * sample_interval;
        snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, files[file_idx]);
        tif = TIFFOpen(full_path, "r");
        if (tif == NULL) continue;
        scanline_size = TIFFScanlineSize(tif);
        buffer = (float*)_TIFFmalloc(scanline_size);
        if (buffer == NULL) { TIFFClose(tif); continue; }
        int start_y = info->height / 4, end_y = 3 * info->height / 4;
        int y_step = (end_y - start_y) / 10;
        if (y_step < 1) y_step = 1;
        for (int y = start_y; y < end_y; y += y_step) {
            if (TIFFReadScanline(tif, buffer, y, 0) >= 0) {
                int start_x = info->width / 4, end_x = 3 * info->width / 4;
                for (int x = start_x; x < end_x; x++) {
                    float val = buffer[x];
                    if (isfinite(val)) {
                        if (val < min_val) min_val = val;
                        if (val > max_val) max_val = val;
                    }
                }
            }
        }
        _TIFFfree(buffer);
        TIFFClose(tif);
    }
    return max_val - min_val;
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
        fprintf(stderr, "Usage: %s <input_dir> <output_dir> [lambda] [iterations]\n", argv[0]);
        fprintf(stderr, "  lambda:     Regularization (default: auto). Larger=less denoising\n");
        fprintf(stderr, "  iterations: CP iterations  (default: 100)\n");
        return 1;
    }

    strncpy(input_dir, argv[1], MAX_PATH_LENGTH - 1); input_dir[MAX_PATH_LENGTH - 1] = '\0';
    strncpy(output_dir, argv[2], MAX_PATH_LENGTH - 1); output_dir[MAX_PATH_LENGTH - 1] = '\0';

    params.lambda = -1.0f;
    params.iterations = 100;

    if (argc > 3) params.lambda = (float)atof(argv[3]);
    if (argc > 4) params.iterations = atoi(argv[4]);

    if (params.iterations < 1 || params.iterations > 10000) {
        fprintf(stderr, "Error: iterations must be between 1 and 10000\n"); return 1;
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

    /* Auto-estimate lambda if not specified */
    if (params.lambda <= 0.0f) {
        if (info.bits_per_sample == 8) {
            params.lambda = 10.0f;
        } else if (info.bits_per_sample == 16) {
            params.lambda = 10.0f;
        } else if (info.bits_per_sample == 32 && info.sample_format == SAMPLEFORMAT_IEEEFP) {
            float dr = get_float_dynamic_range(input_dir, files, file_count, &info);
            /* Scale lambda relative to dynamic range:
             * lambda ~ 10 * (dynamic_range / 256) gives comparable behavior across bit depths */
            params.lambda = (dr > 0.0f) ? 10.0f * (dr / 256.0f) : 10.0f;
        } else {
            params.lambda = 10.0f;
        }
    }

    fprintf(stderr, "TV Denoising 3D (Chambolle-Pock, GPU)\n");
    fprintf(stderr, "  Image: %u x %u x %u, %u-bit\n",
            info.width, info.height, info.depth, info.bits_per_sample);
    fprintf(stderr, "  Lambda: %.4f\n", params.lambda);
    fprintf(stderr, "  Iterations: %d\n", params.iterations);

    overlap_size = calculate_overlap_size();
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
            fprintf(f, "   %% lambda %g  iterations %d\n", params.lambda, params.iterations);
            fclose(f);
        }
    }

    CUDA_CHECK(cudaDeviceReset());
    return 0;
}
