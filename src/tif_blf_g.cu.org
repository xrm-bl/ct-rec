/*
 * 3D Bilateral Filter for Large TIFF Image Stacks - CUDA Streaming Version
 * Handles datasets larger than GPU memory by processing in chunks
 * Optimized for maximum GPU utilization with minimal memory usage
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>

#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
    #include <direct.h>
    #include "tiffio.h"
    #define PATH_SEPARATOR "\\"
    #define mkdir(path, mode) _mkdir(path)
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

/* CUDA error checking macro */
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

/* Constants */
#define MAX_PATH_LENGTH 1024
#define MAX_FILES 10000
#define LOG_BUFFER_SIZE 4096
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#define BLOCK_SIZE_Z 4
#define MAX_KERNEL_SIZE 21
#define GPU_MEMORY_FRACTION 0.8f  /* Use 80% of GPU memory */

/* Structure for filter parameters */
typedef struct {
    int kernel_size;
    float spatial_sigma;
    float intensity_sigma;
} FilterParams;

/* Structure for image metadata */
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

/* Structure for GPU chunk processing */
typedef struct {
    int start_z;
    int end_z;
    int chunk_depth;
    void *h_data;      /* Host pinned memory */
    void *d_data;      /* Device memory */
    void *d_output;    /* Device output memory */
    cudaStream_t stream;
} GPUChunk;

/* Global variables for logging */
static FILE *log_file = NULL;

/* CUDA kernel for 3D bilateral filter - optimized for streaming */
template<typename T>
__global__ void bilateral_filter_3d_kernel_stream(
    const T* __restrict__ input,
    T* __restrict__ output,
    int width, int height, int chunk_depth,
    int kernel_size, 
    float spatial_sigma_sq_inv,
    float intensity_sigma_sq_inv,
    float max_value,
    int global_z_offset,
    int total_depth,
    const T* __restrict__ prev_overlap,
    const T* __restrict__ next_overlap,
    int prev_overlap_depth,
    int next_overlap_depth)
{
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int z = blockIdx.z * blockDim.z + threadIdx.z;
    
    if (x >= width || y >= height || z >= chunk_depth) return;
    
    const int half_kernel = kernel_size / 2;
    const size_t slice_size = (size_t)width * height;
    const size_t center_idx = (size_t)z * slice_size + (size_t)y * width + x;
    const float center_value = (float)input[center_idx];
    
    float weighted_sum = 0.0f;
    float weight_sum = 0.0f;
    
    /* Process kernel */
    for (int kz = -half_kernel; kz <= half_kernel; kz++) {
        const int local_z = z + kz;
        const int global_z = global_z_offset + z + kz;
        
        /* Skip if outside global bounds */
        if (global_z < 0 || global_z >= total_depth) continue;
        
        for (int ky = -half_kernel; ky <= half_kernel; ky++) {
            const int ny = y + ky;
            if (ny < 0 || ny >= height) continue;
            
            for (int kx = -half_kernel; kx <= half_kernel; kx++) {
                const int nx = x + kx;
                if (nx < 0 || nx >= width) continue;
                
                float neighbor_value;
                
                /* Determine which buffer to read from */
                if (local_z < 0 && prev_overlap != NULL) {
                    /* Read from previous overlap */
                    const int prev_z = prev_overlap_depth + local_z;
                    if (prev_z >= 0) {
                        const size_t prev_idx = (size_t)prev_z * slice_size + 
                                               (size_t)ny * width + nx;
                        neighbor_value = (float)prev_overlap[prev_idx];
                    } else {
                        continue;
                    }
                } else if (local_z >= chunk_depth && next_overlap != NULL) {
                    /* Read from next overlap */
                    const int next_z = local_z - chunk_depth;
                    if (next_z < next_overlap_depth) {
                        const size_t next_idx = (size_t)next_z * slice_size + 
                                               (size_t)ny * width + nx;
                        neighbor_value = (float)next_overlap[next_idx];
                    } else {
                        continue;
                    }
                } else if (local_z >= 0 && local_z < chunk_depth) {
                    /* Read from current chunk */
                    const size_t neighbor_idx = (size_t)local_z * slice_size + 
                                               (size_t)ny * width + nx;
                    neighbor_value = (float)input[neighbor_idx];
                } else {
                    continue;
                }
                
                /* Calculate weights */
                const float spatial_dist_sq = (float)(kx*kx + ky*ky + kz*kz);
                const float spatial_weight = __expf(-spatial_dist_sq * spatial_sigma_sq_inv);
                
                const float intensity_diff = neighbor_value - center_value;
                const float intensity_dist_sq = intensity_diff * intensity_diff;
                const float intensity_weight = __expf(-intensity_dist_sq * intensity_sigma_sq_inv);
                
                const float weight = spatial_weight * intensity_weight;
                weighted_sum += neighbor_value * weight;
                weight_sum += weight;
            }
        }
    }
    
    /* Write output */
    if (weight_sum > 0.0f) {
        float result = weighted_sum / weight_sum;
        result = fminf(fmaxf(result, 0.0f), max_value);
        output[center_idx] = (T)result;
    } else {
        output[center_idx] = input[center_idx];
    }
}

/* Function prototypes */
static void init_logging(const char *program_name);
static void close_logging(void);
static void log_message(const char *format, ...);
static int create_directory(const char *path);
static int get_tiff_files(const char *dir_path, char files[][MAX_PATH_LENGTH], int *file_count);
static int compare_strings(const void *a, const void *b);
static int get_image_info(const char *dir_path, const char *filename, ImageInfo *info);
static size_t get_gpu_memory_available(void);
static int calculate_gpu_chunk_size(ImageInfo *info, int kernel_size);
static GPUChunk* allocate_gpu_chunk(ImageInfo *info, int chunk_depth);
static void free_gpu_chunk(GPUChunk *chunk);
static int load_chunk_to_host(const char *dir_path, char files[][MAX_PATH_LENGTH], 
                              ImageInfo *info, GPUChunk *chunk);
static int save_chunk_from_host(const char *dir_path, char files[][MAX_PATH_LENGTH], 
                                ImageInfo *info, GPUChunk *chunk);
static void process_chunk_on_gpu(GPUChunk *chunk, ImageInfo *info, FilterParams *params,
                                GPUChunk *prev_chunk, GPUChunk *next_chunk);
static double get_current_time(void);
static float get_float_dynamic_range(const char *dir_path, char files[][MAX_PATH_LENGTH], 
                                    int file_count, ImageInfo *info);

/* Initialize logging */
static void init_logging(const char *program_name) {
    /*
    char log_path[MAX_PATH_LENGTH];
    char timestamp[64];
    time_t now;
    struct tm *tm_info;
    const char *home_dir;
    
#ifdef _WIN32
    home_dir = getenv("USERPROFILE");
#else
    home_dir = getenv("HOME");
#endif
    
    if (home_dir == NULL) {
        home_dir = ".";
    }
    
    snprintf(log_path, MAX_PATH_LENGTH, "%s%scom-log", home_dir, PATH_SEPARATOR);
    create_directory(log_path);
    
    time(&now);
    tm_info = localtime(&now);
    strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", tm_info);
    
    snprintf(log_path, MAX_PATH_LENGTH, "%s%scom-log%sbilateral_3d_cuda_large_%s.log", 
             home_dir, PATH_SEPARATOR, PATH_SEPARATOR, timestamp);
    
    log_file = fopen(log_path, "w");
    if (log_file == NULL) {
        fprintf(stderr, "Warning: Could not create log file\n");
    }
    
    log_message("=== 3D Bilateral Filter (CUDA Large Data) Started ===");
    log_message("Program: %s", program_name);
    
    // Log CUDA device information
    int device_count;
    cudaGetDeviceCount(&device_count);
    if (device_count > 0) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        log_message("CUDA Device: %s", prop.name);
        log_message("Total Global Memory: %.2f GB", prop.totalGlobalMem / (1024.0 * 1024.0 * 1024.0));
        log_message("Memory Bus Width: %d bits", prop.memoryBusWidth);
        log_message("Memory Bandwidth: %.2f GB/s", 
                    2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6);
    }
    */
}

/* Close logging */
static void close_logging(void) {
    /*
    if (log_file != NULL) {
        log_message("=== 3D Bilateral Filter (CUDA Large Data) Completed ===");
        fclose(log_file);
        log_file = NULL;
    }
    */
}

/* Log message with timestamp */
static void log_message(const char *format, ...) {
/* All logging disabled
    va_list args;
    time_t now;
    struct tm *tm_info;
    char timestamp[32];
    
    if (log_file == NULL) return;
    
    time(&now);
    tm_info = localtime(&now);
    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d %H:%M:%S", tm_info);
    
    fprintf(log_file, "[%s] ", timestamp);
    
    va_start(args, format);
    vfprintf(log_file, format, args);
    va_end(args);
    
    fprintf(log_file, "\n");
    fflush(log_file);
    
    printf("[%s] ", timestamp);
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
    printf("\n");
*/
}

/* Create directory */
static int create_directory(const char *path) {
#ifdef _WIN32
    if (GetFileAttributesA(path) == INVALID_FILE_ATTRIBUTES) {
        if (_mkdir(path) != 0) {
            return -1;
        }
    }
#else
    struct stat st;
    if (stat(path, &st) != 0) {
        if (mkdir(path, 0755) != 0) {
            return -1;
        }
    }
#endif
    return 0;
}

/* String comparison for qsort */
static int compare_strings(const void *a, const void *b) {
    return strcmp((const char *)a, (const char *)b);
}

/* Get list of TIFF files in directory */
static int get_tiff_files(const char *dir_path, char files[][MAX_PATH_LENGTH], int *file_count) {
#ifdef _WIN32
    WIN32_FIND_DATAA find_data;
    HANDLE find_handle;
    char search_path[MAX_PATH_LENGTH];
    
    snprintf(search_path, MAX_PATH_LENGTH, "%s\\*.tif*", dir_path);
    find_handle = FindFirstFileA(search_path, &find_data);
    
    if (find_handle == INVALID_HANDLE_VALUE) {
        return -1;
    }
    
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
    if (dir == NULL) {
        return -1;
    }
    
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

/* Get image information from first file */
static int get_image_info(const char *dir_path, const char *filename, ImageInfo *info) {
    TIFF *tif;
    char full_path[MAX_PATH_LENGTH];
    uint32 width, height;
    uint16 bits_per_sample, samples_per_pixel, sample_format;
    
    snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, filename);
    tif = TIFFOpen(full_path, "r");
    if (tif == NULL) {
        return -1;
    }
    
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);
    
    sample_format = SAMPLEFORMAT_UINT;
    TIFFGetFieldDefaulted(tif, TIFFTAG_SAMPLEFORMAT, &sample_format);
    
    TIFFClose(tif);
    
    info->width = width;
    info->height = height;
    info->bits_per_sample = bits_per_sample;
    info->samples_per_pixel = samples_per_pixel;
    info->sample_format = sample_format;
    info->bytes_per_pixel = (bits_per_sample / 8) * samples_per_pixel;
    info->bytes_per_slice = (size_t)width * height * info->bytes_per_pixel;
    
    return 0;
}

/* Get available GPU memory */
static size_t get_gpu_memory_available(void) {
    size_t free_mem, total_mem;
    
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
    
    /* Use only a fraction of available memory */
    free_mem = (size_t)(free_mem * GPU_MEMORY_FRACTION);
    
    return free_mem / (1024 * 1024); /* Return in MB */
}

/* Calculate optimal GPU chunk size */
static int calculate_gpu_chunk_size(ImageInfo *info, int kernel_size) {
    size_t available_mb = get_gpu_memory_available();
    size_t slice_mb = info->bytes_per_slice / (1024 * 1024);
    int half_kernel = kernel_size / 2;
    int chunk_depth;
    
    /* Need memory for: input chunk, output chunk, and overlap buffers */
    /* Factor of 3 for safety (input, output, working memory) */
    chunk_depth = (int)(available_mb / (3 * slice_mb));
    
    /* Ensure minimum chunk size */
    if (chunk_depth < kernel_size + 2) {
        chunk_depth = kernel_size + 2;
    }
    
    /* Limit maximum chunk size for better overlap handling */
    if (chunk_depth > 128) {
        chunk_depth = 128;
    }
    
    log_message("GPU available memory: %zu MB, Slice size: %zu MB, GPU chunk depth: %d", 
                available_mb, slice_mb, chunk_depth);
    
    return chunk_depth;
}

/* Allocate GPU chunk with streams */
static GPUChunk* allocate_gpu_chunk(ImageInfo *info, int chunk_depth) {
    GPUChunk *chunk;
    size_t chunk_size;
    
    chunk = (GPUChunk*)malloc(sizeof(GPUChunk));
    if (chunk == NULL) {
        return NULL;
    }
    
    chunk->chunk_depth = chunk_depth;
    chunk_size = (size_t)chunk_depth * info->bytes_per_slice;
    
    /* Allocate pinned host memory for async transfers */
    CUDA_CHECK(cudaMallocHost(&chunk->h_data, chunk_size));
    
    /* Allocate device memory */
    CUDA_CHECK(cudaMalloc(&chunk->d_data, chunk_size));
    CUDA_CHECK(cudaMalloc(&chunk->d_output, chunk_size));
    
    /* Create stream for async operations */
    CUDA_CHECK(cudaStreamCreate(&chunk->stream));
    
    return chunk;
}

/* Free GPU chunk */
static void free_gpu_chunk(GPUChunk *chunk) {
    if (chunk != NULL) {
        if (chunk->h_data != NULL) cudaFreeHost(chunk->h_data);
        if (chunk->d_data != NULL) cudaFree(chunk->d_data);
        if (chunk->d_output != NULL) cudaFree(chunk->d_output);
        if (chunk->stream != NULL) cudaStreamDestroy(chunk->stream);
        free(chunk);
    }
}

/* Load chunk from disk to pinned host memory */
static int load_chunk_to_host(const char *dir_path, char files[][MAX_PATH_LENGTH], 
                              ImageInfo *info, GPUChunk *chunk) {
    TIFF *tif;
    char full_path[MAX_PATH_LENGTH];
    int z, y;
    tsize_t scanline_size;
    unsigned char *buffer;
    void *slice_data;
    
    for (z = 0; z < chunk->chunk_depth && (chunk->start_z + z) < info->depth; z++) {
        snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, 
                 files[chunk->start_z + z]);
        tif = TIFFOpen(full_path, "r");
        if (tif == NULL) {
            log_message("Error: Cannot open %s", full_path);
            return -1;
        }
        
        scanline_size = TIFFScanlineSize(tif);
        buffer = (unsigned char*)_TIFFmalloc(scanline_size);
        if (buffer == NULL) {
            TIFFClose(tif);
            return -1;
        }
        
        slice_data = (char*)chunk->h_data + z * info->bytes_per_slice;
        
        for (y = 0; y < info->height; y++) {
            if (TIFFReadScanline(tif, buffer, y, 0) < 0) {
                _TIFFfree(buffer);
                TIFFClose(tif);
                return -1;
            }
            memcpy((char*)slice_data + y * scanline_size, buffer, scanline_size);
        }
        
        _TIFFfree(buffer);
        TIFFClose(tif);
    }
    
    chunk->end_z = chunk->start_z + z - 1;
    
    return 0;
}

/* Save chunk from pinned host memory to disk */
static int save_chunk_from_host(const char *dir_path, char files[][MAX_PATH_LENGTH], 
                                ImageInfo *info, GPUChunk *chunk) {
    TIFF *tif;
    char full_path[MAX_PATH_LENGTH];
    int z, y;
    tsize_t scanline_size;
    void *slice_data;
    
    for (z = 0; z < chunk->chunk_depth && (chunk->start_z + z) < info->depth; z++) {
        snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, 
                 files[chunk->start_z + z]);
        tif = TIFFOpen(full_path, "w");
        if (tif == NULL) {
            log_message("Error: Cannot create %s", full_path);
            return -1;
        }
        
        /* Set TIFF tags */
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
        
        for (y = 0; y < info->height; y++) {
            if (TIFFWriteScanline(tif, (char*)slice_data + y * scanline_size, y, 0) < 0) {
                TIFFClose(tif);
                return -1;
            }
        }
        
        TIFFClose(tif);
    }
    
    return 0;
}

/* Process chunk on GPU with overlap handling */
static void process_chunk_on_gpu(GPUChunk *chunk, ImageInfo *info, FilterParams *params,
                                GPUChunk *prev_chunk, GPUChunk *next_chunk) {
    dim3 block(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
    dim3 grid(
        (info->width + block.x - 1) / block.x,
        (info->height + block.y - 1) / block.y,
        (chunk->chunk_depth + block.z - 1) / block.z
    );
    
    float spatial_sigma_sq_inv = 0.5f / (params->spatial_sigma * params->spatial_sigma);
    float intensity_sigma_sq_inv = 0.5f / (params->intensity_sigma * params->intensity_sigma);
    float max_value;
    
    /* Set max value based on bit depth */
    if (info->bits_per_sample == 8) {
        max_value = 255.0f;
    } else if (info->bits_per_sample == 16) {
        max_value = 65535.0f;
    } else if (info->bits_per_sample == 32) {
        if (info->sample_format == SAMPLEFORMAT_IEEEFP) {
            max_value = FLT_MAX;
        } else {
            max_value = 4294967295.0f;
        }
    } else {
        max_value = FLT_MAX;
    }
    
    /* Prepare overlap data pointers */
    void *d_prev_overlap = NULL;
    void *d_next_overlap = NULL;
    int prev_overlap_depth = 0;
    int next_overlap_depth = 0;
    int half_kernel = params->kernel_size / 2;
    
    /* Copy overlap regions if available */
    if (prev_chunk != NULL && chunk->start_z > 0) {
        prev_overlap_depth = half_kernel;
        size_t overlap_offset = (prev_chunk->chunk_depth - prev_overlap_depth) * info->bytes_per_slice;
        d_prev_overlap = (char*)prev_chunk->d_output + overlap_offset;
    }
    
    if (next_chunk != NULL && chunk->end_z < info->depth - 1) {
        next_overlap_depth = half_kernel;
        d_next_overlap = next_chunk->d_data;
    }
    
    /* Launch kernel based on data type */
    if (info->bits_per_sample == 8) {
        bilateral_filter_3d_kernel_stream<unsigned char><<<grid, block, 0, chunk->stream>>>(
            (unsigned char*)chunk->d_data, (unsigned char*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            params->kernel_size, spatial_sigma_sq_inv, intensity_sigma_sq_inv, max_value,
            chunk->start_z, info->depth,
            (unsigned char*)d_prev_overlap, (unsigned char*)d_next_overlap,
            prev_overlap_depth, next_overlap_depth
        );
    } else if (info->bits_per_sample == 16) {
        bilateral_filter_3d_kernel_stream<unsigned short><<<grid, block, 0, chunk->stream>>>(
            (unsigned short*)chunk->d_data, (unsigned short*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            params->kernel_size, spatial_sigma_sq_inv, intensity_sigma_sq_inv, max_value,
            chunk->start_z, info->depth,
            (unsigned short*)d_prev_overlap, (unsigned short*)d_next_overlap,
            prev_overlap_depth, next_overlap_depth
        );
    } else if (info->bits_per_sample == 32 && info->sample_format == SAMPLEFORMAT_IEEEFP) {
        bilateral_filter_3d_kernel_stream<float><<<grid, block, 0, chunk->stream>>>(
            (float*)chunk->d_data, (float*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            params->kernel_size, spatial_sigma_sq_inv, intensity_sigma_sq_inv, max_value,
            chunk->start_z, info->depth,
            (float*)d_prev_overlap, (float*)d_next_overlap,
            prev_overlap_depth, next_overlap_depth
        );
    }
    
    /* Check for kernel errors */
    CUDA_CHECK(cudaGetLastError());
}

/* Get current time in seconds */
static double get_current_time(void) {
#ifdef _WIN32
    LARGE_INTEGER frequency, counter;
    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&counter);
    return (double)counter.QuadPart / (double)frequency.QuadPart;
#else
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
#endif
}

/* Get dynamic range for 32-bit float images by sampling */
static float get_float_dynamic_range(const char *dir_path, char files[][MAX_PATH_LENGTH], 
                                    int file_count, ImageInfo *info) {
    TIFF *tif;
    char full_path[MAX_PATH_LENGTH];
    tsize_t scanline_size;
    float *buffer;
    int sample_count, sample_interval;
    int i, y, x;
    float min_val = FLT_MAX;
    float max_val = -FLT_MAX;
    
    /* Sample 1/20 of images */
    sample_count = file_count / 20;
    if (sample_count < 1) sample_count = 1;
    if (sample_count > 10) sample_count = 10;  /* Limit to 10 samples */
    
    sample_interval = file_count / sample_count;
    
    log_message("Sampling %d images to determine dynamic range...", sample_count);
    
    for (i = 0; i < sample_count; i++) {
        int file_idx = i * sample_interval;
        snprintf(full_path, MAX_PATH_LENGTH, "%s%s%s", dir_path, PATH_SEPARATOR, files[file_idx]);
        
        tif = TIFFOpen(full_path, "r");
        if (tif == NULL) continue;
        
        scanline_size = TIFFScanlineSize(tif);
        buffer = (float*)_TIFFmalloc(scanline_size);
        if (buffer == NULL) {
            TIFFClose(tif);
            continue;
        }
        
        /* Sample center region of image */
        int start_y = info->height / 4;
        int end_y = 3 * info->height / 4;
        int y_step = (end_y - start_y) / 10;  /* Sample 10 lines per image */
        if (y_step < 1) y_step = 1;
        
        for (y = start_y; y < end_y; y += y_step) {
            if (TIFFReadScanline(tif, buffer, y, 0) >= 0) {
                int start_x = info->width / 4;
                int end_x = 3 * info->width / 4;
                for (x = start_x; x < end_x; x++) {
                    float val = buffer[x];
                    if (isfinite(val)) {  /* Skip NaN and Inf */
                        if (val < min_val) min_val = val;
                        if (val > max_val) max_val = val;
                    }
                }
            }
        }
        
        _TIFFfree(buffer);
        TIFFClose(tif);
    }
    
    float range = max_val - min_val;
    log_message("Detected dynamic range: [%.3f, %.3f], range: %.3f", 
                min_val, max_val, range);
    
    return range;
}

/* Main function */
int main(int argc, char *argv[]) {
    char input_dir[MAX_PATH_LENGTH];
    char output_dir[MAX_PATH_LENGTH];
    FilterParams params;
    char (*files)[MAX_PATH_LENGTH];
    int file_count;
    ImageInfo info;
    GPUChunk *chunks[3];  /* Triple buffering: processing, loading, saving */
    int chunk_size, num_chunks, chunk_idx;
    int half_kernel;
    double start_time, end_time, total_start_time;
    size_t total_size_mb, chunk_bytes;
    int current_chunk, next_chunk, prev_chunk;
    
    int	i;
    
    /* Initialize CUDA */
    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count == 0) {
        fprintf(stderr, "Error: No CUDA-capable devices found\n");
        return 1;
    }
    
    /* Select best device */
    int best_device = 0;
    int max_sm_count = 0;
    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, i));
        int sm_count = prop.multiProcessorCount;
        if (sm_count > max_sm_count) {
            max_sm_count = sm_count;
            best_device = i;
        }
    }
    CUDA_CHECK(cudaSetDevice(best_device));
    
    /* Initialize logging */
    init_logging(argv[0]);
    
    /* Parse command line arguments */
    if (argc < 3) {
        printf("Usage: %s <input_dir> <output_dir> [kernel_size] [spatial_sigma] "
                    "[intensity_sigma]", argv[0]);
//        close_logging();
        return 1;
    }
    
    strncpy(input_dir, argv[1], MAX_PATH_LENGTH - 1);
    input_dir[MAX_PATH_LENGTH - 1] = '\0';
    strncpy(output_dir, argv[2], MAX_PATH_LENGTH - 1);
    output_dir[MAX_PATH_LENGTH - 1] = '\0';
    
    /* Set default parameters */
    params.kernel_size = 5;
    params.spatial_sigma = 2.0f;
    params.intensity_sigma = 50.0f;
    
    /* Parse optional parameters */
    if (argc > 3) params.kernel_size = atoi(argv[3]);
    if (argc > 4) params.spatial_sigma = (float)atof(argv[4]);
    if (argc > 5) params.intensity_sigma = (float)atof(argv[5]);
    
    /* Validate parameters */
    if (params.kernel_size < 3 || params.kernel_size > MAX_KERNEL_SIZE || params.kernel_size % 2 == 0) {
        log_message("Error: Kernel size must be odd and between 3 and %d", MAX_KERNEL_SIZE);
        close_logging();
        return 1;
    }
    
    if (params.spatial_sigma <= 0 || params.intensity_sigma <= 0) {
        log_message("Error: Sigma values must be positive");
        close_logging();
        return 1;
    }
    
    half_kernel = params.kernel_size / 2;
    
    log_message("Input directory: %s", input_dir);
    log_message("Output directory: %s", output_dir);
    
    /* Create output directory */
    if (create_directory(output_dir) != 0) {
        log_message("Error: Cannot create output directory");
        close_logging();
        return 1;
    }
    
    /* Allocate file list */
    files = (char (*)[MAX_PATH_LENGTH])malloc(MAX_FILES * MAX_PATH_LENGTH);
    if (files == NULL) {
        log_message("Error: Memory allocation failed for file list");
        close_logging();
        return 1;
    }
    
    /* Get list of TIFF files */
    if (get_tiff_files(input_dir, files, &file_count) != 0) {
        log_message("Error: Cannot read input directory");
        free(files);
        close_logging();
        return 1;
    }
    
    if (file_count == 0) {
        log_message("Error: No TIFF files found in input directory");
        free(files);
        close_logging();
        return 1;
    }
    
    log_message("Found %d TIFF files", file_count);
    
    /* Get image information */
    if (get_image_info(input_dir, files[0], &info) != 0) {
        log_message("Error: Cannot read image information");
        free(files);
        close_logging();
        return 1;
    }
    
    info.depth = file_count;
    total_size_mb = (size_t)info.bytes_per_slice * file_count / (1024 * 1024);
    
    log_message("Image dimensions: %ux%ux%u", info.width, info.height, info.depth);
    log_message("Bits per sample: %u", info.bits_per_sample);
    log_message("Sample format: %u", info.sample_format);
    log_message("Total data size: %zu MB", total_size_mb);
    
    /* Auto-adjust intensity_sigma based on bit depth if using default value */
    if (argc <= 5 && params.intensity_sigma == 50.0f) {
        if (info.bits_per_sample == 8) {
            params.intensity_sigma = 256.0f / 3.0f;  /* ~85 */
            log_message("Auto-adjusted intensity_sigma to %.1f for 8-bit images", 
                        params.intensity_sigma);
        } else if (info.bits_per_sample == 16) {
            params.intensity_sigma = 65536.0f / 3.0f;  /* ~21845 */
            log_message("Auto-adjusted intensity_sigma to %.1f for 16-bit images", 
                        params.intensity_sigma);
        } else if (info.bits_per_sample == 32 && info.sample_format == SAMPLEFORMAT_IEEEFP) {
            /* For 32-bit float, sample images to determine dynamic range */
            float dynamic_range = get_float_dynamic_range(input_dir, files, file_count, &info);
            if (dynamic_range > 0.0f) {
                params.intensity_sigma = dynamic_range / 3.0f;
                log_message("Auto-adjusted intensity_sigma to %.3f for 32-bit float images", 
                            params.intensity_sigma);
            } else {
                params.intensity_sigma = 0.1f;  /* Fallback for normalized 0-1 range */
                log_message("Using default intensity_sigma %.3f for 32-bit float images", 
                            params.intensity_sigma);
            }
        }
    }
    
    /* Calculate chunk size */
    chunk_size = calculate_gpu_chunk_size(&info, params.kernel_size);
    num_chunks = (file_count + chunk_size - 1) / chunk_size;
    chunk_bytes = (size_t)chunk_size * info.bytes_per_slice;
    
    log_message("Processing in %d chunks of %d slices each", num_chunks, chunk_size);
    log_message("Chunk size: %.2f MB", chunk_bytes / (1024.0 * 1024.0));
    log_message("Kernel size: %d", params.kernel_size);
    log_message("Spatial sigma: %.2f", params.spatial_sigma);
    log_message("Intensity sigma: %.2f", params.intensity_sigma);
    
    total_start_time = get_current_time();
    
    /* Initialize chunks for triple buffering */
    for (int i = 0; i < 3; i++) {
        chunks[i] = allocate_gpu_chunk(&info, chunk_size + 2 * half_kernel);
        if (chunks[i] == NULL) {
            log_message("Error: Failed to allocate GPU chunk %d", i);
            for (int j = 0; j < i; j++) {
                free_gpu_chunk(chunks[j]);
            }
            free(files);
            close_logging();
            return 1;
        }
    }
    
    /* Process chunks with triple buffering */
    current_chunk = 0;
    next_chunk = 1;
    prev_chunk = 2;
    
    for (chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {
        int chunk_start, chunk_end, actual_chunk_depth;
        int overlap_start, overlap_end;
        
        log_message("Processing chunk %d/%d", chunk_idx + 1, num_chunks);
        
        /* Calculate chunk boundaries */
        chunk_start = chunk_idx * chunk_size;
        chunk_end = chunk_start + chunk_size - 1;
        if (chunk_end >= file_count) {
            chunk_end = file_count - 1;
        }
        
        /* Add overlap for kernel */
        overlap_start = chunk_start - half_kernel;
        overlap_end = chunk_end + half_kernel;
        
        if (overlap_start < 0) overlap_start = 0;
        if (overlap_end >= file_count) overlap_end = file_count - 1;
        
        actual_chunk_depth = overlap_end - overlap_start + 1;
        
        /* Set chunk parameters */
        chunks[current_chunk]->start_z = overlap_start;
        chunks[current_chunk]->chunk_depth = actual_chunk_depth;
        
        /* Load current chunk */
        start_time = get_current_time();
        if (load_chunk_to_host(input_dir, files, &info, chunks[current_chunk]) != 0) {
            log_message("Error: Failed to load chunk");
            for (int i = 0; i < 3; i++) {
                free_gpu_chunk(chunks[i]);
            }
            free(files);
            close_logging();
            return 1;
        }
        end_time = get_current_time();
        log_message("Loaded chunk in %.2f seconds", end_time - start_time);
        
        /* Transfer to GPU */
        start_time = get_current_time();
        CUDA_CHECK(cudaMemcpyAsync(chunks[current_chunk]->d_data, chunks[current_chunk]->h_data,
                                   chunk_bytes, cudaMemcpyHostToDevice, chunks[current_chunk]->stream));
        
        /* Process on GPU */
        process_chunk_on_gpu(chunks[current_chunk], &info, &params,
                            (chunk_idx > 0) ? chunks[prev_chunk] : NULL,
                            (chunk_idx < num_chunks - 1) ? chunks[next_chunk] : NULL);
        
        /* Transfer back to host */
        CUDA_CHECK(cudaMemcpyAsync(chunks[current_chunk]->h_data, chunks[current_chunk]->d_output,
                                   chunk_bytes, cudaMemcpyDeviceToHost, chunks[current_chunk]->stream));
        
        /* Synchronize and measure GPU time */
        CUDA_CHECK(cudaStreamSynchronize(chunks[current_chunk]->stream));
        end_time = get_current_time();
        log_message("GPU processing in %.2f seconds", end_time - start_time);
        
        /* Save chunk */
        start_time = get_current_time();
        if (save_chunk_from_host(output_dir, files, &info, chunks[current_chunk]) != 0) {
            log_message("Error: Failed to save chunk");
            for (int i = 0; i < 3; i++) {
                free_gpu_chunk(chunks[i]);
            }
            free(files);
            close_logging();
            return 1;
        }
        end_time = get_current_time();
        log_message("Saved chunk in %.2f seconds", end_time - start_time);
        
        /* Rotate chunk indices for triple buffering */
        int temp = prev_chunk;
        prev_chunk = current_chunk;
        current_chunk = next_chunk;
        next_chunk = temp;
    }
    
    /* Clean up */
    for (int i = 0; i < 3; i++) {
        free_gpu_chunk(chunks[i]);
    }
    free(files);
    
    /* Report total processing time */
    end_time = get_current_time();
    log_message("Total processing time: %.2f seconds", end_time - total_start_time);
    
    /* Calculate and log performance metrics */
    size_t total_pixels = (size_t)info.width * info.height * info.depth;
    double mpixels_per_second = total_pixels / ((end_time - total_start_time) * 1e6);
    log_message("Performance: %.2f MPixels/second", mpixels_per_second);
    
    log_message("Processing completed successfully");
    close_logging();
    
    /* Reset CUDA device */
    CUDA_CHECK(cudaDeviceReset());
    
// append to log file
	FILE		*f;
	if ((f = fopen("cmd-hst.log", "a")) == NULL) {
		return(-1);
	}
	for (i = 0; i<argc; ++i) fprintf(f, "%s ", argv[i]);
    fprintf(f,"  %% Kernel size: %d", params.kernel_size);
    fprintf(f," Spatial sigma: %.2f", params.spatial_sigma);
    fprintf(f," Intensity sigma: %.2f", params.intensity_sigma);
	fprintf(f, "\n");
	fclose(f);

    return 0;
}