/*
 * tif_blf_g.cu - 3D Bilateral Filter for TIFF Image Stacks (CUDA version)
 *
 * Compile:
 *   nvcc -O3 -o tif_blf_g tif_blf_g.cu -ltiff
 *   nvcc -O3 -arch=sm_86 -o tif_blf_g tif_blf_g.cu -ltiff  (for RTX 3090)
 *
 * Usage:
 *   tif_blf_g <input_dir> <output_dir> [kernel_size] [spatial_sigma] [intensity_sigma]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <stdarg.h>
#include <cuda_runtime.h>

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
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#define BLOCK_SIZE_Z 4
#define MAX_KERNEL_SIZE 21
#define GPU_MEMORY_FRACTION 0.7f

typedef struct {
    int kernel_size;
    float spatial_sigma;
    float intensity_sigma;
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

/* CUDA kernel: bilateral filter for valid region of a chunk.
 * All kernel neighbors are guaranteed to be within the chunk
 * because overlap >= max(half_kernel, 3*spatial_sigma).
 */
template<typename T>
__global__ void bilateral_filter_3d_kernel(
    const T* __restrict__ input,
    T* __restrict__ output,
    int width, int height, int chunk_depth,
    int valid_start, int valid_end,
    int kernel_size,
    float spatial_sigma_sq_inv,
    float intensity_sigma_sq_inv,
    float max_value)
{
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    const int z = valid_start + blockIdx.z * blockDim.z + threadIdx.z;

    if (x >= width || y >= height || z > valid_end) return;

    const int half_kernel = kernel_size / 2;
    const size_t slice_size = (size_t)width * height;
    const size_t center_idx = (size_t)z * slice_size + (size_t)y * width + x;
    const float center_value = (float)input[center_idx];

    float weighted_sum = 0.0f;
    float weight_sum = 0.0f;

    for (int kz = -half_kernel; kz <= half_kernel; kz++) {
        const int nz = z + kz;
        if (nz < 0 || nz >= chunk_depth) continue;

        for (int ky = -half_kernel; ky <= half_kernel; ky++) {
            const int ny = y + ky;
            if (ny < 0 || ny >= height) continue;

            for (int kx = -half_kernel; kx <= half_kernel; kx++) {
                const int nx = x + kx;
                if (nx < 0 || nx >= width) continue;

                const size_t neighbor_idx = (size_t)nz * slice_size + (size_t)ny * width + nx;
                const float neighbor_value = (float)input[neighbor_idx];

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

    if (weight_sum > 0.0f) {
        float result = weighted_sum / weight_sum;
        result = fminf(fmaxf(result, 0.0f), max_value);
        output[center_idx] = (T)result;
    } else {
        output[center_idx] = input[center_idx];
    }
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
    uint32 width, height;
    uint16 bits_per_sample, samples_per_pixel, sample_format;
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
    info->width = width;
    info->height = height;
    info->bits_per_sample = bits_per_sample;
    info->samples_per_pixel = samples_per_pixel;
    info->sample_format = sample_format;
    info->bytes_per_pixel = (bits_per_sample / 8) * samples_per_pixel;
    info->bytes_per_slice = (size_t)width * height * info->bytes_per_pixel;
    return 0;
}

static int calculate_overlap_size(FilterParams *params) {
    int half_kernel = params->kernel_size / 2;
    int sigma_range = (int)ceilf(3.0f * params->spatial_sigma);
    int overlap = (half_kernel > sigma_range) ? half_kernel : sigma_range;
    return overlap + 1;
}

static int calculate_gpu_chunk_size(ImageInfo *info, int overlap_size) {
    size_t free_mem, total_mem;
    CUDA_CHECK(cudaMemGetInfo(&free_mem, &total_mem));
    size_t available_mb = (size_t)(free_mem * GPU_MEMORY_FRACTION) / (1024 * 1024);
    size_t slice_mb = info->bytes_per_slice / (1024 * 1024);
    if (slice_mb < 1) slice_mb = 1;
    /* Need d_data + d_output: factor of 2 */
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

    /* Transfer to GPU */
    CUDA_CHECK(cudaMemcpy(chunk->d_data, chunk->h_data, chunk_bytes, cudaMemcpyHostToDevice));

    /* Initialize output buffer with input (overlap regions stay as input) */
    CUDA_CHECK(cudaMemcpy(chunk->d_output, chunk->d_data, chunk_bytes, cudaMemcpyDeviceToDevice));

    /* Compute kernel parameters */
    float spatial_sigma_sq_inv = 0.5f / (params->spatial_sigma * params->spatial_sigma);
    float intensity_sigma_sq_inv = 0.5f / (params->intensity_sigma * params->intensity_sigma);
    float max_value;
    if (info->bits_per_sample == 8) max_value = 255.0f;
    else if (info->bits_per_sample == 16) max_value = 65535.0f;
    else max_value = FLT_MAX;

    /* Launch kernel - process only valid region */
    dim3 block(BLOCK_SIZE_X, BLOCK_SIZE_Y, BLOCK_SIZE_Z);
    dim3 grid(
        (info->width + block.x - 1) / block.x,
        (info->height + block.y - 1) / block.y,
        (valid_depth + block.z - 1) / block.z
    );

    if (info->bits_per_sample == 8) {
        bilateral_filter_3d_kernel<unsigned char><<<grid, block>>>(
            (unsigned char*)chunk->d_data, (unsigned char*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            chunk->valid_start, chunk->valid_end,
            params->kernel_size, spatial_sigma_sq_inv, intensity_sigma_sq_inv, max_value);
    } else if (info->bits_per_sample == 16) {
        bilateral_filter_3d_kernel<unsigned short><<<grid, block>>>(
            (unsigned short*)chunk->d_data, (unsigned short*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            chunk->valid_start, chunk->valid_end,
            params->kernel_size, spatial_sigma_sq_inv, intensity_sigma_sq_inv, max_value);
    } else if (info->bits_per_sample == 32 && info->sample_format == SAMPLEFORMAT_IEEEFP) {
        bilateral_filter_3d_kernel<float><<<grid, block>>>(
            (float*)chunk->d_data, (float*)chunk->d_output,
            info->width, info->height, chunk->chunk_depth,
            chunk->valid_start, chunk->valid_end,
            params->kernel_size, spatial_sigma_sq_inv, intensity_sigma_sq_inv, max_value);
    }

    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaDeviceSynchronize());

    /* Transfer back to host */
    CUDA_CHECK(cudaMemcpy(chunk->h_data, chunk->d_output, chunk_bytes, cudaMemcpyDeviceToHost));
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
        int start_y = info->height / 4;
        int end_y = 3 * info->height / 4;
        int y_step = (end_y - start_y) / 10;
        if (y_step < 1) y_step = 1;
        for (int y = start_y; y < end_y; y += y_step) {
            if (TIFFReadScanline(tif, buffer, y, 0) >= 0) {
                int start_x = info->width / 4;
                int end_x = 3 * info->width / 4;
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

    /* Initialize CUDA */
    int device_count;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count == 0) {
        fprintf(stderr, "Error: No CUDA-capable devices found\n");
        return 1;
    }

    /* Select best device by SM count */
    int best_device = 0;
    int max_sm_count = 0;
    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp prop;
        CUDA_CHECK(cudaGetDeviceProperties(&prop, i));
        if (prop.multiProcessorCount > max_sm_count) {
            max_sm_count = prop.multiProcessorCount;
            best_device = i;
        }
    }
    CUDA_CHECK(cudaSetDevice(best_device));

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input_dir> <output_dir> [kernel_size] [spatial_sigma] [intensity_sigma]\n", argv[0]);
        return 1;
    }

    strncpy(input_dir, argv[1], MAX_PATH_LENGTH - 1);
    input_dir[MAX_PATH_LENGTH - 1] = '\0';
    strncpy(output_dir, argv[2], MAX_PATH_LENGTH - 1);
    output_dir[MAX_PATH_LENGTH - 1] = '\0';

    params.kernel_size = 5;
    params.spatial_sigma = 2.0f;
    params.intensity_sigma = 50.0f;

    if (argc > 3) params.kernel_size = atoi(argv[3]);
    if (argc > 4) params.spatial_sigma = (float)atof(argv[4]);
    if (argc > 5) params.intensity_sigma = (float)atof(argv[5]);

    if (params.kernel_size < 3 || params.kernel_size > MAX_KERNEL_SIZE || params.kernel_size % 2 == 0) {
        fprintf(stderr, "Error: Kernel size must be odd and between 3 and %d\n", MAX_KERNEL_SIZE);
        return 1;
    }
    if (params.spatial_sigma <= 0 || params.intensity_sigma <= 0) {
        fprintf(stderr, "Error: Sigma values must be positive\n");
        return 1;
    }

    if (create_directory(output_dir) != 0) {
        fprintf(stderr, "Error: Cannot create output directory\n");
        return 1;
    }

    files = (char (*)[MAX_PATH_LENGTH])malloc((size_t)MAX_FILES * MAX_PATH_LENGTH);
    if (files == NULL) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        return 1;
    }

    if (get_tiff_files(input_dir, files, &file_count) != 0) {
        fprintf(stderr, "Error: Cannot read input directory\n");
        free(files); return 1;
    }
    if (file_count == 0) {
        fprintf(stderr, "Error: No TIFF files found\n");
        free(files); return 1;
    }
    if (get_image_info(input_dir, files[0], &info) != 0) {
        fprintf(stderr, "Error: Cannot read image information\n");
        free(files); return 1;
    }
    info.depth = file_count;

    /* Auto-adjust intensity_sigma based on bit depth */
    if (argc <= 5 && params.intensity_sigma == 50.0f) {
        if (info.bits_per_sample == 8) {
            params.intensity_sigma = 256.0f / 3.0f;
        } else if (info.bits_per_sample == 16) {
            params.intensity_sigma = 65536.0f / 3.0f;
        } else if (info.bits_per_sample == 32 && info.sample_format == SAMPLEFORMAT_IEEEFP) {
            float dynamic_range = get_float_dynamic_range(input_dir, files, file_count, &info);
            if (dynamic_range > 0.0f) {
                params.intensity_sigma = dynamic_range / 3.0f;
            } else {
                params.intensity_sigma = 0.1f;
            }
        }
    }

    overlap_size = calculate_overlap_size(&params);
    chunk_size = calculate_gpu_chunk_size(&info, overlap_size);
    valid_chunk_size = chunk_size - 2 * overlap_size;
    if (valid_chunk_size < 1) {
        valid_chunk_size = 1;
        chunk_size = 1 + 2 * overlap_size;
    }
    num_chunks = (file_count + valid_chunk_size - 1) / valid_chunk_size;

    chunk = allocate_chunk(&info, chunk_size);
    if (chunk == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for chunk\n");
        free(files); return 1;
    }

    for (chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {
        int valid_global_start = chunk_idx * valid_chunk_size;
        int valid_global_end = valid_global_start + valid_chunk_size - 1;
        if (valid_global_end >= file_count) valid_global_end = file_count - 1;
        int chunk_global_start = valid_global_start - overlap_size;
        int chunk_global_end = valid_global_end + overlap_size;
        if (chunk_global_start < 0) chunk_global_start = 0;
        if (chunk_global_end >= file_count) chunk_global_end = file_count - 1;
        int actual_chunk_depth = chunk_global_end - chunk_global_start + 1;

        chunk->start_z = chunk_global_start;
        chunk->end_z = chunk_global_end;
        chunk->chunk_depth = actual_chunk_depth;
        chunk->valid_start = valid_global_start - chunk_global_start;
        chunk->valid_end = valid_global_end - chunk_global_start;

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

    free_chunk(chunk);
    free(files);
    CUDA_CHECK(cudaDeviceReset());

// append to log file
	FILE		*f;
	int			i;
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