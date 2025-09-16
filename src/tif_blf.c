/*
 * 3D Bilateral Filter for Large TIFF Image Stacks - Streaming Version
 * Handles datasets larger than available RAM by processing in chunks
 * C89 compliant implementation
 * Supports Windows (Visual C++) and Linux (GCC)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

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

#ifdef _OPENMP
    #include <omp.h>
#endif

/* Constants */
#define MAX_PATH_LENGTH 1024
#define MAX_FILES 10000
#define LOG_BUFFER_SIZE 4096
#define CHUNK_SIZE_MB 512  /* Process in 512MB chunks */

/* Structure for filter parameters */
typedef struct {
    int kernel_size;
    double spatial_sigma;
    double intensity_sigma;
} FilterParams;

/* Structure for image metadata */
typedef struct {
    uint32 width;
    uint32 height;
    uint32 depth;
    uint16 bits_per_sample;
    uint16 samples_per_pixel;
    uint16 sample_format;
    size_t bytes_per_pixel;
    size_t bytes_per_slice;
} ImageInfo;

/* Structure for chunk processing */
typedef struct {
    int start_z;
    int end_z;
    int chunk_depth;
    void *data;
} ChunkData;

/* Global variables for logging */
static FILE *log_file = NULL;

/* Function prototypes */
static void init_logging(const char *program_name);
static void close_logging(void);
static void log_message(const char *format, ...);
static int create_directory(const char *path);
static int get_tiff_files(const char *dir_path, char files[][MAX_PATH_LENGTH], int *file_count);
static int compare_strings(const void *a, const void *b);
static int get_image_info(const char *dir_path, const char *filename, ImageInfo *info);
static size_t get_available_memory(void);
static int calculate_chunk_size(ImageInfo *info, int kernel_size);
static ChunkData* allocate_chunk(ImageInfo *info, int chunk_depth);
static void free_chunk(ChunkData *chunk);
static int load_chunk(const char *dir_path, char files[][MAX_PATH_LENGTH], 
                      ImageInfo *info, ChunkData *chunk);
static int save_chunk(const char *dir_path, char files[][MAX_PATH_LENGTH], 
                      ImageInfo *info, ChunkData *chunk);
static void apply_bilateral_filter_chunk(ChunkData *chunk, ImageInfo *info, 
                                        FilterParams *params, ChunkData *prev_chunk, 
                                        ChunkData *next_chunk);
static double gaussian_weight(double distance, double sigma);
static void get_pixel_value(void *data, ImageInfo *info, int x, int y, int z, double *value);
static void set_pixel_value(void *data, ImageInfo *info, int x, int y, int z, double value);
static double get_current_time(void);
static double get_float_dynamic_range(const char *dir_path, char files[][MAX_PATH_LENGTH], 
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
    
    snprintf(log_path, MAX_PATH_LENGTH, "%s%scom-log%sbilateral_3d_large_%s.log", 
             home_dir, PATH_SEPARATOR, PATH_SEPARATOR, timestamp);
    
    log_file = fopen(log_path, "w");
    if (log_file == NULL) {
        fprintf(stderr, "Warning: Could not create log file\n");
    }
    
    log_message("=== 3D Bilateral Filter (Large Data) Started ===");
    log_message("Program: %s", program_name);
    log_message("Chunk size: %d MB", CHUNK_SIZE_MB);
    */
}

/* Close logging */
static void close_logging(void) {
    /*
    if (log_file != NULL) {
        log_message("=== 3D Bilateral Filter (Large Data) Completed ===");
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
    
    /* Set default sample format if not specified */
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

/* Get available system memory */
static size_t get_available_memory(void) {
    size_t available_mb;
    
#ifdef _WIN32
    MEMORYSTATUSEX mem_status;
    mem_status.dwLength = sizeof(mem_status);
    GlobalMemoryStatusEx(&mem_status);
    available_mb = (size_t)(mem_status.ullAvailPhys / (1024 * 1024));
#else
    struct sysinfo info;
    sysinfo(&info);
    available_mb = (size_t)((info.freeram * info.mem_unit) / (1024 * 1024));
#endif
    
    /* Use only 70% of available memory to be safe */
    available_mb = (available_mb * 70) / 100;
    
    /* But not less than CHUNK_SIZE_MB */
    if (available_mb < CHUNK_SIZE_MB) {
        available_mb = CHUNK_SIZE_MB;
    }
    
    return available_mb;
}

/* Calculate optimal chunk size based on available memory and kernel size */
static int calculate_chunk_size(ImageInfo *info, int kernel_size) {
    size_t available_mb = get_available_memory();
    size_t slice_mb = info->bytes_per_slice / (1024 * 1024);
    int half_kernel = kernel_size / 2;
    int chunk_depth;
    
    /* We need space for chunk + overlap regions */
    /* Minimum chunk depth should be kernel_size to ensure proper filtering */
    chunk_depth = (int)(available_mb / (3 * slice_mb)); /* 3x for input, output, temp */
    
    if (chunk_depth < kernel_size) {
        chunk_depth = kernel_size;
    }
    
    /* Add overlap for kernel */
    chunk_depth = chunk_depth - 2 * half_kernel;
    if (chunk_depth < 1) {
        chunk_depth = 1;
    }
    
    /* log_message("Available memory: %zu MB, Slice size: %zu MB, Chunk depth: %d", 
                available_mb, slice_mb, chunk_depth); */
    
    return chunk_depth;
}

/* Allocate memory for a chunk */
static ChunkData* allocate_chunk(ImageInfo *info, int chunk_depth) {
    ChunkData *chunk;
    size_t chunk_size;
    
    chunk = (ChunkData*)malloc(sizeof(ChunkData));
    if (chunk == NULL) {
        return NULL;
    }
    
    chunk->chunk_depth = chunk_depth;
    chunk_size = (size_t)chunk_depth * info->bytes_per_slice;
    
    chunk->data = malloc(chunk_size);
    if (chunk->data == NULL) {
        free(chunk);
        return NULL;
    }
    
    return chunk;
}

/* Free chunk memory */
static void free_chunk(ChunkData *chunk) {
    if (chunk != NULL) {
        if (chunk->data != NULL) {
            free(chunk->data);
        }
        free(chunk);
    }
}

/* Load a chunk of slices from disk */
static int load_chunk(const char *dir_path, char files[][MAX_PATH_LENGTH], 
                      ImageInfo *info, ChunkData *chunk) {
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
        
        slice_data = (char*)chunk->data + z * info->bytes_per_slice;
        
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

/* Save a chunk of slices to disk */
static int save_chunk(const char *dir_path, char files[][MAX_PATH_LENGTH], 
                      ImageInfo *info, ChunkData *chunk) {
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
        
        /* Only set sample format if it's a valid value */
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
        
        slice_data = (char*)chunk->data + z * info->bytes_per_slice;
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

/* Gaussian weight function */
static double gaussian_weight(double distance, double sigma) {
    return exp(-(distance * distance) / (2.0 * sigma * sigma));
}

/* Get pixel value */
static void get_pixel_value(void *data, ImageInfo *info, int x, int y, int z, double *value) {
    size_t offset;
    void *pixel_data;
    
    offset = z * info->bytes_per_slice + y * info->width * info->bytes_per_pixel + 
             x * info->bytes_per_pixel;
    pixel_data = (char*)data + offset;
    
    if (info->bits_per_sample == 8) {
        *value = (double)(*(unsigned char*)pixel_data);
    } else if (info->bits_per_sample == 16) {
        *value = (double)(*(unsigned short*)pixel_data);
    } else if (info->bits_per_sample == 32 && info->sample_format == SAMPLEFORMAT_IEEEFP) {
        *value = (double)(*(float*)pixel_data);
    }
}

/* Set pixel value */
static void set_pixel_value(void *data, ImageInfo *info, int x, int y, int z, double value) {
    size_t offset;
    void *pixel_data;
    
    offset = z * info->bytes_per_slice + y * info->width * info->bytes_per_pixel + 
             x * info->bytes_per_pixel;
    pixel_data = (char*)data + offset;
    
    if (info->bits_per_sample == 8) {
        *(unsigned char*)pixel_data = (unsigned char)(value + 0.5);
    } else if (info->bits_per_sample == 16) {
        *(unsigned short*)pixel_data = (unsigned short)(value + 0.5);
    } else if (info->bits_per_sample == 32 && info->sample_format == SAMPLEFORMAT_IEEEFP) {
        *(float*)pixel_data = (float)value;
    }
}

/* Apply bilateral filter to a chunk with overlap handling */
static void apply_bilateral_filter_chunk(ChunkData *chunk, ImageInfo *info, 
                                        FilterParams *params, ChunkData *prev_chunk, 
                                        ChunkData *next_chunk) {
    void *output_data;
    int x, y, z, kx, ky, kz;
    int half_kernel;
    double center_value, neighbor_value;
    double spatial_dist, intensity_dist;
    double spatial_weight, intensity_weight, weight;
    double weighted_sum, weight_sum;
    double dx, dy, dz;
    int nx, ny, nz;
    int local_z;
    
    /* Allocate temporary output buffer */
    output_data = malloc((size_t)chunk->chunk_depth * info->bytes_per_slice);
    if (output_data == NULL) {
        log_message("Error: Memory allocation failed for output buffer");
        return;
    }
    
    half_kernel = params->kernel_size / 2;
    
    /* Process each voxel in the chunk */
#ifdef _OPENMP
    #pragma omp parallel for private(y, x, kz, ky, kx, center_value, weighted_sum, weight_sum, nz, ny, nx, neighbor_value, local_z, dz, dy, dx, spatial_dist, spatial_weight, intensity_dist, intensity_weight, weight)
#endif
    for (z = 0; z < chunk->chunk_depth; z++) {
        for (y = 0; y < info->height; y++) {
            for (x = 0; x < info->width; x++) {
                get_pixel_value(chunk->data, info, x, y, z, &center_value);
                weighted_sum = 0.0;
                weight_sum = 0.0;
                
                /* Apply kernel */
                for (kz = -half_kernel; kz <= half_kernel; kz++) {
                    nz = z + kz;
                    local_z = nz;
                    
                    /* Handle overlap with previous chunk */
                    if (nz < 0 && prev_chunk != NULL) {
                        local_z = prev_chunk->chunk_depth + nz;
                        /* Check bounds in previous chunk */
                        if (local_z >= 0) {
                            for (ky = -half_kernel; ky <= half_kernel; ky++) {
                                ny = y + ky;
                                if (ny < 0 || ny >= info->height) continue;
                                
                                for (kx = -half_kernel; kx <= half_kernel; kx++) {
                                    nx = x + kx;
                                    if (nx < 0 || nx >= info->width) continue;
                                    
                                    get_pixel_value(prev_chunk->data, info, nx, ny, local_z, 
                                                   &neighbor_value);
                                    
                                    dx = (double)kx;
                                    dy = (double)ky;
                                    dz = (double)kz;
                                    spatial_dist = sqrt(dx*dx + dy*dy + dz*dz);
                                    spatial_weight = gaussian_weight(spatial_dist, 
                                                                   params->spatial_sigma);
                                    
                                    intensity_dist = fabs(neighbor_value - center_value);
                                    intensity_weight = gaussian_weight(intensity_dist, 
                                                                     params->intensity_sigma);
                                    
                                    weight = spatial_weight * intensity_weight;
                                    weighted_sum += neighbor_value * weight;
                                    weight_sum += weight;
                                }
                            }
                        }
                    }
                    /* Handle overlap with next chunk */
                    else if (nz >= chunk->chunk_depth && next_chunk != NULL) {
                        local_z = nz - chunk->chunk_depth;
                        /* Check bounds in next chunk */
                        if (local_z < next_chunk->chunk_depth) {
                            for (ky = -half_kernel; ky <= half_kernel; ky++) {
                                ny = y + ky;
                                if (ny < 0 || ny >= info->height) continue;
                                
                                for (kx = -half_kernel; kx <= half_kernel; kx++) {
                                    nx = x + kx;
                                    if (nx < 0 || nx >= info->width) continue;
                                    
                                    get_pixel_value(next_chunk->data, info, nx, ny, local_z, 
                                                   &neighbor_value);
                                    
                                    dx = (double)kx;
                                    dy = (double)ky;
                                    dz = (double)kz;
                                    spatial_dist = sqrt(dx*dx + dy*dy + dz*dz);
                                    spatial_weight = gaussian_weight(spatial_dist, 
                                                                   params->spatial_sigma);
                                    
                                    intensity_dist = fabs(neighbor_value - center_value);
                                    intensity_weight = gaussian_weight(intensity_dist, 
                                                                     params->intensity_sigma);
                                    
                                    weight = spatial_weight * intensity_weight;
                                    weighted_sum += neighbor_value * weight;
                                    weight_sum += weight;
                                }
                            }
                        }
                    }
                    /* Handle current chunk */
                    else if (nz >= 0 && nz < chunk->chunk_depth) {
                        for (ky = -half_kernel; ky <= half_kernel; ky++) {
                            ny = y + ky;
                            if (ny < 0 || ny >= info->height) continue;
                            
                            for (kx = -half_kernel; kx <= half_kernel; kx++) {
                                nx = x + kx;
                                if (nx < 0 || nx >= info->width) continue;
                                
                                get_pixel_value(chunk->data, info, nx, ny, nz, &neighbor_value);
                                
                                dx = (double)kx;
                                dy = (double)ky;
                                dz = (double)kz;
                                spatial_dist = sqrt(dx*dx + dy*dy + dz*dz);
                                spatial_weight = gaussian_weight(spatial_dist, 
                                                               params->spatial_sigma);
                                
                                intensity_dist = fabs(neighbor_value - center_value);
                                intensity_weight = gaussian_weight(intensity_dist, 
                                                                 params->intensity_sigma);
                                
                                weight = spatial_weight * intensity_weight;
                                weighted_sum += neighbor_value * weight;
                                weight_sum += weight;
                            }
                        }
                    }
                }
                
                /* Set filtered value */
                if (weight_sum > 0.0) {
                    set_pixel_value(output_data, info, x, y, z, weighted_sum / weight_sum);
                } else {
                    set_pixel_value(output_data, info, x, y, z, center_value);
                }
            }
        }
    }
    
    /* Copy filtered data back */
    memcpy(chunk->data, output_data, (size_t)chunk->chunk_depth * info->bytes_per_slice);
    
    /* Free temporary buffer */
    free(output_data);
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
static double get_float_dynamic_range(const char *dir_path, char files[][MAX_PATH_LENGTH], 
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
    
    double range = (double)(max_val - min_val);
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
    ChunkData *current_chunk, *prev_chunk, *next_chunk;
    int chunk_size, num_chunks, chunk_idx;
    int half_kernel;
    double start_time, end_time, total_start_time;
    size_t total_size_mb;
    
	int	i;
	
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
    params.spatial_sigma = 2.0;
    params.intensity_sigma = 50.0;
    
    /* Parse optional parameters */
    if (argc > 3) params.kernel_size = atoi(argv[3]);
    if (argc > 4) params.spatial_sigma = atof(argv[4]);
    if (argc > 5) params.intensity_sigma = atof(argv[5]);
    
    /* Validate parameters */
    if (params.kernel_size < 3 || params.kernel_size > 21 || params.kernel_size % 2 == 0) {
        log_message("Error: Kernel size must be odd and between 3 and 21");
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
    
    /* Auto-adjust intensity_sigma based on bit depth if using default value */
    if (argc <= 5 && params.intensity_sigma == 50.0) {
        if (info.bits_per_sample == 8) {
            params.intensity_sigma = 256.0 / 3.0;  /* ~85 */
            log_message("Auto-adjusted intensity_sigma to %.1f for 8-bit images", 
                        params.intensity_sigma);
        } else if (info.bits_per_sample == 16) {
            params.intensity_sigma = 65536.0 / 3.0;  /* ~21845 */
            log_message("Auto-adjusted intensity_sigma to %.1f for 16-bit images", 
                        params.intensity_sigma);
        } else if (info.bits_per_sample == 32 && info.sample_format == SAMPLEFORMAT_IEEEFP) {
            /* For 32-bit float, sample images to determine dynamic range */
            double dynamic_range = get_float_dynamic_range(input_dir, files, file_count, &info);
            if (dynamic_range > 0.0) {
                params.intensity_sigma = dynamic_range / 3.0;
                log_message("Auto-adjusted intensity_sigma to %.3f for 32-bit float images", 
                            params.intensity_sigma);
            } else {
                params.intensity_sigma = 0.1;  /* Fallback for normalized 0-1 range */
                log_message("Using default intensity_sigma %.3f for 32-bit float images", 
                            params.intensity_sigma);
            }
        }
    }
    
    log_message("Image dimensions: %ux%ux%u", info.width, info.height, info.depth);
    log_message("Bits per sample: %u", info.bits_per_sample);
    log_message("Sample format: %u", info.sample_format);
    log_message("Total data size: %zu MB", total_size_mb);
    
    /* Calculate chunk size */
    chunk_size = calculate_chunk_size(&info, params.kernel_size);
    num_chunks = (file_count + chunk_size - 1) / chunk_size;
    
    log_message("Processing in %d chunks of %d slices each", num_chunks, chunk_size);
    log_message("Kernel size: %d", params.kernel_size);
    log_message("Spatial sigma: %.2f", params.spatial_sigma);
    log_message("Intensity sigma: %.2f", params.intensity_sigma);
    
    total_start_time = get_current_time();
    
    /* Initialize chunk pointers */
    prev_chunk = NULL;
    current_chunk = NULL;
    next_chunk = NULL;
    
    /* Process chunks */
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
        
        /* Add overlap for kernel (except at boundaries) */
        overlap_start = chunk_start - half_kernel;
        overlap_end = chunk_end + half_kernel;
        
        if (overlap_start < 0) overlap_start = 0;
        if (overlap_end >= file_count) overlap_end = file_count - 1;
        
        actual_chunk_depth = overlap_end - overlap_start + 1;
        
        /* Rotate chunk pointers */
        if (prev_chunk != NULL) {
            free_chunk(prev_chunk);
        }
        prev_chunk = current_chunk;
        current_chunk = next_chunk;
        next_chunk = NULL;
        
        /* Allocate and load current chunk if not already loaded */
        if (current_chunk == NULL) {
            current_chunk = allocate_chunk(&info, actual_chunk_depth);
            if (current_chunk == NULL) {
                log_message("Error: Memory allocation failed for chunk");
                if (prev_chunk != NULL) free_chunk(prev_chunk);
                free(files);
                close_logging();
                return 1;
            }
            
            current_chunk->start_z = overlap_start;
            
            start_time = get_current_time();
            if (load_chunk(input_dir, files, &info, current_chunk) != 0) {
                log_message("Error: Failed to load chunk");
                free_chunk(current_chunk);
                if (prev_chunk != NULL) free_chunk(prev_chunk);
                free(files);
                close_logging();
                return 1;
            }
            end_time = get_current_time();
            log_message("Loaded chunk in %.2f seconds", end_time - start_time);
        }
        
        /* Pre-load next chunk for overlap if not last chunk */
        if (chunk_idx < num_chunks - 1) {
            int next_overlap_start = chunk_end + 1 - half_kernel;
            int next_overlap_end = chunk_end + chunk_size + half_kernel;
            int next_actual_depth;
            
            if (next_overlap_start < 0) next_overlap_start = 0;
            if (next_overlap_end >= file_count) next_overlap_end = file_count - 1;
            
            next_actual_depth = next_overlap_end - next_overlap_start + 1;
            
            next_chunk = allocate_chunk(&info, next_actual_depth);
            if (next_chunk != NULL) {
                next_chunk->start_z = next_overlap_start;
                load_chunk(input_dir, files, &info, next_chunk);
            }
        }
        
        /* Apply bilateral filter */
        start_time = get_current_time();
        apply_bilateral_filter_chunk(current_chunk, &info, &params, prev_chunk, next_chunk);
        end_time = get_current_time();
        log_message("Filtered chunk in %.2f seconds", end_time - start_time);
        
        /* Save chunk (skip overlap regions) */
        start_time = get_current_time();
        if (save_chunk(output_dir, files, &info, current_chunk) != 0) {
            log_message("Error: Failed to save chunk");
            free_chunk(current_chunk);
            if (prev_chunk != NULL) free_chunk(prev_chunk);
            if (next_chunk != NULL) free_chunk(next_chunk);
            free(files);
            close_logging();
            return 1;
        }
        end_time = get_current_time();
        log_message("Saved chunk in %.2f seconds", end_time - start_time);
    }
    
    /* Clean up remaining chunks */
    if (prev_chunk != NULL) free_chunk(prev_chunk);
    if (current_chunk != NULL) free_chunk(current_chunk);
    if (next_chunk != NULL) free_chunk(next_chunk);
    
    /* Free file list */
    free(files);
    
    /* Report total processing time */
    end_time = get_current_time();
    log_message("Total processing time: %.2f seconds", end_time - total_start_time);
    log_message("Processing completed successfully");
    
    close_logging();

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