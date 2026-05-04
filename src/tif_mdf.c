/*
 * tif_mdf.c - 3D Median Filter for TIFF Image Stacks (CPU version)
 *
 * Compile:
 *   Windows: cl /O2 /openmp tif_mdf.c libtiff.lib
 *   Linux:   gcc -O2 -fopenmp -o tif_mdf tif_mdf.c -ltiff -lm
 *
 * Usage:
 *   tif_mdf <input_dir> <output_dir> [kernel_size]
 *   Default kernel_size: 3
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

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

#ifdef _OPENMP
    #include <omp.h>
#endif

#define MAX_PATH_LENGTH 1024
#define MAX_FILES 100000
#define CHUNK_SIZE_MB 512

typedef struct {
    int kernel_size;
} FilterParams;

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

typedef struct {
    int start_z;
    int end_z;
    int chunk_depth;
    int valid_start;
    int valid_end;
    void *data;
} ChunkData;

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
    info->width = width; info->height = height;
    info->bits_per_sample = bits_per_sample; info->samples_per_pixel = samples_per_pixel;
    info->sample_format = sample_format;
    info->bytes_per_pixel = (bits_per_sample / 8) * samples_per_pixel;
    info->bytes_per_slice = (size_t)width * height * info->bytes_per_pixel;
    return 0;
}

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
    available_mb = (available_mb * 70) / 100;
    if (available_mb < CHUNK_SIZE_MB) available_mb = CHUNK_SIZE_MB;
    return available_mb;
}

static int calculate_overlap_size(FilterParams *params) {
    return params->kernel_size / 2 + 1;
}

static int calculate_chunk_size(ImageInfo *info, int overlap_size) {
    size_t available_mb = get_available_memory();
    size_t slice_mb = info->bytes_per_slice / (1024 * 1024);
    int chunk_depth, min_chunk_depth;
    if (slice_mb < 1) slice_mb = 1;
    chunk_depth = (int)(available_mb / (2 * slice_mb));
    min_chunk_depth = 4 * overlap_size;
    if (chunk_depth < min_chunk_depth) chunk_depth = min_chunk_depth;
    return chunk_depth;
}

static ChunkData* allocate_chunk(ImageInfo *info, int chunk_depth) {
    ChunkData *chunk;
    size_t chunk_size;
    chunk = (ChunkData*)malloc(sizeof(ChunkData));
    if (chunk == NULL) return NULL;
    chunk->chunk_depth = chunk_depth;
    chunk_size = (size_t)chunk_depth * info->bytes_per_slice;
    chunk->data = malloc(chunk_size);
    if (chunk->data == NULL) { free(chunk); return NULL; }
    return chunk;
}

static void free_chunk(ChunkData *chunk) {
    if (chunk != NULL) {
        if (chunk->data != NULL) free(chunk->data);
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
        slice_data = (char*)chunk->data + z * info->bytes_per_slice;
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
        slice_data = (char*)chunk->data + z * info->bytes_per_slice;
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

static void set_pixel_value(void *data, ImageInfo *info, int x, int y, int z, double value) {
    size_t offset;
    void *pixel_data;
    offset = z * info->bytes_per_slice + y * info->width * info->bytes_per_pixel +
             x * info->bytes_per_pixel;
    pixel_data = (char*)data + offset;
    if (info->bits_per_sample == 8) {
        if (value < 0.0) value = 0.0;
        if (value > 255.0) value = 255.0;
        *(unsigned char*)pixel_data = (unsigned char)(value + 0.5);
    } else if (info->bits_per_sample == 16) {
        if (value < 0.0) value = 0.0;
        if (value > 65535.0) value = 65535.0;
        *(unsigned short*)pixel_data = (unsigned short)(value + 0.5);
    } else if (info->bits_per_sample == 32 && info->sample_format == SAMPLEFORMAT_IEEEFP) {
        *(float*)pixel_data = (float)value;
    }
}

static int compare_double(const void *a, const void *b) {
    double va = *(const double*)a;
    double vb = *(const double*)b;
    if (va < vb) return -1;
    if (va > vb) return 1;
    return 0;
}

static void apply_median_filter_chunk(ChunkData *chunk, ImageInfo *info,
                                     FilterParams *params) {
    void *output_data;
    int x, y, z, kx, ky, kz;
    int half_kernel;
    int nx, ny, nz;
    size_t total_size;
    int max_neighbors;

    total_size = (size_t)chunk->chunk_depth * info->bytes_per_slice;
    output_data = malloc(total_size);
    if (output_data == NULL) return;
    memcpy(output_data, chunk->data, total_size);

    half_kernel = params->kernel_size / 2;
    max_neighbors = params->kernel_size * params->kernel_size * params->kernel_size;

#ifdef _OPENMP
    #pragma omp parallel private(y, x, kz, ky, kx, nz, ny, nx)
    {
#endif
    double *neighbors = (double*)malloc(max_neighbors * sizeof(double));

#ifdef _OPENMP
    #pragma omp for
#endif
    for (z = chunk->valid_start; z <= chunk->valid_end; z++) {
        for (y = 0; y < (int)info->height; y++) {
            for (x = 0; x < (int)info->width; x++) {
                int count = 0;
                for (kz = -half_kernel; kz <= half_kernel; kz++) {
                    nz = z + kz;
                    if (nz < 0 || nz >= chunk->chunk_depth) continue;
                    for (ky = -half_kernel; ky <= half_kernel; ky++) {
                        ny = y + ky;
                        if (ny < 0 || ny >= (int)info->height) continue;
                        for (kx = -half_kernel; kx <= half_kernel; kx++) {
                            nx = x + kx;
                            if (nx < 0 || nx >= (int)info->width) continue;
                            get_pixel_value(chunk->data, info, nx, ny, nz, &neighbors[count]);
                            count++;
                        }
                    }
                }
                qsort(neighbors, count, sizeof(double), compare_double);
                set_pixel_value(output_data, info, x, y, z, neighbors[count / 2]);
            }
        }
    }

    free(neighbors);

#ifdef _OPENMP
    }
#endif

    memcpy(chunk->data, output_data, total_size);
    free(output_data);
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
    int valid_global_start, valid_global_end;
    int chunk_global_start, chunk_global_end;
    int actual_chunk_depth;

    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input_dir> <output_dir> [kernel_size]\n", argv[0]);
        fprintf(stderr, "  Default kernel_size: 3 (odd, 3-21)\n");
        return 1;
    }
    strncpy(input_dir, argv[1], MAX_PATH_LENGTH - 1); input_dir[MAX_PATH_LENGTH - 1] = '\0';
    strncpy(output_dir, argv[2], MAX_PATH_LENGTH - 1); output_dir[MAX_PATH_LENGTH - 1] = '\0';

    params.kernel_size = 3;
    if (argc > 3) params.kernel_size = atoi(argv[3]);
    if (params.kernel_size < 3 || params.kernel_size > 21 || params.kernel_size % 2 == 0) {
        fprintf(stderr, "Error: Kernel size must be odd and between 3 and 21\n"); return 1;
    }

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
    chunk_size = calculate_chunk_size(&info, overlap_size);
    valid_chunk_size = chunk_size - 2 * overlap_size;
    if (valid_chunk_size < 1) { valid_chunk_size = 1; chunk_size = 1 + 2 * overlap_size; }
    num_chunks = (file_count + valid_chunk_size - 1) / valid_chunk_size;

    chunk = allocate_chunk(&info, chunk_size);
    if (chunk == NULL) { fprintf(stderr, "Error: Memory allocation failed\n"); free(files); return 1; }

    for (chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {
        valid_global_start = chunk_idx * valid_chunk_size;
        valid_global_end = valid_global_start + valid_chunk_size - 1;
        if (valid_global_end >= file_count) valid_global_end = file_count - 1;
        chunk_global_start = valid_global_start - overlap_size;
        chunk_global_end = valid_global_end + overlap_size;
        if (chunk_global_start < 0) chunk_global_start = 0;
        if (chunk_global_end >= file_count) chunk_global_end = file_count - 1;
        actual_chunk_depth = chunk_global_end - chunk_global_start + 1;
        chunk->start_z = chunk_global_start; chunk->end_z = chunk_global_end;
        chunk->chunk_depth = actual_chunk_depth;
        chunk->valid_start = valid_global_start - chunk_global_start;
        chunk->valid_end = valid_global_end - chunk_global_start;
        if (load_chunk(input_dir, files, &info, chunk) != 0) {
            fprintf(stderr, "Error: Failed to load chunk %d\n", chunk_idx);
            free_chunk(chunk); free(files); return 1;
        }
        apply_median_filter_chunk(chunk, &info, &params);
        if (save_chunk_valid_region(output_dir, files, &info, chunk) != 0) {
            fprintf(stderr, "Error: Failed to save chunk %d\n", chunk_idx);
            free_chunk(chunk); free(files); return 1;
        }
    }

    free_chunk(chunk);
    free(files);

    /* append to log file */
    {
        FILE *f;
        int i;
        if ((f = fopen("cmd-hst.log", "a")) != NULL) {
            for (i = 0; i < argc; ++i) fprintf(f, "%s ", argv[i]);
            fprintf(f, "\n");
            fprintf(f, "   %% kernel_size %d\n", params.kernel_size);
            fclose(f);
        }
    }

    return 0;
}