#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "sort_filter_omp.h"

/*----------------------------------------------------------------------*/
/* ring removal, Vo et al., Opt. Express 26, 28396-28412 (2018)         */
/* most of the code were written by Claude 3.7 sonnet                   */

/*
How to set environment parameter

tcsh .tcshrc 
setenv KERNEL_SIZE 5
setenv OMP_NUM_THREADS 40

bash .bashrc 
export KERNEL_SIZE=5
export OMP_NUM_THREADS=40

windows
set KERNEL_SIZE=5
set OMP_NUM_THREADS=40
*/

/* Structure to store pixel value and its original position */
typedef struct {
    float value;
    unsigned int original_pos;
} PixelInfo;

/* Comparison function for PixelInfo structure sorting */
static int compare_pixel_info(const void *a, const void *b) {
    PixelInfo *pa = (PixelInfo *)a;
    PixelInfo *pb = (PixelInfo *)b;
    return (pa->value > pb->value) - (pa->value < pb->value);
}

/* Comparison function for float values */
static int compare_floats(const void *a, const void *b) {
    float fa = *(const float*)a;
    float fb = *(const float*)b;
    return (fa > fb) - (fa < fb);
}

/* Apply median filter to a pixel */
float apply_median_filter(float *row, unsigned int Nx, unsigned int x, int kernel_size) {
    int half_size;
    float *kernel_values;
    float median_value;
    int count = 0;
    int i;
    int pos;
    
    half_size = kernel_size / 2;
    kernel_values = (float *)malloc(kernel_size * sizeof(float));
    
    /* Collect values within kernel */
    for (i = -half_size; i <= half_size; i++) {
        pos = x + i;
        
        /* Mirror pixels at boundaries */
        if (pos < 0) {
            pos = -pos;
        } else if (pos >= Nx) {
            pos = 2 * Nx - pos - 1;
        }
        
        kernel_values[count++] = row[pos];
    }
    
    /* Sort kernel values */
    qsort(kernel_values, kernel_size, sizeof(float), compare_floats);
    
    /* Get median value */
    median_value = kernel_values[kernel_size / 2];
    
    /* Free memory */
    free(kernel_values);
    
    return median_value;
}

/* Process a single column - extract, sort, and store original positions */
void process_column(float *image_data, float *sorted_data, PixelInfo *column_info, 
                    int col, unsigned int Nx, unsigned int Ny) {
    unsigned int row;
    
    /* Extract column data and record original positions */
    for (row = 0; row < Ny; row++) {
        column_info[row].value = image_data[row * Nx + col];
        column_info[row].original_pos = row;
    }
    
    /* Sort column data */
    qsort(column_info, Ny, sizeof(PixelInfo), compare_pixel_info);
    
    /* Store sorted values */
    for (row = 0; row < Ny; row++) {
        sorted_data[row * Nx + col] = column_info[row].value;
    }
}

/* Restore a single column to original positions */
void restore_column(float *image_data, float *filtered_data, float *result_data,
                   PixelInfo *column_info, int col, unsigned int Nx, unsigned int Ny) {
    unsigned int row;
    unsigned int sorted_row;
    unsigned int original_row;
    
    /* Re-create original position information */
    for (row = 0; row < Ny; row++) {
        column_info[row].value = image_data[row * Nx + col];
        column_info[row].original_pos = row;
    }
    
    /* Sort to match the order of filtered data */
    qsort(column_info, Ny, sizeof(PixelInfo), compare_pixel_info);
    
    /* Restore filtered values to original positions */
    for (sorted_row = 0; sorted_row < Ny; sorted_row++) {
        original_row = column_info[sorted_row].original_pos;
        result_data[original_row * Nx + col] = filtered_data[sorted_row * Nx + col];
    }
}
// Get kernel size from environment variable
int get_kernel_size_from_env() {
    char *env_kernel_size = getenv("KERNEL_SIZE");
    if (env_kernel_size) {
        int size = atoi(env_kernel_size);
        if (size > -1) {
//            printf("Using kernel size %d from environment variable KERNEL_SIZE\n", size);
            return size;
        } else {
            printf("Invalid value (%s) in environment variable KERNEL_SIZE. Using default value\n", env_kernel_size);
        }
    }
    return 5; // Default value
}

// Get number of threads from environment variable
int get_num_threads_from_env() {
    char *env_num_threads = getenv("OMP_NUM_THREADS");
    if (env_num_threads) {
        int threads = atoi(env_num_threads);
        if (threads > 0) {
//            printf("Using %d threads from environment variable OMP_NUM_THREADS\n", threads);
            return threads;
        } else {
            printf("Invalid value (%s) in environment variable OMP_NUM_THREADS. Using default\n", env_num_threads);
        }
    }
    return 40; // Default to 40 threads
}


/* Main OpenMP processing function */
int sort_filter_restore_omp(float *image_data, float *result_data, 
                           unsigned int Nx, unsigned int Ny, int kernel_size, int num_threads) {
    float *sorted_data = NULL;
    float *filtered_data = NULL;
    int i, j;
    
    /* Return if kernel size is 0 or 1 and ensure kernel size is odd */
    if (kernel_size < 2) {
//        kernel_size = 3;
		for (i = 0; i < (int)Ny; i++) {
			for (j = 0; j < (int)Nx; j++) {
				*(result_data+i * Nx + j) = *(image_data+i * Nx + j);
			}
		}
		return 0;
    } else if (kernel_size % 2 == 0) {
        kernel_size++;
        printf("Kernel size adjusted to %d to ensure odd number\n", kernel_size);
    }
    
//    printf("Using kernel size: %d\n", kernel_size);
//    printf("Image size: %dx%d\n", Nx, Ny);
//    printf("Using %d threads\n", num_threads);
    
    /* Set number of threads */
    omp_set_num_threads(num_threads);
    
    /* Allocate memory for intermediate results */
    sorted_data = (float *)malloc(Nx * Ny * sizeof(float));
    filtered_data = (float *)malloc(Nx * Ny * sizeof(float));
    
    if (!sorted_data || !filtered_data) {
        fprintf(stderr, "Failed to allocate memory\n");
        if (sorted_data) free(sorted_data);
        if (filtered_data) free(filtered_data);
        return 1;
    }
    
    /* Process each column in parallel */
//    printf("Sorting columns...\n");
    #pragma omp parallel private(i)
    {
        PixelInfo *column_info;
        column_info = (PixelInfo *)malloc(Ny * sizeof(PixelInfo));
        
        if (column_info) {
            /* Process columns in parallel - C89 style */
            #pragma omp for
            for (i = 0; i < (int)Nx; i++) {
                process_column(image_data, sorted_data, column_info, i, Nx, Ny);
            }
            
            free(column_info);
        }
    }
    
//    printf("Applying median filter...\n");
    /* Apply median filter in parallel - C89 style */
    #pragma omp parallel for private(i, j)
    for (i = 0; i < (int)Ny; i++) {
        for (j = 0; j < (int)Nx; j++) {
            filtered_data[i * Nx + j] = apply_median_filter(sorted_data + i * Nx, Nx, j, kernel_size);
        }
    }
    
//    printf("Restoring to original positions...\n");
    /* Restore to original positions in parallel */
    #pragma omp parallel private(i)
    {
        PixelInfo *column_info;
        column_info = (PixelInfo *)malloc(Ny * sizeof(PixelInfo));
        
        if (column_info) {
            /* Process columns in parallel - C89 style */
            #pragma omp for
            for (i = 0; i < (int)Nx; i++) {
                restore_column(image_data, filtered_data, result_data, column_info, i, Nx, Ny);
            }
            
            free(column_info);
        }
    }
    
    /* Free memory */
    free(sorted_data);
    free(filtered_data);
    
//    printf("OpenMP processing completed successfully\n");
    
    return 0;
}