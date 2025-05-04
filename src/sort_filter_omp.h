#ifndef SORT_FILTER_OMP_H
#define SORT_FILTER_OMP_H

// Function to sort columns, apply median filter, and restore original positions using OpenMP
int sort_filter_restore_omp(float *image_data, float *result_data, 
                           unsigned int Nx, unsigned int Ny, int kernel_size, int num_threads);

int get_kernel_size_from_env();
int get_num_threads_from_env();


#endif // SORT_FILTER_OMP_H
