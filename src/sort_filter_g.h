#ifndef SORT_FILTER_G_H
#define SORT_FILTER_G_H

#ifdef __cplusplus
extern "C" {
#endif

/* GPU (CUDA) version of the sorting-based ring removal.
 * Same signature as sort_filter_restore_omp() so it is a drop-in
 * replacement. num_threads is accepted for compatibility but unused. */
int sort_filter_restore_gpu(float *image_data, float *result_data,
                            unsigned int Nx, unsigned int Ny,
                            int kernel_size, int num_threads);

/* Same environment helpers as the OpenMP version, provided here so that
 * ct_rec.c can keep calling them unchanged when built for the GPU. */
int get_kernel_size_from_env(void);
int get_num_threads_from_env(void);

#ifdef __cplusplus
}
#endif

#endif /* SORT_FILTER_G_H */
