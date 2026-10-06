#ifndef SORT_FILTER_G_H
#define SORT_FILTER_G_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* GPU (CUDA) version of the sorting-based ring removal.
 * Same signature as sort_filter_restore_omp() so it is a drop-in
 * replacement. num_threads is accepted for compatibility but unused.
 * result_data may be the same buffer as image_data (in-place): the input
 * is copied to the device before the result is copied back. */
int sort_filter_restore_gpu(float *image_data, float *result_data,
                            unsigned int Nx, unsigned int Ny,
                            int kernel_size, int num_threads);

/* Same environment helpers as the OpenMP version, provided here so that
 * ct_rec.c can keep calling them unchanged when built for the GPU. */
int get_kernel_size_from_env(void);
int get_num_threads_from_env(void);

/* Device buffers are allocated on the first call and kept for later calls
 * of the same or smaller size (3 * Nx * Ny * 4 bytes stay allocated).
 * Call this to free them early; they are also released at process exit. */
void sort_filter_gpu_release(void);

/* Pinned (page-locked) host memory for image_data / result_data: makes the
 * host<->device copies inside sort_filter_restore_gpu() run at full PCIe
 * speed.  Falls back to malloc() when cudaMallocHost fails; *pinned tells
 * which one was used and must be passed back to sort_filter_gpu_host_free(). */
void *sort_filter_gpu_host_alloc(size_t bytes, int *pinned);
void  sort_filter_gpu_host_free(void *ptr, int pinned);

#ifdef __cplusplus
}
#endif

#endif /* SORT_FILTER_G_H */
