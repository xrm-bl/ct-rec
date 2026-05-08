/*
 * cuda13_compat.h - CUDA 12.2 / 13.x Compatibility Layer
 *
 * Provides:
 *   1. ENABLE_SMEM_SPILLING  - opt-in shared memory register spilling (CUDA 13.0+)
 *   2. safe_prefetch()       - cudaMemPrefetchAsync wrapper for API change (CUDA 12.2+)
 *   3. Deprecation-safe device query helper
 *
 * Usage: #include "cuda13_compat.h" at the top of each .cu file
 */

#ifndef CUDA13_COMPAT_H
#define CUDA13_COMPAT_H

#include <cuda_runtime.h>

/*
 * Shared-memory register spilling (CUDA 13.0+ / PTXAS 9.0+)
 *
 * Insert ENABLE_SMEM_SPILLING(); as the first statement inside any
 * __global__ kernel that has high register pressure.  On older
 * toolkits the macro compiles to nothing.
 *
 * NOTE: The kernel must have sufficient unused shared memory for this
 * to be beneficial.  If the kernel already uses most of the SM's shared
 * memory budget, enabling spilling may *reduce* occupancy.
 */
#if CUDART_VERSION >= 13000
  #define ENABLE_SMEM_SPILLING()  asm volatile(".pragma \"enable_smem_spilling\";")
#else
  #define ENABLE_SMEM_SPILLING()  ((void)0)
#endif

/*
 * cudaMemPrefetchAsync compatibility wrapper.
 *
 * API evolution:
 *   CUDA <= 12.1 : (ptr, count, int device,    cudaStream_t)        4 args
 *   CUDA 12.2+   : (ptr, count, cudaMemLocation, cudaStream_t)      4 args (device->location)
 *   CUDA 13.0+   : (ptr, count, cudaMemLocation, unsigned flags, cudaStream_t)  5 args
 */
static inline cudaError_t safe_prefetch(const void *devPtr, size_t count,
                                        int dstDevice, cudaStream_t stream)
{
#if CUDART_VERSION >= 13000
    /* CUDA 13.0+: 5-argument form with flags */
    cudaMemLocation loc;
    loc.type = cudaMemLocationTypeDevice;
    loc.id   = dstDevice;
    return cudaMemPrefetchAsync(devPtr, count, loc, 0, stream);
#elif CUDART_VERSION >= 12020
    /* CUDA 12.2-12.x: 4-argument form with cudaMemLocation */
    cudaMemLocation loc;
    loc.type = cudaMemLocationTypeDevice;
    loc.id   = dstDevice;
    return cudaMemPrefetchAsync(devPtr, count, loc, stream);
#else
    /* CUDA <= 12.1: 4-argument form with int device */
    return cudaMemPrefetchAsync(devPtr, count, dstDevice, stream);
#endif
}

/*
 * Best-GPU selection helper.
 *
 * cudaDeviceProp::multiProcessorCount is still valid in CUDA 13.x,
 * but several other members were removed.  This helper isolates
 * the query so that callers need not touch cudaDeviceProp directly.
 */
static inline int cuda_select_best_gpu(void)
{
    int device_count = 0, best = 0, max_sm = 0;
    cudaGetDeviceCount(&device_count);
    for (int i = 0; i < device_count; i++) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, i);
        if (prop.multiProcessorCount > max_sm) {
            max_sm = prop.multiProcessorCount;
            best = i;
        }
    }
    return best;
}

#endif /* CUDA13_COMPAT_H */