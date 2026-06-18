/*----------------------------------------------------------------------*/
/* sort_filter_g.cu  (Thrust-free version)                               */
/*                                                                       */
/* CUDA implementation of the sorting-based ring (stripe) removal,       */
/* Vo et al., Opt. Express 26, 28396-28412 (2018).                       */
/*                                                                       */
/* Drop-in GPU replacement for sort_filter_omp.c, identical public       */
/* interface to sort_filter_restore_omp().                               */
/*                                                                       */
/* This version uses NO Thrust / CUB (CCCL). The per-column sort is done */
/* with a self-contained bitonic sort: one thread-block per column,      */
/* sorting that column in shared memory while carrying the original row  */
/* index. Ties are broken by the original index, so the result is a      */
/* deterministic, reproducible total order.                              */
/*                                                                       */
/* Because there is no C++ library dependency, it compiles with the      */
/* plain nvcc flags (no /Zc:preprocessor) and links without -lstdc++.    */
/*----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include <cuda_runtime.h>

#include "sort_filter_g.h"

/* Largest supported median kernel */
#define MAX_K 129

/* Threads per block for the per-column sort */
#define SORT_THREADS 256

/* Static shared-memory limit (bytes) below which no opt-in is needed */
#define SMEM_STATIC_LIMIT 49152   /* 48 KB */

/*----------------------------------------------------------------------*/
#define CUDA_CHECK(call)                                                   \
    do {                                                                   \
        cudaError_t err__ = (call);                                        \
        if (err__ != cudaSuccess) {                                        \
            fprintf(stderr, "CUDA error %s:%d: %s\n",                      \
                    __FILE__, __LINE__, cudaGetErrorString(err__));        \
            return 1;                                                      \
        }                                                                  \
    } while (0)

/*----------------------------------------------------------------------*/
/* Smallest power of two >= v (and >= 2). */
static unsigned int next_pow2(unsigned int v)
{
    unsigned int p = 1;
    while (p < v) p <<= 1;
    return p < 2 ? 2 : p;
}

/*----------------------------------------------------------------------*/
/* Per-column bitonic sort.  One block per column (blockIdx.x == column).*/
/* Shared memory holds pad floats (values) followed by pad ints (rows).  */
/* pad is a power of two >= Ny; positions [Ny,pad) are +FLT_MAX padding  */
/* that sorts to the top and is never stored back.                       */
__global__ void sort_columns_bitonic(const float *img,
                                     float *sval, unsigned int *perm,
                                     unsigned int Nx, unsigned int Ny,
                                     unsigned int pad)
{
    extern __shared__ unsigned char smem[];
    float *sv = (float *)smem;
    int   *si = (int   *)(sv + pad);

    unsigned int c = blockIdx.x;
    if (c >= Nx) return;

    /* load column c into shared memory (with padding) */
    for (unsigned int i = threadIdx.x; i < pad; i += blockDim.x) {
        if (i < Ny) {
            sv[i] = img[(size_t)i * Nx + c];
            si[i] = (int)i;
        } else {
            sv[i] = FLT_MAX;     /* padding sorts to the top */
            si[i] = -1;
        }
    }
    __syncthreads();

    /* bitonic sort network; total order = (value, then original index) */
    for (unsigned int k = 2; k <= pad; k <<= 1) {
        for (unsigned int j = k >> 1; j > 0; j >>= 1) {
            for (unsigned int i = threadIdx.x; i < pad; i += blockDim.x) {
                unsigned int ixj = i ^ j;
                if (ixj > i) {
                    int ascending = ((i & k) == 0);
                    #define GT(x, y) ( sv[x] >  sv[y] || \
                                      (sv[x] == sv[y] && si[x] > si[y]) )
                    int sw = ascending ? GT(i, ixj) : GT(ixj, i);
                    #undef GT
                    if (sw) {
                        float tf = sv[i];   sv[i]   = sv[ixj]; sv[ixj] = tf;
                        int   ti = si[i];   si[i]   = si[ixj]; si[ixj] = ti;
                    }
                }
            }
            __syncthreads();
        }
    }

    /* store the real ranks [0,Ny) in column-major layout */
    for (unsigned int i = threadIdx.x; i < Ny; i += blockDim.x) {
        sval[(size_t)c * Ny + i] = sv[i];
        perm[(size_t)c * Ny + i] = (unsigned int)si[i];
    }
}

/*----------------------------------------------------------------------*/
/* Median filter (horizontal, across columns at fixed rank) fused with   */
/* the scatter back to original row positions.                           */
/* Boundary handling mirrors the CPU apply_median_filter() exactly.      */
__global__ void median_scatter(const float *sval,
                               const unsigned int *perm,
                               float *result,
                               unsigned int Nx, unsigned int Ny, int ks)
{
    size_t e     = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    size_t total = (size_t)Nx * (size_t)Ny;
    if (e >= total) return;

    unsigned int c = (unsigned int)(e / Ny);   /* column */
    unsigned int i = (unsigned int)(e % Ny);   /* rank   */

    int   half = ks / 2;
    float w[MAX_K];
    int   cnt = 0;
    int   k;

    for (k = -half; k <= half; k++) {
        long pos = (long)c + k;
        if (pos < 0) {
            pos = -pos;
        } else if (pos >= (long)Nx) {
            pos = 2L * (long)Nx - pos - 1;
        }
        if (pos < 0)         pos = 0;
        if (pos >= (long)Nx) pos = (long)Nx - 1;

        w[cnt++] = sval[(size_t)pos * Ny + i];
    }

    /* insertion sort of the (small) kernel window */
    for (k = 1; k < ks; k++) {
        float key = w[k];
        int   b   = k - 1;
        while (b >= 0 && w[b] > key) { w[b + 1] = w[b]; b--; }
        w[b + 1] = key;
    }

    float med = w[ks / 2];

    result[(size_t)perm[e] * Nx + c] = med;
}

/*----------------------------------------------------------------------*/
/* Environment-variable helpers (same behaviour as the CPU version).     */
extern "C" int get_kernel_size_from_env(void)
{
    char *env_kernel_size = getenv("KERNEL_SIZE");
    if (env_kernel_size) {
        int size = atoi(env_kernel_size);
        if (size > -1) {
            return size;
        } else {
            printf("Invalid value (%s) in environment variable KERNEL_SIZE. "
                   "Using default value\n", env_kernel_size);
        }
    }
    return 5;
}

extern "C" int get_num_threads_from_env(void)
{
    char *env_num_threads = getenv("OMP_NUM_THREADS");
    if (env_num_threads) {
        int threads = atoi(env_num_threads);
        if (threads > 0) {
            return threads;
        } else {
            printf("Invalid value (%s) in environment variable OMP_NUM_THREADS. "
                   "Using default\n", env_num_threads);
        }
    }
    return 40; /* kept only for interface compatibility (unused on GPU) */
}

/*----------------------------------------------------------------------*/
/* Main GPU processing function.  Same signature as the OpenMP version;  */
/* num_threads is accepted for compatibility but unused on the GPU.      */
extern "C" int sort_filter_restore_gpu(float *image_data, float *result_data,
                                       unsigned int Nx, unsigned int Ny,
                                       int kernel_size, int num_threads)
{
    (void)num_threads;

    size_t total = (size_t)Nx * (size_t)Ny;

    if (kernel_size < 2) {
        memcpy(result_data, image_data, total * sizeof(float));
        return 0;
    }
    if (kernel_size % 2 == 0) {
        kernel_size++;
        printf("Kernel size adjusted to %d to ensure odd number\n", kernel_size);
    }
    if (kernel_size > MAX_K - 1) {
        fprintf(stderr, "kernel_size %d too large (max %d)\n",
                kernel_size, MAX_K - 1);
        return 1;
    }

    unsigned int pad     = next_pow2(Ny);
    size_t       shbytes = (size_t)pad * (sizeof(float) + sizeof(int));

    /* Make sure the column fits in shared memory (opt-in if > 48 KB). */
    if (shbytes > SMEM_STATIC_LIMIT) {
        int maxOptin = 0;
        CUDA_CHECK(cudaDeviceGetAttribute(
            &maxOptin, cudaDevAttrMaxSharedMemoryPerBlockOptin, 0));
        if ((size_t)maxOptin < shbytes) {
            fprintf(stderr,
                "Ny=%u needs %zu bytes of shared memory per block, but the "
                "device allows at most %d. Reduce Ny or use the Thrust-based "
                "build for very tall sinograms.\n", Ny, shbytes, maxOptin);
            return 1;
        }
        CUDA_CHECK(cudaFuncSetAttribute(
            sort_columns_bitonic,
            cudaFuncAttributeMaxDynamicSharedMemorySize, (int)shbytes));
    }

    float        *d_image = NULL, *d_result = NULL, *d_sval = NULL;
    unsigned int *d_perm  = NULL;

    CUDA_CHECK(cudaMalloc((void **)&d_image,  total * sizeof(float)));
    CUDA_CHECK(cudaMalloc((void **)&d_sval,   total * sizeof(float)));
    CUDA_CHECK(cudaMalloc((void **)&d_perm,   total * sizeof(unsigned int)));

    CUDA_CHECK(cudaMemcpy(d_image, image_data, total * sizeof(float),
                          cudaMemcpyHostToDevice));

    /* 1) sort every column (one block per column) */
    sort_columns_bitonic<<<Nx, SORT_THREADS, shbytes>>>(
        d_image, d_sval, d_perm, Nx, Ny, pad);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaFree(d_image));
    d_image = NULL;

    /* 2) median filter + scatter back to original positions */
    CUDA_CHECK(cudaMalloc((void **)&d_result, total * sizeof(float)));

    int block = 256;
    int grid  = (int)((total + block - 1) / block);
    median_scatter<<<grid, block>>>(d_sval, d_perm, d_result, Nx, Ny, kernel_size);
    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaMemcpy(result_data, d_result, total * sizeof(float),
                          cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaFree(d_sval));
    CUDA_CHECK(cudaFree(d_perm));
    CUDA_CHECK(cudaFree(d_result));

    return 0;
}
