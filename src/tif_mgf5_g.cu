/*
 * tif_mgf5_g.cu - 2D filter (median + gaussian) for a single TIFF, CUDA
 *                 version specialised for small median kernels (3x3 / 5x5).
 *
 * Drop-in replacement for tif_mgf_g.cu:
 *   tif_mgf5_g <input_file> <output_file> [median_kernel_size] [gaussian_sigma]
 * Same CLI, same edge mode (mirror), same order (median first, then the
 * separable gaussian in float), so the output pixels are identical to
 * tif_mgf_g.cu.  As in tif_mgf_g.cu the gaussian accumulates in single
 * precision, so results may differ from the CPU build (tif_mgf / tif_mgf5)
 * by at most 1 LSB after rounding; the median pass is exact.
 *
 * Why it is faster
 *   1. tif_mgf_g.cu declares "float win[625]" (25x25) in every median
 *      thread.  625 floats = 2500 bytes per thread cannot live in
 *      registers, so the array is placed in local memory (off-chip) and the
 *      insertion sort then does its data-dependent swaps through it.  That
 *      single array dominates the kernel's run time regardless of the
 *      kernel size actually requested.
 *      Here the window size is a template parameter, so 3x3 / 5x5 use a
 *      9- or 25-element array with compile-time indices, which stays in
 *      registers, and the pixel type is kept (no float conversion needed).
 *   2. The sort is replaced by a branch-free sorting network - 19
 *      compare-exchanges for 3x3 (Smith 1996), 99 for 5x5 - which suits
 *      SIMT execution because all threads of a warp run the same
 *      instruction sequence regardless of the data.  The median of a window
 *      is by definition the same value, so the result does not change.
 *   3. Median kernel size 1 is the identity and is now skipped entirely.
 *
 * Median kernel sizes above 5 still work: they use the generic local-array
 * path (the same code tif_mgf_g.cu always uses), so this program can be
 * used as a straight replacement.
 *
 * Even kernel sizes behave as in tif_mgf_g.cu: k=2 -> 3x3, k=4 -> 5x5
 * (the window is (2*(k/2)+1)^2).
 *
 * Compile:
 *   nvcc -O3 -o tif_mgf5_g tif_mgf5_g.cu -ltiff                     (Linux)
 *   nvcc -O3 -o tif_mgf5_g.exe tif_mgf5_g.cu libtiff.lib jpeg.lib lzma.lib zs.lib   (Windows)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <cuda_runtime.h>
#include "cuda13_compat.h"

#ifdef _WIN32
    #define WIN32_LEAN_AND_MEAN
    #include <windows.h>
    #include "tiffio.h"
#else
    #include <stdint.h>
    #include <tiffio.h>
#endif

#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA error at %s:%d - %s\n", __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(1); \
    } \
} while(0)

#define BLOCK_X 16
#define BLOCK_Y 16
#define MEDIAN_MAX 625   /* 25 x 25, generic path only */
#define GAUSS_MAX  64    /* kernel size capped at 51 in calculate_kernel_size() */

/* 1D Gaussian kernel in constant memory (host fills it before the passes). */
__constant__ float c_gk[GAUSS_MAX];

/* ---- device helpers ------------------------------------------------------ */

#define MSORT(a,b) do { if ((b) < (a)) { T t_=(a); (a)=(b); (b)=t_; } } while(0)

/* median of 9 (3x3): 19 compare-exchanges */
template<typename T>
__device__ __forceinline__ T med9(T *p)
{
    MSORT(p[1],p[2]); MSORT(p[4],p[5]); MSORT(p[7],p[8]);
    MSORT(p[0],p[1]); MSORT(p[3],p[4]); MSORT(p[6],p[7]);
    MSORT(p[1],p[2]); MSORT(p[4],p[5]); MSORT(p[7],p[8]);
    MSORT(p[0],p[3]); MSORT(p[5],p[8]); MSORT(p[4],p[7]);
    MSORT(p[3],p[6]); MSORT(p[1],p[4]); MSORT(p[2],p[5]);
    MSORT(p[4],p[7]); MSORT(p[4],p[2]); MSORT(p[6],p[4]);
    MSORT(p[4],p[2]);
    return p[4];
}

/* median of 25 (5x5): 99 compare-exchanges */
template<typename T>
__device__ __forceinline__ T med25(T *p)
{
    MSORT(p[0],p[1]);   MSORT(p[3],p[4]);   MSORT(p[2],p[4]);
    MSORT(p[2],p[3]);   MSORT(p[6],p[7]);   MSORT(p[5],p[7]);
    MSORT(p[5],p[6]);   MSORT(p[9],p[10]);  MSORT(p[8],p[10]);
    MSORT(p[8],p[9]);   MSORT(p[12],p[13]); MSORT(p[11],p[13]);
    MSORT(p[11],p[12]); MSORT(p[15],p[16]); MSORT(p[14],p[16]);
    MSORT(p[14],p[15]); MSORT(p[18],p[19]); MSORT(p[17],p[19]);
    MSORT(p[17],p[18]); MSORT(p[21],p[22]); MSORT(p[20],p[22]);
    MSORT(p[20],p[21]); MSORT(p[23],p[24]); MSORT(p[2],p[5]);
    MSORT(p[3],p[6]);   MSORT(p[0],p[6]);   MSORT(p[0],p[3]);
    MSORT(p[4],p[7]);   MSORT(p[1],p[7]);   MSORT(p[1],p[4]);
    MSORT(p[11],p[14]); MSORT(p[8],p[14]);  MSORT(p[8],p[11]);
    MSORT(p[12],p[15]); MSORT(p[9],p[15]);  MSORT(p[9],p[12]);
    MSORT(p[13],p[16]); MSORT(p[10],p[16]); MSORT(p[10],p[13]);
    MSORT(p[20],p[23]); MSORT(p[17],p[23]); MSORT(p[17],p[20]);
    MSORT(p[21],p[24]); MSORT(p[18],p[24]); MSORT(p[18],p[21]);
    MSORT(p[19],p[22]); MSORT(p[8],p[17]);  MSORT(p[9],p[18]);
    MSORT(p[0],p[18]);  MSORT(p[0],p[9]);   MSORT(p[10],p[19]);
    MSORT(p[1],p[19]);  MSORT(p[1],p[10]);  MSORT(p[11],p[20]);
    MSORT(p[2],p[20]);  MSORT(p[2],p[11]);  MSORT(p[12],p[21]);
    MSORT(p[3],p[21]);  MSORT(p[3],p[12]);  MSORT(p[13],p[22]);
    MSORT(p[4],p[22]);  MSORT(p[4],p[13]);  MSORT(p[14],p[23]);
    MSORT(p[5],p[23]);  MSORT(p[5],p[14]);  MSORT(p[15],p[24]);
    MSORT(p[6],p[24]);  MSORT(p[6],p[15]);  MSORT(p[7],p[16]);
    MSORT(p[7],p[19]);  MSORT(p[13],p[21]); MSORT(p[15],p[23]);
    MSORT(p[7],p[13]);  MSORT(p[7],p[15]);  MSORT(p[1],p[9]);
    MSORT(p[3],p[11]);  MSORT(p[5],p[17]);  MSORT(p[11],p[17]);
    MSORT(p[9],p[17]);  MSORT(p[4],p[10]);  MSORT(p[6],p[12]);
    MSORT(p[7],p[14]);  MSORT(p[4],p[6]);   MSORT(p[4],p[7]);
    MSORT(p[12],p[14]); MSORT(p[10],p[14]); MSORT(p[6],p[7]);
    MSORT(p[10],p[12]); MSORT(p[6],p[10]);  MSORT(p[6],p[17]);
    MSORT(p[12],p[17]); MSORT(p[7],p[17]);  MSORT(p[7],p[10]);
    MSORT(p[12],p[18]); MSORT(p[7],p[12]);  MSORT(p[10],p[18]);
    MSORT(p[12],p[20]); MSORT(p[10],p[20]); MSORT(p[10],p[12]);
    return p[12];
}

__device__ __forceinline__ void insertion_sort(float *arr, int n) {
    for (int i = 1; i < n; i++) {
        float key = arr[i];
        int j = i - 1;
        while (j >= 0 && arr[j] > key) { arr[j + 1] = arr[j]; j--; }
        arr[j + 1] = key;
    }
}

/* mirror reflection of an index, as in tif_mgf.c / tif_mgf_g.cu */
__device__ __forceinline__ int mirror(int v, int n) {
    if (v < 0)  v = -v;
    if (v >= n) v = 2 * n - v - 2;
    if (v < 0)  v = 0;          /* window wider than the image */
    if (v >= n) v = n - 1;
    return v;
}

/* Pick the network for the window size at compile time.  A plain
 * "(K == 3) ? med9(p) : med25(p)" would also compile the 25-element network
 * against a 9-element array, so the choice is made by specialisation. */
template<typename T, int K> struct MedNet;
template<typename T> struct MedNet<T, 3> {
    __device__ __forceinline__ static T of(T *p) { return med9<T>(p); }
};
template<typename T> struct MedNet<T, 5> {
    __device__ __forceinline__ static T of(T *p) { return med25<T>(p); }
};

/* K x K median (K = 3 or 5), mirror boundary.  The window is a fixed-size
 * array indexed by compile-time constants, so it stays in registers. */
template<typename T, int K>
__global__ void median_small_kernel(const T* __restrict__ in, T* __restrict__ out,
                                   int width, int height) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= width || y >= height) return;
    const int half = K / 2;

    T p[K * K];
#pragma unroll
    for (int j = 0; j < K; j++) {
        const int ny = mirror(y + j - half, height);
        const T *row = in + (size_t)ny * width;
#pragma unroll
        for (int i = 0; i < K; i++)
            p[j * K + i] = row[mirror(x + i - half, width)];
    }
    out[(size_t)y * width + x] = MedNet<T, K>::of(p);
}

/* Generic median for kernel sizes above 5 (window in local memory).
 * Same code path as tif_mgf_g.cu, kept so that k > 5 still works. */
template<typename T>
__global__ void median_generic_kernel(const T* __restrict__ in, T* __restrict__ out,
                                      int width, int height, int kernel_size) {
    ENABLE_SMEM_SPILLING();
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= width || y >= height) return;
    const int half = kernel_size / 2;
    float win[MEDIAN_MAX];
    int count = 0;
    for (int j = -half; j <= half; j++) {
        for (int i = -half; i <= half; i++) {
            win[count++] = (float)in[(size_t)mirror(y + j, height) * width
                                     + mirror(x + i, width)];
        }
    }
    insertion_sort(win, count);
    out[(size_t)y * width + x] = (T)win[count / 2];
}

template<typename T>
__global__ void to_float_kernel(const T* __restrict__ in, float* __restrict__ out,
                                int width, int height) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= width || y >= height) return;
    size_t idx = (size_t)y * width + x;
    out[idx] = (float)in[idx];
}

template<typename T>
__global__ void from_float_kernel(const float* __restrict__ in, T* __restrict__ out,
                                  int width, int height, float maxval) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= width || y >= height) return;
    size_t idx = (size_t)y * width + x;
    float v = in[idx] + 0.5f;          /* round, as in tif_mgf.c convert_from_float */
    if (v < 0.0f) v = 0.0f;
    if (v > maxval) v = maxval;
    out[idx] = (T)v;
}

/* Separable Gaussian, horizontal pass, mirror boundary on x. */
__global__ void gauss_h_kernel(const float* __restrict__ in, float* __restrict__ out,
                               int width, int height, int kernel_size) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= width || y >= height) return;
    const int half = kernel_size / 2;
    const float *row = in + (size_t)y * width;
    float s = 0.0f;
    for (int i = 0; i < kernel_size; i++)
        s += row[mirror(x + i - half, width)] * c_gk[i];
    out[(size_t)y * width + x] = s;
}

/* Separable Gaussian, vertical pass, mirror boundary on y. */
__global__ void gauss_v_kernel(const float* __restrict__ in, float* __restrict__ out,
                               int width, int height, int kernel_size) {
    const int x = blockIdx.x * blockDim.x + threadIdx.x;
    const int y = blockIdx.y * blockDim.y + threadIdx.y;
    if (x >= width || y >= height) return;
    const int half = kernel_size / 2;
    float s = 0.0f;
    for (int i = 0; i < kernel_size; i++)
        s += in[(size_t)mirror(y + i - half, height) * width + x] * c_gk[i];
    out[(size_t)y * width + x] = s;
}

/* ---- host helpers (identical math to tif_mgf.c) -------------------------- */

static int calculate_kernel_size(double sigma) {
    int size;
    if (sigma <= 0.0) return 0;
    size = (int)ceil(3.0 * sigma) * 2 + 1;
    if (size % 2 == 0) size++;
    if (size > 51) size = 51;
    return size;
}

static void create_gaussian_kernel_1d(double *kernel, int size, double sigma) {
    int center = size / 2;
    double sum = 0.0;
    double two_sigma_sq = 2.0 * sigma * sigma;
    double norm = 1.0 / sqrt(2.0 * 3.14159265358979323846 * sigma * sigma);
    for (int i = 0; i < size; i++) {
        int d = i - center;
        kernel[i] = norm * exp(-(double)(d * d) / two_sigma_sq);
        sum += kernel[i];
    }
    for (int i = 0; i < size; i++) kernel[i] /= sum;
}

/* Copy non-geometry metadata from input to output (as in tif_mgf_g.cu).
 * Guards every field, so a missing tag is skipped rather than propagating a
 * NULL/garbage value.  Call before closing the input TIFF. */
static void copy_tiff_metadata(TIFF *in, TIFF *out) {
    uint32_t u32, count;
    uint16_t u16, u16a, u16b, *u16ptr;
    float f;
    double dbl;
    char *str;
    void *data;
    if (TIFFGetField(in, TIFFTAG_IMAGEDESCRIPTION, &str)) TIFFSetField(out, TIFFTAG_IMAGEDESCRIPTION, str);
    if (TIFFGetField(in, TIFFTAG_MAKE,             &str)) TIFFSetField(out, TIFFTAG_MAKE,             str);
    if (TIFFGetField(in, TIFFTAG_MODEL,            &str)) TIFFSetField(out, TIFFTAG_MODEL,            str);
    if (TIFFGetField(in, TIFFTAG_SOFTWARE,         &str)) TIFFSetField(out, TIFFTAG_SOFTWARE,         str);
    if (TIFFGetField(in, TIFFTAG_DATETIME,         &str)) TIFFSetField(out, TIFFTAG_DATETIME,         str);
    if (TIFFGetField(in, TIFFTAG_ARTIST,           &str)) TIFFSetField(out, TIFFTAG_ARTIST,           str);
    if (TIFFGetField(in, TIFFTAG_HOSTCOMPUTER,     &str)) TIFFSetField(out, TIFFTAG_HOSTCOMPUTER,     str);
    if (TIFFGetField(in, TIFFTAG_COPYRIGHT,        &str)) TIFFSetField(out, TIFFTAG_COPYRIGHT,        str);
    if (TIFFGetField(in, TIFFTAG_DOCUMENTNAME,     &str)) TIFFSetField(out, TIFFTAG_DOCUMENTNAME,     str);
    if (TIFFGetField(in, TIFFTAG_PAGENAME,         &str)) TIFFSetField(out, TIFFTAG_PAGENAME,         str);
    if (TIFFGetField(in, TIFFTAG_TARGETPRINTER,    &str)) TIFFSetField(out, TIFFTAG_TARGETPRINTER,    str);
    if (TIFFGetField(in, TIFFTAG_SUBFILETYPE, &u32)) TIFFSetField(out, TIFFTAG_SUBFILETYPE, u32);
    if (TIFFGetField(in, TIFFTAG_ORIENTATION,      &u16)) TIFFSetField(out, TIFFTAG_ORIENTATION,      u16);
    if (TIFFGetField(in, TIFFTAG_RESOLUTIONUNIT,   &u16)) TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT,   u16);
    if (TIFFGetField(in, TIFFTAG_XRESOLUTION, &f)) TIFFSetField(out, TIFFTAG_XRESOLUTION, f);
    if (TIFFGetField(in, TIFFTAG_YRESOLUTION, &f)) TIFFSetField(out, TIFFTAG_YRESOLUTION, f);
    if (TIFFGetField(in, TIFFTAG_XPOSITION,   &f)) TIFFSetField(out, TIFFTAG_XPOSITION,   f);
    if (TIFFGetField(in, TIFFTAG_YPOSITION,   &f)) TIFFSetField(out, TIFFTAG_YPOSITION,   f);
    if (TIFFGetField(in, TIFFTAG_SMINSAMPLEVALUE, &dbl)) TIFFSetField(out, TIFFTAG_SMINSAMPLEVALUE, dbl);
    if (TIFFGetField(in, TIFFTAG_SMAXSAMPLEVALUE, &dbl)) TIFFSetField(out, TIFFTAG_SMAXSAMPLEVALUE, dbl);
    if (TIFFGetField(in, TIFFTAG_STONITS,         &dbl)) TIFFSetField(out, TIFFTAG_STONITS,         dbl);
    if (TIFFGetField(in, TIFFTAG_PAGENUMBER, &u16a, &u16b)) TIFFSetField(out, TIFFTAG_PAGENUMBER, u16a, u16b);
    if (TIFFGetField(in, TIFFTAG_ICCPROFILE, &count, &data)) TIFFSetField(out, TIFFTAG_ICCPROFILE, count, data);
    if (TIFFGetField(in, TIFFTAG_EXTRASAMPLES, &u16, &u16ptr)) TIFFSetField(out, TIFFTAG_EXTRASAMPLES, u16, u16ptr);
}

/* ---- GPU pipeline -------------------------------------------------------- */

/* Runs median (optional) then gaussian (optional) on the host buffer in place.
 * T is the pixel type (unsigned char / unsigned short). */
template<typename T>
static void run_pipeline(void *host, int width, int height,
                         int median_size, double sigma, float maxval) {
    const size_t npix = (size_t)width * height;
    const size_t nbytes = npix * sizeof(T);
    dim3 block(BLOCK_X, BLOCK_Y);
    dim3 grid((width + block.x - 1) / block.x, (height + block.y - 1) / block.y);

    T *d_a = NULL, *d_b = NULL;
    CUDA_CHECK(cudaMalloc(&d_a, nbytes));
    CUDA_CHECK(cudaMalloc(&d_b, nbytes));
    CUDA_CHECK(cudaMemcpy(d_a, host, nbytes, cudaMemcpyHostToDevice));

    /* Median (integer, exact).  half = median_size/2, so 2 -> 3x3 and
     * 4 -> 5x5, and 1 (half = 0) is the identity and is skipped. */
    const int half = median_size / 2;
    if (half > 0) {
        if (half == 1)
            median_small_kernel<T, 3><<<grid, block>>>(d_a, d_b, width, height);
        else if (half == 2)
            median_small_kernel<T, 5><<<grid, block>>>(d_a, d_b, width, height);
        else
            median_generic_kernel<T><<<grid, block>>>(d_a, d_b, width, height,
                                                      2 * half + 1);
        CUDA_CHECK(cudaGetLastError());
        T *tmp = d_a; d_a = d_b; d_b = tmp;   /* result now in d_a */
    }

    /* Separable Gaussian (float) */
    if (sigma > 0.0) {
        int ks = calculate_kernel_size(sigma);
        if (ks > 0) {
            double kd[GAUSS_MAX];
            float  kf[GAUSS_MAX];
            create_gaussian_kernel_1d(kd, ks, sigma);
            for (int i = 0; i < ks; i++) kf[i] = (float)kd[i];
            CUDA_CHECK(cudaMemcpyToSymbol(c_gk, kf, ks * sizeof(float)));

            float *d_f = NULL, *d_ft = NULL;
            CUDA_CHECK(cudaMalloc(&d_f,  npix * sizeof(float)));
            CUDA_CHECK(cudaMalloc(&d_ft, npix * sizeof(float)));
            to_float_kernel<T><<<grid, block>>>(d_a, d_f, width, height);
            CUDA_CHECK(cudaGetLastError());
            gauss_h_kernel<<<grid, block>>>(d_f, d_ft, width, height, ks);   /* horizontal */
            CUDA_CHECK(cudaGetLastError());
            gauss_v_kernel<<<grid, block>>>(d_ft, d_f, width, height, ks);   /* vertical */
            CUDA_CHECK(cudaGetLastError());
            from_float_kernel<T><<<grid, block>>>(d_f, d_a, width, height, maxval);
            CUDA_CHECK(cudaGetLastError());
            CUDA_CHECK(cudaFree(d_f));
            CUDA_CHECK(cudaFree(d_ft));
        }
    }

    CUDA_CHECK(cudaDeviceSynchronize());
    CUDA_CHECK(cudaMemcpy(host, d_a, nbytes, cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_a));
    CUDA_CHECK(cudaFree(d_b));
}

static int process_tiff_file(const char *input_file, const char *output_file,
                             int median_size, double sigma) {
    TIFF *in_tiff = TIFFOpen(input_file, "r");
    if (!in_tiff) { fprintf(stderr, "Error: Cannot open input TIFF file\n"); return 1; }

    uint32_t width = 0, height = 0;
    uint16_t bits_per_sample = 0;
    TIFFGetField(in_tiff, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(in_tiff, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(in_tiff, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);

    if (bits_per_sample != 8 && bits_per_sample != 16) {
        fprintf(stderr, "Error: Only 8-bit and 16-bit images are supported\n");
        TIFFClose(in_tiff);
        return 1;
    }

    tsize_t scanline_size = TIFFScanlineSize(in_tiff);
    void *data = malloc((size_t)height * scanline_size);
    if (!data) { fprintf(stderr, "Error: Cannot allocate memory\n"); TIFFClose(in_tiff); return 1; }

    for (uint32_t row = 0; row < height; row++) {
        if (TIFFReadScanline(in_tiff, (char*)data + (size_t)row * scanline_size, row, 0) < 0) {
            fprintf(stderr, "Error: Cannot read scanline %u\n", row);
            free(data); TIFFClose(in_tiff); return 1;
        }
    }

    /* GPU processing */
    if (bits_per_sample == 16)
        run_pipeline<unsigned short>(data, width, height, median_size, sigma, 65535.0f);
    else
        run_pipeline<unsigned char>(data, width, height, median_size, sigma, 255.0f);

    /* Write output (same tags as tif_mgf_g.cu) */
    TIFF *out_tiff = TIFFOpen(output_file, "w");
    if (!out_tiff) {
        fprintf(stderr, "Error: Cannot create output TIFF file\n");
        free(data); TIFFClose(in_tiff); return 1;
    }
    copy_tiff_metadata(in_tiff, out_tiff);
    TIFFClose(in_tiff);

    TIFFSetField(out_tiff, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField(out_tiff, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField(out_tiff, TIFFTAG_BITSPERSAMPLE, bits_per_sample);
    TIFFSetField(out_tiff, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(out_tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(out_tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(out_tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out_tiff, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize(out_tiff, 0));
    TIFFSetField(out_tiff, TIFFTAG_ARTIST, "tif_mgf5_g_libtiff");

    for (uint32_t row = 0; row < height; row++) {
        if (TIFFWriteScanline(out_tiff, (char*)data + (size_t)row * scanline_size, row, 0) < 0) {
            fprintf(stderr, "Error: Cannot write scanline %u\n", row);
            free(data); TIFFClose(out_tiff); return 1;
        }
    }
    TIFFClose(out_tiff);
    free(data);
    return 0;
}

int main(int argc, char *argv[]) {
    if (argc < 3 || argc > 5) {
        fprintf(stderr, "Usage: %s <input_file> <output_file> [median_kernel_size] [gaussian_sigma]\n", argv[0]);
        fprintf(stderr, "  median_kernel_size: 0-25 (default: 1, 0 to skip; 3 and 5 take the fast path)\n");
        fprintf(stderr, "  gaussian_sigma: 0.0-100.0 (default: 0.5, 0.0 to skip)\n");
        return 1;
    }

    const char *input_file  = argv[1];
    const char *output_file = argv[2];
    int    median_kernel_size = 1;
    double gaussian_sigma     = 0.5;

    if (argc > 3) {
        median_kernel_size = atoi(argv[3]);
        if (median_kernel_size < 0 || median_kernel_size > 25) {
            fprintf(stderr, "Error: Invalid median kernel size (must be 0-25)\n");
            return 1;
        }
    }
    if (argc > 4) {
        gaussian_sigma = atof(argv[4]);
        if (gaussian_sigma < 0.0 || gaussian_sigma > 100.0) {
            fprintf(stderr, "Error: Invalid gaussian sigma (must be 0.0-100.0)\n");
            return 1;
        }
    }

    int device_count = 0;
    CUDA_CHECK(cudaGetDeviceCount(&device_count));
    if (device_count == 0) { fprintf(stderr, "Error: No CUDA-capable devices found\n"); return 1; }
    CUDA_CHECK(cudaSetDevice(cuda_select_best_gpu()));

    int result = process_tiff_file(input_file, output_file, median_kernel_size, gaussian_sigma);

    /* append to log file (same format as tif_mgf.c) */
    {
        FILE *flog = fopen("cmd-hst.log", "a");
        if (flog) {
            for (int i = 0; i < argc; ++i) fprintf(flog, "%s ", argv[i]);
            fprintf(flog, "\t");
            fprintf(flog, "   %% median_kernel_size %d  gaussian_sigma %g\n",
                    median_kernel_size, gaussian_sigma);
            fclose(flog);
        }
    }

    if (result == 0) CUDA_CHECK(cudaDeviceReset());
    return result;
}
