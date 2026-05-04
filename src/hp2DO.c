/*
 * hp2DO.c - Rotation axis determination from 0/180 degree projections
 *
 * Supports both img (SPring-8 HIS format) and 16bit monochrome TIFF input.
 * The input format is auto-detected from the extension of the first argument:
 *   *.tif  -> TIFF mode
 *   otherwise -> img mode
 *
 * Usage:
 *   hp2DO D000 I000 E000 D180 I180 E180 {y1 y2} {outputImage} {search%}
 *
 * Compile (Linux):
 *   gcc -O2 -fopenmp -o hp2DO hp2DO.c -ltiff -lm
 *
 * Compile (Windows, MSVC):
 *   cl /O2 /openmp hp2DO.c libtiff.lib /Fehp2DO.exe
 */

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>
#include "tiffio.h"

/* ---------------------------------------------------------------
 * Common utilities
 * --------------------------------------------------------------- */

static void Abort(const char *path, const char *msg)
{
    fprintf(stderr, "%s : %s\n", path, msg);
    exit(1);
}

static void Error(const char *msg)
{
    fputs(msg, stderr);
    fputc('\n', stderr);
    exit(1);
}

#define MA(cnt, ptr) malloc((size_t)(cnt) * sizeof(*(ptr)))

static double Log_val(double I0, double I)
{
    return (I0 <= 0.0 || I <= 0.0) ? 0.0 : log(I0 / I);
}

static double Sqrt_val(double d)
{
    return (d < 0.0) ? 0.0 : sqrt(d);
}

static double CalculateRMS(int N1, int dx, int x1, int x2,
                           double *p000, double *p180)
{
    int x;
    double d, sum = 0.0;

    for (x = x1; x < x2; x++) {
        d = p000[x + dx] - p180[N1 - x];
        sum += (d * d);
    }
    return Sqrt_val(sum / (double)(x2 - x1));
}

/* Check if filename ends with .tif/.tiff */
static int IsTiffFile(const char *path)
{
    const char *dot = strrchr(path, '.');
    if (dot == NULL) return 0;
    if (strcmp(dot, ".tif") == 0 || strcmp(dot, ".TIF") == 0 ||
        strcmp(dot, ".tiff") == 0 || strcmp(dot, ".TIFF") == 0)
        return 1;
    return 0;
}

/* ---------------------------------------------------------------
 * IMG (HIS) format: read entire file into unsigned short **
 * Applies 1D Gaussian filter (radius=1, sigma=1.0) to each scanline.
 * --------------------------------------------------------------- */

#define GetByte(file) fgetc(file)

static int GetWord(FILE *file)
{
    int lo, hi;
    return ((lo = GetByte(file)) == EOF ||
            (hi = GetByte(file)) == EOF) ? -1 : (hi << 8) | lo;
}

static unsigned short **ReadImgFile(const char *path, int *Nx, int *Ny,
                                    int firstFile)
{
    FILE *file;
    int len, x, y, word;
    unsigned short **data;

    file = fopen(path, "rb");
    if (file == NULL) Abort(path, "file not open.");

    if (GetByte(file) != 'I' ||
        GetByte(file) != 'M') Abort(path, "bad magic number.");

    if ((len = GetWord(file)) < 0) Abort(path, "bad comment length.");

    if (firstFile) {
        if ((*Nx = GetWord(file)) <= 0 ||
            (*Ny = GetWord(file)) <= 0) Abort(path, "bad image size.");
    } else {
        if (*Nx != GetWord(file) ||
            *Ny != GetWord(file)) Abort(path, "image size not match.");
    }

    fseek(file, 4L, 1);
    if (GetWord(file) != 2) Abort(path, "bad image type.");
    fseek(file, (long)(len + 50), 1);

    /* Allocate image memory */
    data = (unsigned short **)MA(*Ny, data);
    if (data == NULL) Abort(path, "memory allocation error.");
    for (y = 0; y < *Ny; y++) {
        data[y] = (unsigned short *)MA(*Nx, *data);
        if (data[y] == NULL) Abort(path, "memory allocation error.");
    }

    /* Read all scanlines */
    for (y = 0; y < *Ny; y++) {
        for (x = 0; x < *Nx; x++) {
            if ((word = GetWord(file)) < 0)
                Abort(path, "unexpected end of file.");
            data[y][x] = (unsigned short)word;
        }
    }

    fclose(file);
    return data;
}

/* ---------------------------------------------------------------
 * TIFF format: read entire file into unsigned short **
 * --------------------------------------------------------------- */

static unsigned short **ReadTiffFile(const char *path, int *Nx, int *Ny,
                                     int firstFile)
{
    TIFF *tif;
    uint32_t width, height;
    uint16_t bps, spp;
    unsigned short **data;
    int y;

    tif = TIFFOpen(path, "r");
    if (tif == NULL) Abort(path, "file not open.");

    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);

    if (bps != 16) Abort(path, "not 16bit image.");
    if (spp != 1) Abort(path, "not monochrome image.");

    if (firstFile) {
        *Nx = (int)width;
        *Ny = (int)height;
    } else {
        if (*Nx != (int)width || *Ny != (int)height)
            Abort(path, "image size not match.");
    }

    data = (unsigned short **)MA(*Ny, data);
    if (data == NULL) Abort(path, "memory allocation error.");

    for (y = 0; y < *Ny; y++) {
        data[y] = (unsigned short *)MA(*Nx, *data);
        if (data[y] == NULL) Abort(path, "memory allocation error.");
    }

    for (y = 0; y < (int)height; y++) {
        if (TIFFReadScanline(tif, data[y], y, 0) < 0)
            Abort(path, "error reading scanline.");
    }

    TIFFClose(tif);
    return data;
}

/* ---------------------------------------------------------------
 * Common: free 2D image, calculate projection, store TIFF output
 * --------------------------------------------------------------- */

static void FreeImageData(unsigned short **data, int Ny)
{
    int y;
    for (y = 0; y < Ny; y++) free(data[y]);
    free(data);
}

/*
 * Apply 1D Gaussian filter (radius=1, sigma=1.0) to a scanline.
 * src: input unsigned short array, dst: output double array, Nx: length.
 */
static void GaussianFilter1D(unsigned short *src, double *dst, int Nx)
{
    static const double K0 = 0.6065306597633104;  /* exp(-0.5) */
    static const double K1 = 1.0;                 /* exp(0)    */
    int x, i;

    for (x = 0; x < Nx; x++) {
        double filtered = 0.0;
        double total_weight = 0.0;
        for (i = -1; i <= 1; i++) {
            int idx = x + i;
            if (idx >= 0 && idx < Nx) {
                double w = (i == 0) ? K1 : K0;
                filtered += w * (double)src[idx];
                total_weight += w;
            }
        }
        dst[x] = filtered / total_weight;
    }
}

/*
 * Calculate projection: apply Gaussian filter to D, I, E scanlines,
 * then compute p[x] = -log( (filtE - filtD) / (filtI - filtD) ).
 * fD, fI, fE are pre-allocated work buffers of length Nx.
 */
static void CalcProjection(unsigned short **dataD, unsigned short **dataI,
                           unsigned short **dataE, int Nx, int y,
                           double *p, double *fD, double *fI, double *fE)
{
    int x;

    GaussianFilter1D(dataD[y], fD, Nx);
    GaussianFilter1D(dataI[y], fI, Nx);
    GaussianFilter1D(dataE[y], fE, Nx);

    for (x = 0; x < Nx; x++) {
        double I0 = fI[x] - fD[x];
        double I  = fE[x] - fD[x];
        p[x] = Log_val(I0, I);
    }
}

static void StoreTiffFile(const char *path, int width, int height,
                          unsigned short **cell, const char *desc)
{
    TIFF *tif;
    int y;
    uint16_t *buf;

    tif = TIFFOpen(path, "w");
    if (tif == NULL) {
        fprintf(stderr, "%s : cannot create output file.\n", path);
        return;
    }

    TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, (uint32_t)width);
    TIFFSetField(tif, TIFFTAG_IMAGELENGTH, (uint32_t)height);
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, 1);
    TIFFSetField(tif, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    if (desc != NULL && strlen(desc) > 0)
        TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, desc);

    buf = (uint16_t *)_TIFFmalloc(width * sizeof(uint16_t));
    if (buf == NULL) {
        fprintf(stderr, "%s : memory allocation error.\n", path);
        TIFFClose(tif);
        return;
    }

    for (y = 0; y < height; y++) {
        memcpy(buf, cell[y], width * sizeof(uint16_t));
        TIFFWriteScanline(tif, buf, y, 0);
    }

    _TIFFfree(buf);
    TIFFClose(tif);
}

/* ---------------------------------------------------------------
 * Main
 * --------------------------------------------------------------- */

#define PathD000 argv[1]
#define PathI000 argv[2]
#define PathE000 argv[3]
#define PathD180 argv[4]
#define PathI180 argv[5]
#define PathE180 argv[6]

#define BPS      16
#define DESC_LEN 256

#define DEFAULT_SEARCH_PERCENT 20.0

int main(int argc, char **argv)
{
    unsigned short **dataD000, **dataI000, **dataE000,
                   **dataD180, **dataI180, **dataE180;
    int Nx, Ny, z1, z2, N1, N2, Nc, y;
    int use_tiff;
    int search_range;
    int save_image;
    int effective_argc;
    double **array = NULL;
    double D, A, B, RMS1, RMS2, dRMS,
           sCC = 0.0, sCCy = 0.0, sCCR0 = 0.0,
           sCCy2 = 0.0, sCCR02 = 0.0, sCCyR0 = 0.0;
    double *R0_arr, *RMS0_arr, *CC_arr;
    int *dx0_arr;
    char desc[DESC_LEN],
         *cmae = "cell memory allocation error.";
    int num_threads = 16;
    double search_percent = DEFAULT_SEARCH_PERCENT;

    if (argc < 7 || argc > 11)
        Error("usage : hp2DO D000 I000 E000 D180 I180 E180 {y1 y2} {image} {search%}");

    /* Detect format from first argument */
    use_tiff = IsTiffFile(PathD000);

    /* Suppress libtiff warnings for cleaner output */
    TIFFSetWarningHandler(NULL);

    /* Parse optional search% (last argument, numeric) */
    effective_argc = argc;
    if (argc >= 8) {
        double tmp = atof(argv[argc - 1]);
        if (tmp > 0.0 && tmp <= 50.0) {
            if (effective_argc == 8 || effective_argc == 10) {
                search_percent = tmp;
                effective_argc--;
            }
        }
    }

    /* ---------------------------------------------------------------
     * Load all input files into memory
     * (only this section differs between img and tif)
     * --------------------------------------------------------------- */
    if (use_tiff) {
        fprintf(stderr, "TIFF mode: Loading files...\n");
        dataD000 = ReadTiffFile(PathD000, &Nx, &Ny, 1);
        dataI000 = ReadTiffFile(PathI000, &Nx, &Ny, 0);
        dataE000 = ReadTiffFile(PathE000, &Nx, &Ny, 0);
        dataD180 = ReadTiffFile(PathD180, &Nx, &Ny, 0);
        dataI180 = ReadTiffFile(PathI180, &Nx, &Ny, 0);
        dataE180 = ReadTiffFile(PathE180, &Nx, &Ny, 0);
    } else {
        fprintf(stderr, "img mode: Loading files...\n");
        dataD000 = ReadImgFile(PathD000, &Nx, &Ny, 1);
        dataI000 = ReadImgFile(PathI000, &Nx, &Ny, 0);
        dataE000 = ReadImgFile(PathE000, &Nx, &Ny, 0);
        dataD180 = ReadImgFile(PathD180, &Nx, &Ny, 0);
        dataI180 = ReadImgFile(PathI180, &Nx, &Ny, 0);
        dataE180 = ReadImgFile(PathE180, &Nx, &Ny, 0);
    }
    fprintf(stderr, "Loading complete.\n");

    /* Parse slice range */
    if (effective_argc == 7 || effective_argc == 8) {
        z1 = 0;
        z2 = Ny - 1;
    } else {
        if ((z1 = atoi(argv[7])) < 0 ||
            (z2 = atoi(argv[8])) >= Ny || z1 >= z2)
            Error("bad slice range.");
    }

    printf("%d\t%d\n", Nx, Ny);

    N1 = Nx - 1;
    N2 = Nx / 2;

    /* Calculate search range */
    search_range = (int)(Nx * search_percent / 100.0);
    if (search_range < 1) search_range = 1;
    if (search_range > N2) search_range = N2;
    fprintf(stderr, "Search range: +/-%d pixels (%.1f%% of width %d)\n",
            search_range, search_percent, Nx);

    /* Allocate result arrays */
    R0_arr   = (double *)MA(Ny, R0_arr);
    RMS0_arr = (double *)MA(Ny, RMS0_arr);
    CC_arr   = (double *)MA(Ny, CC_arr);
    dx0_arr  = (int    *)MA(Ny, dx0_arr);
    if (R0_arr == NULL || RMS0_arr == NULL || CC_arr == NULL || dx0_arr == NULL)
        Error("result array memory allocation error.");

    save_image = (effective_argc == 8 || effective_argc == 10);
    Nc = search_range * 2 + 1;
    if (save_image) {
        if ((array = (double **)MA(Ny, array)) == NULL) Error(cmae);
        for (y = 0; y < Ny; y++)
            if ((array[y] = (double *)MA(Nc, *array)) == NULL) Error(cmae);
    }

    /* Set number of threads */
    omp_set_num_threads(num_threads);

    /* ---------------------------------------------------------------
     * Main computation (parallel, identical for both formats)
     * --------------------------------------------------------------- */
    #pragma omp parallel
    {
        double *p000, *p180, *fD, *fI, *fE;

        p000 = (double *)MA(Nx, p000);
        p180 = (double *)MA(Nx, p180);
        fD   = (double *)MA(Nx, fD);
        fI   = (double *)MA(Nx, fI);
        fE   = (double *)MA(Nx, fE);
        if (p000 == NULL || p180 == NULL ||
            fD == NULL || fI == NULL || fE == NULL)
            Error("thread-local memory allocation error.");

        #pragma omp for schedule(dynamic)
        for (y = 0; y < Ny; y++) {
            int dx0_l, x1_l, x2_l, dx, x;
            double RMS0, RMS, R0, CC;

            CalcProjection(dataD000, dataI000, dataE000, Nx, y, p000,
                           fD, fI, fE);
            CalcProjection(dataD180, dataI180, dataE180, Nx, y, p180,
                           fD, fI, fE);

            RMS0 = CalculateRMS(N1, dx0_l = 0, x1_l = 0, x2_l = Nx,
                                p000, p180);
            if (save_image) array[y][search_range] = RMS0;

            for (dx = 1; dx <= search_range; dx++) {
                if ((RMS = CalculateRMS(N1, -dx, dx, Nx,
                                        p000, p180)) < RMS0) {
                    RMS0 = RMS; dx0_l = -dx; x1_l = dx; x2_l = Nx;
                }
                if (save_image) array[y][search_range + dx] = RMS;

                if ((RMS = CalculateRMS(N1, dx, 0, Nx - dx,
                                        p000, p180)) < RMS0) {
                    RMS0 = RMS; dx0_l = dx; x1_l = 0; x2_l = Nx - dx;
                }
                if (save_image) array[y][search_range - dx] = RMS;
            }
            R0 = (-(double)(N1 + dx0_l) / 2.0);

            CC = 0.0;
            for (x = x1_l; x < x2_l; x++)
                CC += (p000[x + dx0_l] * p180[N1 - x]);

            R0_arr[y]   = R0;
            RMS0_arr[y] = RMS0;
            CC_arr[y]   = CC;
            dx0_arr[y]  = dx0_l;
        }

        free(p000);
        free(p180);
        free(fD);
        free(fI);
        free(fE);
    }

    /* ---------------------------------------------------------------
     * Output results (sequential)
     * --------------------------------------------------------------- */
    for (y = 0; y < Ny; y++) {
        printf("%d\t%lg\t%lf\t%lf\n", y, R0_arr[y], RMS0_arr[y], CC_arr[y]);

        if (y >= z1 && y <= z2) {
            double CC = CC_arr[y];
            double R0 = R0_arr[y];
            sCC    += CC;
            sCCy   += (CC * (double)y);
            sCCR0  += (CC * R0);
            sCCy2  += (CC * (double)y * (double)y);
            sCCR02 += (CC * R0 * R0);
            sCCyR0 += (CC * (double)y * R0);
        }
    }

    D = sCCy2 * sCC - sCCy * sCCy;
    A = (sCC * sCCyR0 - sCCy * sCCR0) / D;
    B = (sCCy2 * sCCR0 - sCCy * sCCyR0) / D;
    printf("%lf\t%lf\t%lf", A, B,
           Sqrt_val((sCCR02 - A * sCCyR0 - B * sCCR0) / sCC));
    if (effective_argc == 9 || effective_argc == 10)
        printf("\t%d\t%d", z1, z2);
    putchar('\n');

    /* ---------------------------------------------------------------
     * Save output image (optional)
     * --------------------------------------------------------------- */
    if (save_image) {
        const char *out_path = argv[effective_argc - 1];
        unsigned short **cell_t;

        RMS1 = RMS2 = RMS0_arr[0];
        for (y = 0; y < Ny; y++) {
            int dx;
            for (dx = 0; dx < Nc; dx++) {
                if (array[y][dx] < RMS1) RMS1 = array[y][dx];
                if (array[y][dx] > RMS2) RMS2 = array[y][dx];
            }
        }

        dRMS = (RMS2 - RMS1) / ((double)(1 << BPS) - 1.0);

        cell_t = (unsigned short **)MA(Ny, cell_t);
        if (cell_t == NULL) Error(cmae);
        for (y = 0; y < Ny; y++) {
            cell_t[y] = (unsigned short *)MA(Nc, *cell_t);
            if (cell_t[y] == NULL) Error(cmae);
        }

        for (y = 0; y < Ny; y++) {
            int dx;
            for (dx = 0; dx < Nc; dx++)
                cell_t[y][dx] = (unsigned short)floor(
                    (array[y][dx] - RMS1) / dRMS + 0.5);
        }

        sprintf(desc, "%d\t%lf\t%lf", Nx, RMS1, RMS2);
        StoreTiffFile(out_path, Nc, Ny, cell_t, desc);

        for (y = 0; y < Ny; y++) free(cell_t[y]);
        free(cell_t);

        for (y = 0; y < Ny; y++) free(array[y]);
        free(array);
    }

    /* Free result arrays */
    free(R0_arr);
    free(RMS0_arr);
    free(CC_arr);
    free(dx0_arr);

    /* Free image data */
    FreeImageData(dataD000, Ny);
    FreeImageData(dataI000, Ny);
    FreeImageData(dataE000, Ny);
    FreeImageData(dataD180, Ny);
    FreeImageData(dataI180, Ny);
    FreeImageData(dataE180, Ny);

    return 0;
}
