/*
 * tf2DO.c - 16bit Monochrome TIFF version of hp2DO
 * 
 * Rotation axis determination from 0/180 degree projections
 * Input: 16bit monochrome TIFF files
 * Output: Text data and optional TIFF image
 *
 * Compile (Windows, Visual Studio):
 *   cl /O2 /openmp tf2DO.c libtiff.lib /Fetf2DO.exe
 *
 * Compile (Windows, MinGW-w64):
 *   gcc -O2 -fopenmp -o tf2DO.exe tf2DO.c -ltiff
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

/* Cell type for output image (16bit) */
typedef unsigned short Cell;

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

/*
 * Read entire 16bit monochrome TIFF file into memory
 * Returns allocated 2D array [Ny][Nx] of unsigned short
 */
static unsigned short **ReadTiffFile(const char *path, int *Nx, int *Ny, int firstFile)
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

    /* Allocate memory for entire image */
    data = (unsigned short **)MA(*Ny, data);
    if (data == NULL) Abort(path, "memory allocation error.");

    for (y = 0; y < *Ny; y++) {
        data[y] = (unsigned short *)MA(*Nx, *data);
        if (data[y] == NULL) Abort(path, "memory allocation error.");
    }

    /* Read all scanlines */
    for (y = 0; y < (int)height; y++) {
        if (TIFFReadScanline(tif, data[y], y, 0) < 0) {
            Abort(path, "error reading scanline.");
        }
    }

    TIFFClose(tif);
    return data;
}

/*
 * Free 2D image array
 */
static void FreeTiffData(unsigned short **data, int Ny)
{
    int y;
    for (y = 0; y < Ny; y++) free(data[y]);
    free(data);
}

static double Log(double I0, double I)
{
    return (I0 <= 0.0 || I <= 0.0) ? 0.0 : log(I0 / I);
}

/*
 * Calculate projection data from pre-loaded D(dark), I(I0), E(exposure) arrays
 */
static void CalcProjection(unsigned short **dataD, unsigned short **dataI,
                           unsigned short **dataE, int Nx, int y, double *p)
{
    int x;
    double I0, I;

    for (x = 0; x < Nx; x++) {
        I0 = (double)dataI[y][x] - (double)dataD[y][x];
        I  = (double)dataE[y][x] - (double)dataD[y][x];
        p[x] = Log(I0, I);
    }
}

static double Sqrt(double d)
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
    return Sqrt(sum / (double)(x2 - x1));
}

/*
 * Save result as 16bit monochrome TIFF
 */
static void StoreTiffFile(const char *path, int width, int height,
                          Cell **cell, const char *desc)
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
    if (desc != NULL && strlen(desc) > 0) {
        TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, desc);
    }

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

#define PathD000 argv[1]
#define PathI000 argv[2]
#define PathE000 argv[3]
#define PathD180 argv[4]
#define PathI180 argv[5]
#define PathE180 argv[6]

#define BPS 16
#define DESC_LEN 256

/* Search range as percentage of image width (default 10%) */
#define DEFAULT_SEARCH_PERCENT 10.0

int main(int argc, char **argv)
{
    unsigned short **dataD000, **dataI000, **dataE000,
                   **dataD180, **dataI180, **dataE180;
    int Nx, Ny, z1, z2, N1, N2, Nc, y;
    int search_range;  /* Actual search range in pixels */
    double **array, D, A, B, RMS1, RMS2, dRMS,
           sCC = 0.0, sCCy = 0.0, sCCR0 = 0.0, sCCy2 = 0.0, sCCR02 = 0.0, sCCyR0 = 0.0;
    double *R0_arr, *RMS0_arr, *CC_arr;
    int *dx0_arr;
    char desc[DESC_LEN],
         *cmae = "cell memory allocation error.";
    Cell **cell;
    int num_threads = 16;
    double search_percent = DEFAULT_SEARCH_PERCENT;

    if (argc < 7 || argc > 11)
        Error("usage : tf2DO D000 I000 E000 D180 I180 E180 {y1 y2} {TIFF} {search%}");

    /* Suppress libtiff warnings/errors to stderr for cleaner output */
    TIFFSetWarningHandler(NULL);

    /* Check for search percentage as last argument (if numeric) */
    if (argc >= 8) {
        double tmp = atof(argv[argc - 1]);
        if (tmp > 0.0 && tmp <= 50.0) {
            /* Last argument looks like a percentage */
            if (argc == 8 || argc == 10) {
                search_percent = tmp;
                argc--;  /* Adjust argc to exclude search% from other parsing */
            }
        }
    }

    /* Read all TIFF files into memory */
    fprintf(stderr, "Loading TIFF files...\n");
    dataD000 = ReadTiffFile(PathD000, &Nx, &Ny, 1);
    dataI000 = ReadTiffFile(PathI000, &Nx, &Ny, 0);
    dataE000 = ReadTiffFile(PathE000, &Nx, &Ny, 0);
    dataD180 = ReadTiffFile(PathD180, &Nx, &Ny, 0);
    dataI180 = ReadTiffFile(PathI180, &Nx, &Ny, 0);
    dataE180 = ReadTiffFile(PathE180, &Nx, &Ny, 0);
    fprintf(stderr, "Loading complete.\n");

    if (argc == 7 || argc == 8) {
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

    /* Calculate search range from percentage */
    search_range = (int)(Nx * search_percent / 100.0);
    if (search_range < 1) search_range = 1;
    if (search_range > N2) search_range = N2;
    fprintf(stderr, "Search range: +/-%d pixels (%.1f%% of width %d)\n", 
            search_range, search_percent, Nx);

    /* Allocate result arrays for parallel processing */
    R0_arr = (double *)MA(Ny, R0_arr);
    RMS0_arr = (double *)MA(Ny, RMS0_arr);
    CC_arr = (double *)MA(Ny, CC_arr);
    dx0_arr = (int *)MA(Ny, dx0_arr);
    if (R0_arr == NULL || RMS0_arr == NULL || CC_arr == NULL || dx0_arr == NULL)
        Error("result array memory allocation error.");

    if (argc == 8 || argc == 10) {
        if ((array = (double **)MA(Ny, array)) == NULL) Error(cmae);

        Nc = search_range * 2 + 1;  /* Use search_range instead of N2 */
        for (y = 0; y < Ny; y++)
            if ((array[y] = (double *)MA(Nc, *array)) == NULL) Error(cmae);
    }

    /* Set number of threads */
    omp_set_num_threads(num_threads);

    /* Parallel processing loop */
    #pragma omp parallel
    {
        double *p000, *p180;
        int tid = omp_get_thread_num();

        /* Thread-local buffers */
        p000 = (double *)MA(Nx, p000);
        p180 = (double *)MA(Nx, p180);
        if (p000 == NULL || p180 == NULL) {
            fprintf(stderr, "Thread %d: line memory allocation error.\n", tid);
        }

        #pragma omp for schedule(dynamic)
        for (y = 0; y < Ny; y++) {
            int dx0, x1, x2, dx, x;
            double RMS0, RMS, R0, CC;

            CalcProjection(dataD000, dataI000, dataE000, Nx, y, p000);
            CalcProjection(dataD180, dataI180, dataE180, Nx, y, p180);

            RMS0 = CalculateRMS(N1, dx0 = 0, x1 = 0, x2 = Nx, p000, p180);
            if (argc == 8 || argc == 10) array[y][search_range] = RMS0;

            /* Search only within +/- search_range from center */
            for (dx = 1; dx <= search_range; dx++) {
                if ((RMS = CalculateRMS(N1, -dx, dx, Nx, p000, p180)) < RMS0) {
                    RMS0 = RMS;
                    dx0 = (-dx);
                    x1 = dx;
                    x2 = Nx;
                }
                if (argc == 8 || argc == 10) array[y][search_range + dx] = RMS;

                if ((RMS = CalculateRMS(N1, dx, 0, Nx - dx, p000, p180)) < RMS0) {
                    RMS0 = RMS;
                    dx0 = dx;
                    x1 = 0;
                    x2 = Nx - dx;
                }
                if (argc == 8 || argc == 10) array[y][search_range - dx] = RMS;
            }
            R0 = (-(double)(N1 + dx0) / 2.0);

            CC = 0.0;
            for (x = x1; x < x2; x++) CC += (p000[x + dx0] * p180[N1 - x]);

            /* Store results */
            R0_arr[y] = R0;
            RMS0_arr[y] = RMS0;
            CC_arr[y] = CC;
            dx0_arr[y] = dx0;
        }

        free(p000);
        free(p180);
    }

    /* Output results and calculate sums (sequential) */
    for (y = 0; y < Ny; y++) {
        printf("%d\t%lg\t%lf\t%lf\n", y, R0_arr[y], RMS0_arr[y], CC_arr[y]);

        if (y >= z1 && y <= z2) {
            double CC = CC_arr[y];
            double R0 = R0_arr[y];
            sCC += CC;
            sCCy += (CC * (double)y);
            sCCR0 += (CC * R0);
            sCCy2 += (CC * (double)y * (double)y);
            sCCR02 += (CC * R0 * R0);
            sCCyR0 += (CC * (double)y * R0);
        }
    }

    D = sCCy2 * sCC - sCCy * sCCy;
    A = (sCC * sCCyR0 - sCCy * sCCR0) / D;
    B = (sCCy2 * sCCR0 - sCCy * sCCyR0) / D;
    printf("%lf\t%lf\t%lf", A, B, Sqrt((sCCR02 - A * sCCyR0 - B * sCCR0) / sCC));
    if (argc == 9 || argc == 10) printf("\t%d\t%d", z1, z2);
    putchar('\n');

    if (argc == 8 || argc == 10) {
        RMS1 = RMS2 = RMS0_arr[0];
        for (y = 0; y < Ny; y++) {
            int dx;
            for (dx = 0; dx < Nc; dx++) {
                if (array[y][dx] < RMS1) RMS1 = array[y][dx];
                if (array[y][dx] > RMS2) RMS2 = array[y][dx];
            }
        }
        dRMS = (RMS2 - RMS1) / ((double)(1 << BPS) - 1.0);

        cell = (Cell **)MA(Ny, cell);
        if (cell == NULL) Error(cmae);
        for (y = 0; y < Ny; y++) {
            cell[y] = (Cell *)MA(Nc, *cell);
            if (cell[y] == NULL) Error(cmae);
        }

        for (y = 0; y < Ny; y++) {
            int dx;
            for (dx = 0; dx < Nc; dx++) {
                cell[y][dx] = (Cell)floor((array[y][dx] - RMS1) / dRMS + 0.5);
            }
        }

        sprintf(desc, "%d\t%lf\t%lf", Nx, RMS1, RMS2);
        StoreTiffFile(argv[argc - 1], Nc, Ny, cell, desc);

        /* Free cell memory */
        for (y = 0; y < Ny; y++) free(cell[y]);
        free(cell);

        /* Free array memory */
        for (y = 0; y < Ny; y++) free(array[y]);
        free(array);
    }

    /* Free result arrays */
    free(R0_arr);
    free(RMS0_arr);
    free(CC_arr);
    free(dx0_arr);

    /* Free TIFF data */
    FreeTiffData(dataD000, Ny);
    FreeTiffData(dataI000, Ny);
    FreeTiffData(dataE000, Ny);
    FreeTiffData(dataD180, Ny);
    FreeTiffData(dataI180, Ny);
    FreeTiffData(dataE180, Ny);

    return 0;
}