/*
 * sif_tiff.c - libtiff-based image writer
 *
 * Created: 2026-09-24
 * AI: Claude Fable 5 (Anthropic)
 *
 * Drop-in replacement for sif.c/sif_fast.c (StoreImageFile), backed
 * by libtiff 4.x. Link this file INSTEAD of sif_fast.c, together
 * with $(TIFFLIB).
 *
 * Writes uncompressed single-sample grayscale TIFF, 8 or 16 bits
 * per sample, using ct_write_raw_strips (~1MB strips, see
 * tifwrite.h / 20260717_LIBTIFF-MIGRATION.md). Pixel values are
 * identical to sif_fast.c output; only the strip layout and tag
 * order in the container differ.
 *
 * Scope: BPS 8 or 16 (other depths raise a clear error; sif.c could
 * pack 1..32 bit, but no ct-rec pipeline stores those).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "tiffio.h"
#include "cell.h"
#include "sif.h"
#include "tifwrite.h"

static void Error(char *path, char *msg)
{
    (void)fprintf(stderr, "%s : %s\n", path, msg);
    exit(1);
}

void StoreImageFile(char *path, int Nx, int Ny, int BPS,
                    Cell **cell, char *desc)
{
    TIFF *image;
    size_t noe = (size_t)Nx * (size_t)Ny;
    void *buf;
    int y, x;

    if (BPS != 8 && BPS != 16)
        Error(path, "bits per sample not supported (8/16 only).");

    /* stage into one contiguous buffer (cell rows may be scattered) */
    if ((buf = malloc(noe * (size_t)(BPS / 8))) == NULL)
        Error(path, "no allocatable memory for image buffer.");

    if (BPS == 16) {
        uint16_t *p = (uint16_t *)buf;

        for (y = 0; y < Ny; y++, p += Nx)
            for (x = 0; x < Nx; x++) p[x] = (uint16_t)cell[y][x];
    } else {
        uint8_t *p = (uint8_t *)buf;

        for (y = 0; y < Ny; y++, p += Nx)
            for (x = 0; x < Nx; x++) p[x] = (uint8_t)cell[y][x];
    }

    if ((image = TIFFOpen(path, "w")) == NULL)
        Error(path, "file not stored.");

    TIFFSetField(image, TIFFTAG_IMAGEWIDTH, (uint32_t)Nx);
    TIFFSetField(image, TIFFTAG_IMAGELENGTH, (uint32_t)Ny);
    TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, (uint16_t)BPS);
    TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
    TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, (uint16_t)1);
    TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    /* libtiff 4.x crashes on a NULL description - only set when given */
    if (desc != NULL)
        TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, desc);

    ct_write_raw_strips(image, buf, (uint32_t)Nx, (uint32_t)Ny,
                        (size_t)(BPS / 8));

    TIFFClose(image);
    free(buf);
}
