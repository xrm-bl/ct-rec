#ifndef CT_TIFWRITE_H
#define CT_TIFWRITE_H
/*
 * ct_write_raw_strips - write an uncompressed, contiguous row-major image as
 * multi-row raw strips (~1 MB each) for fast saves.
 *
 * Replaces the historical "1 row = 1 strip" loop
 *     for (i = 0; i < wY; i++)
 *         TIFFWriteRawStrip(image, i, data + i*wX, wX * sizeof(T));
 * which forced RowsPerStrip = 1 (slow: one strip + one write per row).
 *
 * The bytes written are identical (contiguous rows); only the strip grouping
 * changes, so readers get bit-identical pixels. See 20260717_LIBTIFF-MIGRATION.md.
 *
 * Requires <tiffio.h> to be included first.  Sets TIFFTAG_ROWSPERSTRIP itself,
 * so callers need not (and should not) set it separately.
 */
static void ct_write_raw_strips(TIFF *image, const void *data,
                                uint32_t wX, uint32_t wY, size_t elemsize)
{
    size_t   row_bytes = (size_t)wX * elemsize;
    uint32_t R, s, ns;

    if (row_bytes == 0) row_bytes = 1;
    /* ~1 MB per strip (benchmarked sweet spot for uncompressed writes). */
    R = (row_bytes >= 1048576) ? 1u : (uint32_t)(1048576u / row_bytes);
    if (R < 1u) R = 1u;

    TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, R);

    ns = (wY + R - 1u) / R;
    for (s = 0; s < ns; s++) {
        uint32_t r0   = s * R;
        uint32_t rows = (r0 + R <= wY) ? R : (wY - r0);
        TIFFWriteRawStrip(image, s,
                          (void *)((const char *)data + (size_t)r0 * row_bytes),
                          (tmsize_t)((size_t)rows * row_bytes));
    }
}

/*
 * Compile-time switch -DLZW: the 8/16-bit TIFF output of act_tg / ofact_srec
 * (rh / ro) is written LZW-compressed instead of as raw strips.  Pixel values,
 * ImageDescription and the other tags are unchanged; only the strip encoding
 * differs, and any libtiff-based reader (TIFFReadScanline /
 * TIFFReadEncodedStrip), ImageJ and rif.c decode it.  Without -DLZW the
 * output is byte-identical to before.
 */
#ifdef LZW
#define CT_TIF_COMPRESSION	COMPRESSION_LZW
#else
#define CT_TIF_COMPRESSION	COMPRESSION_NONE
#endif

/*
 * ct_write_strips - ct_write_raw_strips with the compression given.
 * COMPRESSION_NONE writes raw strips exactly as ct_write_raw_strips does;
 * any other scheme goes through TIFFWriteEncodedStrip with the same ~1 MB
 * strip grouping.  Sets TIFFTAG_COMPRESSION and TIFFTAG_ROWSPERSTRIP itself,
 * so callers set neither.
 */
static void ct_write_strips(TIFF *image, const void *data,
                            uint32_t wX, uint32_t wY, size_t elemsize,
                            uint16_t compression)
{
    size_t   row_bytes = (size_t)wX * elemsize;
    uint32_t R, s, ns;

    TIFFSetField(image, TIFFTAG_COMPRESSION, compression);
    if (compression == COMPRESSION_NONE) {
        ct_write_raw_strips(image, data, wX, wY, elemsize);
        return;
    }
    if (row_bytes == 0) row_bytes = 1;
    R = (row_bytes >= 1048576) ? 1u : (uint32_t)(1048576u / row_bytes);
    if (R < 1u) R = 1u;
    TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, R);

    ns = (wY + R - 1u) / R;
    for (s = 0; s < ns; s++) {
        uint32_t r0   = s * R;
        uint32_t rows = (r0 + R <= wY) ? R : (wY - r0);
        if (TIFFWriteEncodedStrip(image, s,
                (void *)((const char *)data + (size_t)r0 * row_bytes),
                (tmsize_t)((size_t)rows * row_bytes)) < 0) {
            fprintf(stderr, "ct_write_strips: TIFFWriteEncodedStrip failed (strip %u)\n", s);
            exit(1);
        }
    }
}

#endif /* CT_TIFWRITE_H */
