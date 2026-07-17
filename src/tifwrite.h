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
 * changes, so readers get bit-identical pixels. See LIBTIFF-MIGRATION.md.
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

#endif /* CT_TIFWRITE_H */
