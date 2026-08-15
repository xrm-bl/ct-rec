# The three corrections applied on top of plain CBP

2026-08-15

Reconstruction in this suite is, at heart, parallel-beam convolution back
projection (CBP / filtered back projection): normalization by dark and I0,
logarithm, convolution with a ramp-type filter, back projection. That part is
textbook.

On top of it, three corrections are applied as standard. Their implementation,
principle and references are given below.

---

## 1. Ring (stripe) artifact removal - sorting-based filter (Vo's Algorithm 3)

**Implementation**: [src/sort_filter_omp.c](src/sort_filter_omp.c) (CPU/OpenMP),
[src/sort_filter_g.cu](src/sort_filter_g.cu) (GPU). Applied in every
reconstruction program immediately after the sinogram is built. Environment
variable `KERNEL_SIZE` (default 5; 0 disables it; around 11-17 is recommended
when the sample has no prominent structure).

**Principle**: each column of the sinogram (one detector channel) is sorted by
intensity along the projection-angle direction; a median filter is then applied
to the sorted image along the horizontal (channel) direction; finally the values
are put back into their original places using the recorded index matrix.
Sorting rearranges sample-borne structure into a monotonic ramp along the angle
axis, whereas a gain error fixed to one channel survives as an offset of the
whole column. The median along the channel direction therefore removes only the
latter and barely touches the sample signal. It is particularly strong on
partial stripes and does not generate extra stripe artifacts of its own.

**References**: Vo, Atwood, Drakopoulos, "Superior techniques for eliminating ring
artifacts in X-ray micro-tomography", *Opt. Express* 26(22), 28396-28412 (2018),
doi:10.1364/OE.26.028396 - this implementation corresponds to Algorithm 3 (the
sorting-based technique) of that paper. For comparison with earlier methods,
Münch et al., *Opt. Express* 17(10), 8567-8591 (2009) (wavelet-FFT method).

---

## 2. Truncation (cupping) correction - extrapolation pad on the projection data

**Implementation**: shared by the whole CBP layer - `DetectPad()` / `FillPad()`
in [src/cbp_thread.c](src/cbp_thread.c) (and `cbp_thread_int.c` /
`cbp_thread_nai.c` / `cbp_thread_avx.c`), and [src/cbp.cu](src/cbp.cu) on the
GPU. Environment variable `PAD_THRESH` (default 0.3 = automatic detection ON;
`PAD_THRESH=0` forces it OFF).

**Principle**: when the sample overfills the field of view, the projection does
not fall to zero at its ends, and the long tail of the ramp filter,
-1/(2 pi^2 r^2), picks up that step, producing cupping (the periphery lifted,
the centre depressed). This implementation judges the data truncated only when
the mean amplitude of the outermost sinogram columns exceeds `PAD_THRESH` times
the overall mean, and then appends N/2 pixels to each side of the filter input
(N = the number of horizontal pixels per projection, i.e. the sinogram width)
before the convolution. The appended values are the end values p[0] and p[N-1]
multiplied by a raised-cosine window w(n) = (1 + cos(pi n / (N/2))) / 2,
n = 1...N/2, decaying smoothly to zero from the height of the edge (plain zero
padding would leave the step, which is the very cause of the cupping). After the
inverse FFT only the leading N samples are taken, so the reconstructed region
stays the original N×N, and a sample that fits inside the field of view has ends
near zero and is left untouched.

**References**: Ohnesorge, Flohr, Schwarz, Heiken, Bae, "Efficient correction for
CT image artifacts caused by objects extending outside the scan field of view",
*Med. Phys.* 27(1), 39-46 (2000).
In addition, for the variant that goes further and extends the field of view,
Hsieh et al., "A novel reconstruction algorithm to extend the CT scan
field-of-view", *Med. Phys.* 31(9), 2385-2391 (2004).
In the area of synchrotron local tomography, Kyrieleis et al., *Nucl. Instrum.
Methods A* 607, 677-684 (2009).

---

## 3. Low-transmission guard - floor applied to the logarithm

**Implementation**: [src/blacklim.h](src/blacklim.h) (`BlackLog()`, called by each
sinogram builder). Environment variables `CT_REC_BLACK_THRESH` (default 2.0
counts above dark) and `CT_REC_BLACK_FRAC` (default 0.5). Technical note:
[20260806_low_transmission_guard.md](20260806_low_transmission_guard.md).

**Principle**: p = log(I0/(I - dark)) produces +Inf or NaN once I - dark reaches
zero or below, which happens for thick or dense samples. The guard imposes a
per-pixel floor on the signal, capping the absorbance at log(I0/threshold). In
addition, when the fraction of clipped pixels in one projection exceeds
`CT_REC_BLACK_FRAC`, that projection is judged black and replaced by the mean of
the good projections (the missing-angle case).

**References**: the physical background is photon starvation. Hsieh, "Adaptive
streak artifact reduction in computed tomography resulting from excessive x-ray
photon noise", *Med. Phys.* 25(11), 2139-2147 (1998) is the classic source.
On the bias of the logarithm at low counts, Whiting et al., "Properties of
preprocessed sinogram data in x-ray computed tomography", *Med. Phys.* 33(9),
3290-3303 (2006). The projection-replacement rule (the per-line test on `FRAC`)
is specific to this software and has no corresponding publication.
