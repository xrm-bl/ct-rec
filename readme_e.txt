Description of Image Reconstruction Software and Normalization
Based on Nakano's Software

Uesugi

2026.10.16  ver. 2.8
2026.10.08  ver. 2.7
2026.09.22  ver. 2.6
2026.09.09  ver. 2.5
2026.08.06  ver. 2.4
2026.07.30  ver. 2.3
2026.07.02  ver. 2.2
2026.06.30  ver. 2.1
2026.06.30  ver. 2.0
2026.05.04  ver. 1.7

[ver 2.8 changes]
  - New his_rec_P_F, his_srec_P_F, his_DO and his_sg (2g):
    reconstruction, rotation-axis estimation and sinograms directly from
    a.his.  The
    beamline's conv.bat (how the frames become dark.img and q????.img:
    img_ave / ren / copy) is interpreted in memory and output.log is read
    as usual, so his2img and conv.bat need not be run and no .img file is
    created.  The angle column of output.log tells a 180-degree scan from
    a 360-degree offset scan (sample views reaching 350 degrees, the test
    ofct_rec has always applied), and the matching body is run:
    his_rec = ct_rec / ofct_rec (single slice, run beside a.his; without
    center the axis is estimated from that layer, for an offset scan
    with a one-layer version of ofct_DO's search),
    his_srec = hp_tg / ofct_srec (whole volume, run one level up; hp_tg
    argument form for both), his_DO = hp2DO / ofct_DO (axis estimate),
    his_sg = sinog / of_sinog (sinogram of one layer, run beside a.his).
    Z scans (several scans in one a.his): RHIS_SCAN=k selects the scan;
    rc-check-his / gen-all-his (exe/*.bat, bin/*.sh) run the axis check
    and the whole-volume reconstruction over all scans without spl and
    conv.bat; his_info reports the number of scans and the centre layer.
    his_tb_rec_P_F / his_tb_srec_P_F (or RHIS_TB=1 with any his_*
    program) use the integer mean of all sample views in place of the
    I0 images.  his_srec_g_F_n2i reconstructs K = 2 or 4 interleaved
    view subsets into rec/0 .. rec/K-1 for Noise2Inverse (GPU only).
  - ofct_srec / ofct_rec / of_sinog: the rotation axis position of an
    offset scan may be fractional, as for a 180-degree scan.  The opposing view is
    resampled by twice the fraction (linear interpolation) and the
    rotation centre moves by the fraction; an integer value gives exactly
    the old result.  ofct_rec used to accept a fractional center but
    stitched the views at the truncated integer; this is fixed.  his_rec's
    single-layer axis estimate returns the fractional value (his_sg uses
    the same estimate).
  - ofct_srec: tilted rotation axis (L1 C1 L2 C2 form, 3b).  The axis
    position is interpolated per slice, every slice is written
    at the size of the widest one (narrower slices centred, margins hold
    the row's end values), and the ImageDescription records Rc+0.5 of the
    slice (field 2, the convention of ofct_rec).  his_srec's L1 C1 L2 C2 form
    on an offset scan uses it.  With C1 == C2 the output is identical to
    before.
  - CPU builds (cbp_thread_int.c): the fixed-point detector positions
    of the back-projection were 32-bit and overflowed for reconstruction
    widths above 5792 pixels, which only the offset programs reach
    (ofct_srec / ofct_rec / his_* on a 360-degree scan with a 4096-pixel
    detector, e.g. Rc 3809 -> 7618 pixels): the CPU image was
    geometrically wrong while the GPU build was fine.  The positions are
    now 64-bit; results for smaller widths are bit-identical to before.
    The GPU builds were never affected.
    The averages are computed exactly like img_ave (integer sum,
    truncating division), so the results are bit-identical to the
    existing chains (verified on the CPU and GPU builds; only the ARTIST
    tag of his_srec differs).  Only 16-bit his files; the number of
    dark / I0 frames and the timing of the I0 groups are taken from
    conv.bat and output.log, never assumed.  The existing programs are
    unchanged (their file-scope names were made static so that two of
    them can be linked into one his_* program).

[ver 2.7 changes]
  - New act_tg_g_F (2f) and ofact_srec_g_F (3d): whole-volume
    reconstruction with direct 16-bit + 8-bit output.  They are hp_tg /
    ofct_srec up to the reconstruction, but instead of 32-bit rec*.tif
    every slice is written at the same time as rec16/rh%05d.tif and
    rec8/ro%05d.tif, normalized with two ranges given on the command line
    (rec16 LACmin16 LACmax16 rec8 LACmin8 LACmax8).  Quantisation,
    even-size trimming and the 8-field ImageDescription are those of
    tif_f2i, so the files are identical to "hp_tg / ofct_srec -> tif_f2i"
    (verified on the CPU and GPU builds).  Written for the automatic CT
    system, but usable for ordinary CT measurements as well; note that
    the two normalization ranges must always be given explicitly - there
    is no automatic range.  GPU builds only (_r / _s / _c).
  - act_rec_g_r / ofact_rec_g_r (ver 2.6, single-slice, direct 8/16-bit)
    are retired; they are no longer built or distributed.  Use act_tg /
    ofact_srec, or ct_rec / ofct_rec followed by tif_f2i.
  - hp_tg: in the pipelined build (ver 2.6 update) the rotation-centre
    field of the ImageDescription could belong to the next slice for the
    tilted-axis form (L1 C1 L2 C2); pixel values were not affected.  Fixed.
  - Windows build: sources with Japanese comments now carry a UTF-8 BOM
    and cl is called with /utf-8, so the C4819 code-page warnings are
    gone.  No effect on Linux or on the generated code.
  - ct_rec / ofct_rec (2026-10-13 update): when every projection of a
    slice is judged black (no X-rays at all, see 1e), a placeholder image
    is written under the normal name and size instead of nothing.  It
    shows the file name, "fail!" and "all projections black" drawn in
    large letters (background 0, text 1.0), and its ImageDescription ends
    with "FAILED: all projections black".  The exit status stays non-zero
    (1) so batch files still see the failure.  ofct_rec used to write an
    all-zero image with exit status 0 in this case; it now behaves like
    ct_rec.  Slices with at least one good projection are unchanged.

[ver 2.6 changes]
  - act_rec_g_r / ofact_rec_g_r: single-slice reconstruction for the
    automatic CT system with direct 8/16-bit output (retired in ver 2.7).
  - All GPU tools: clearer messages on GPU memory exhaustion.  The
    free/total GPU memory and a hint are appended, and cuFFT errors now
    show the error-code name.  The reconstruction (cbp) layer estimates
    its total device demand before allocating and warns when it exceeds
    the free memory (on Windows the driver silently falls back to
    system RAM instead of failing, so this warning is the only visible
    sign of the shortage).
  - act_spl2 / act_spl (Linux): the default upper limit of concurrent
    background filter jobs (environment variable ACT_SPL_JOBS) is raised
    from 8 to 16.  The limit itself is verified to work (16 by default,
    3 with ACT_SPL_JOBS=3); if it seems ineffective, the binary predates
    2026-07-23 or tif_mgf / gf_sd is not on the PATH.  Windows is
    unchanged (start /b, no limit).
  - New note 20260919_ct_rec_gpu_memory.md: GPU memory usage of ct_rec_g_*
    and an estimate of how many instances fit on one GPU.
  - hp_tg / ofct_srec: less per-slice overhead, same results.  Every slice
    used to malloc/free about 1.4 GB of work buffers (ring removal input and
    output, output image), which cost ~185 million page faults per run and
    most of the "sys" time, and five single-threaded copies of the sinogram
    and tomogram.  Now the ring removal runs in place on the CBP projection
    buffer (pinned on the GPU build, so the host<->device copies run at full
    speed), hp_tg writes the TIFF directly from the reconstruction buffer,
    ofct_srec keeps its float conversion buffer, and the remaining copies
    (projections in, black-projection check, tomogram out) are OpenMP
    parallel.  sort_filter_g.cu keeps its device buffers between calls
    instead of cudaMalloc/cudaFree x4 per call (new sort_filter_gpu_release();
    result_data may now equal image_data).  The order of every arithmetic
    operation is unchanged, so the output is byte-identical (verified on the
    CPU build with a synthetic data set including black projections).
    On top of that, hp_tg / ofct_srec now run a one-slice-ahead pipeline
    (cbp_pipe.h): the projection buffer is double-buffered, and while the
    GPU reconstructs slice z the host already prepares slice z+1
    (projection copy, black-projection check, ring removal on its own CUDA
    stream); the TIFF copy/write of slice z overlaps the GPU work of z+1.
    cbp.cu gained AllocCBPProjection() / SelectCBPProjection() for this; the
    existing CBP() / BeginCBP() / EndCBP() are unchanged, and the angle table
    is now uploaded before the FFTs so that BeginCBP() returns without
    waiting for the filter kernels.  Note
    for ct_rec_g_* run many-at-a-time: the ring-removal device buffers
    (3 x Nx x Nt x 4 bytes) now stay allocated until the process exits.
  - p_rec: two fixes and the same speed-ups as hp_tg.  (1) The read loop
    skipped the first p-file, so projection 0 was always empty (one view
    out of Nt missing).  Every file is now used, which changes the result
    by exactly that one view; files may be numbered from 0 or from 1.
    (2) The outlier test applied the integer abs() to float data, so it
    never rejected NaN/Inf and could leave the previous slice's values in
    place; it now uses fabsf() and zeroes rejected pixels (threshold 100
    unchanged, so normal data is unaffected).  Reading is parallel
    (HPTG_READ_THREADS) and fetches only the strips of the rows a band
    needs; ring removal runs in place, the output buffer is reused, unit
    conversion and min/max are OpenMP, the TIFF is written by a separate
    thread, and the one-slice-ahead pipeline of hp_tg is used
    (HPTG_PIPELINE=0 for the sequential order).  Verified identical to the
    old version (apart from fix 1) on the CPU and GPU builds.

[ver 2.5 changes]
  - The projection-image generator ct_prj_f gained optional Paganin
    single-distance phase retrieval (parallel beam).  It is active only
    when a parameter file pr.par is given as the third argument (see 6h);
    without it the behaviour and the results are unchanged.
  - P_SIZE in pr.par (the effective pixel size) is in MICROMETRES.  The
    old mkpms/mkpmbg chain took it in cm, so reusing an old file is a
    4-decade mistake; values below 1 nm are rejected as suspected
    cm-unit files.
  - Boundaries are mirror-padded (width = 10 x the kernel decay length
    rc = sqrt(R2 d/mu), at least 32 pixels).  Verified on synthetic data
    to agree with the legacy mkpms + mkpmbg + pm_sino chain after its
    background subtraction (difference ~4e-7, the 32-bit float
    quantisation level).  Details: 20260909_paganin_phase_retrieval.md
    (Japanese).
  - The FFT is built in (radix-2, single precision, real-packed, OpenMP
    parallel), no external library.  About 0.1 s per 2048 x 2048
    projection on 16 threads.  Thread count: environment variable
    PAGANIN_THREADS (default: the number of logical CPUs, independent
    of OMP_NUM_THREADS).
  - New tif2zar: numbered-TIFF to OME-Zarr converter (6s).  Turns a
    reconstructed volume into the multiscale, blosc-zstd-compressed
    format for Fiji + BigDataViewer/MoBIE in one streaming pass.
    Details: 20260909_tif2zar.md (Japanese).

[ver 2.4 changes]
  - Added a per-pixel guard against insufficient transmittance. When a thick
    or dense sample drives the transmitted signal to (or below) the dark
    level, log(I0/(I-dark)) returned +Inf/NaN and the GPU ring removal
    aborted with
      CUDA error sort_filter_g.cu:266: an illegal memory access was encountered
    (the CPU build silently corrupted the sort of those columns instead).
    Fixed at the source.
  - The meaning of the environment variable CT_REC_BLACK_THRESH changed from
    "black-projection test on the line average" to "per-pixel signal floor"
    (default 1 -> 2). The unit (counts above dark) is unchanged. Pixels at or
    below the floor are clamped, capping the absorbance at
    log(I0/CT_REC_BLACK_THRESH), so Inf/NaN can no longer occur. Pixels above
    the floor are untouched.
  - The black-projection (missing-angle) decision is now made by the new
    variable CT_REC_BLACK_FRAC (fraction of clipped pixels, default 0.5).
    Plate-like samples trigger at the same transmittance as before. If an
    ultra-dense sample covering nearly the whole field of view loses too many
    projections, raise it to about 0.8. See 1e for details.
  - Both variables are validated (bad values fall back to the default).
    A normal run prints nothing; the values in effect are shown only when a
    projection is judged black, and the end-of-run summary (clipped pixels,
    minimum signal, maximum absorbance) only when pixels were clipped.
  - Also fixed on the way: an off-by-one in the air-reference lookup (it read
    uninitialised memory, so reconstruction values change slightly), a
    division by zero when every projection is judged black, and the
    unguarded log in the rotation-centre search. See also the technical note
    20260806_low_transmission_guard.md.

[ver 2.3 changes]
  - Added a truncation (cupping) correction to the CBP layer shared by every
    reconstruction program. It is controlled by the environment variable
    PAD_THRESH; when unset the default is 0.3 (auto-detection ON).
    PAD_THRESH=0 forces it OFF (results identical to the previous version).
    See 1f.
  - The default CPU thread count changed from a fixed 8 to "number of logical
    cores - 1" (environment variable CBP_THREADS). The included GPU
    executables are now built with CUDA Toolkit 13.2 (CUDA 13 requires
    Turing / sm_75 or newer).
  - ofct_DO gained view decimation (stride); ofct_DO_g gained view-pair
    averaging (navg) and decimation (stride), so the axis search can be run
    in much less time.
  - hp_tg / ofct_srec now read projections in parallel (environment variable
    HPTG_READ_THREADS, default 16).
  - Added ct_rec_loop, a wrapper that reproduces hp_tg's whole-volume
    reconstruction by running ct_rec in parallel (2d). Both a Windows .bat
    and a Linux .sh are provided, and dark.tif input is supported.
  - Documented act_spl2 / act_spl, the split-and-save tools for continuously
    acquired tif data (6q).
  - Added rec2rec, re-projection / re-reconstruction of finished CT slices (6r).
  - Added tif_mgf_g, a GPU version of the median + gaussian filter (6p).
  - The bundled libtiff was updated from 3.6.0 to 4.6.0. Reading and writing
    8/16/32-bit tiff is fully compatible with the previous version.
  - The full list of environment variables (with defaults) is now collected in
    20260723_CT_env_vars.md.

[ver 2.2 changes]
  - Added a GPU version of the offset-CT rotation-axis finder, ofct_DO_g
    (ofct_DO.cu). Host-side processing is identical to ofct_DO, so the
    estimated center/Oy match the CPU version; the per-view-pair MSD is
    computed on the GPU. The CPU ofct_DO is retained.

[ver 2.1 changes]
  - Merged ct_rec and tf_rec. ct_rec now auto-detects the input format
    from the dark frame (dark.img -> .img, otherwise dark.tif -> .tif),
    using ct_rec_c.c. The TIFF-only tf_rec is therefore retired:
      tf_rec_P_F -> ct_rec_P_F
  - The offset-CT single-slice otf_rec gained the same img/tif auto-detection
    and is renamed (the TIFF-only otf_rec is retired):
      otf_rec_P_F -> ofct_rec_P_F
  - ct_prj_f (projection-image generator) also gained img/tif auto-detection
    (ct_prj_f.c + tf_prj_f.c -> ct_prj_f_c.c); the TIFF-only tf_prj_f is retired.
  - ofct_xy (offset-CT rotation-axis finder) replaced by ofct_DO (img/tif
    auto-detection; prints a ready-to-run ofct_srec command).
  - Offset-CT reconstruction auto-detects available memory and chunks
    automatically; the manual memory-check step (old 3b) was removed.

[ver 2.0 changes]
  - Merged hp_tg and tf_tg. hp_tg now auto-detects the input format from
    the dark frame in the directory (dark.img -> .img, otherwise
    dark.tif -> .tif), using rhp_c.c (a merge of rhp.c[.img] and rtf.c[.tif]).
  - The TIFF-only programs are therefore retired and unified into:
      tf_tg_P_F     -> hp_tg_P_F
      oftf_srec_P_F -> ofct_srec_P_F
      oftf_xy       -> ofct_xy

0. Please contact the author for bug reports or feature requests.

1. General Concepts
   a. Input
      Basically img format. For single-slice (ct_rec), continuous reconstruction (hp_tg) and offset CT (ofct_srec / ofct_DO), the input format is auto-detected: dark.img -> .img, otherwise dark.tif -> .tif.

   b. Output
      CT images are output as 32-bit TIFF files: rec?????.tif (5-digit numbering)
      The following information is embedded in the TIFF tags of each image,
      in this order:
      pixel size, rotation axis position, number of projections,
      rotation angle offset, minimum and maximum values of the image.
      During normalization, the min/max values used for normalization are
      appended to the existing tags.
      For continuous reconstruction and normalization, a log is saved to
      cmd-hst.log upon completion. The rotation-axis finder (ofct_DO), the
      conversion / averaging utilities and the tif_* filters also append to
      cmd-hst.log. The tif_* filters record the command and its parameters
      as a single tab-separated line.

   c. Program Suffixes
      The reconstruction software has suffixes such as _t_c.
      These specify the processor and reconstruction filter to use.
      _P: Processor
          _t: Use CPU multi-threading. Controlled by the environment
              variable CBP_THREADS. The default is the number of logical
              cores of the machine minus 1 (minimum 1); it used to be a
              fixed 8.
          _g: Use GPGPU. The included executables are compiled with
              CUDA Toolkit 13.2 (CUDA 13 requires Turing / sm_75 or newer).
      _F: Filter
          _c: Chesler filter
          _s: Shepp-Logan filter
          _r: Ramachandran (HAN) filter
      There are two exceptions to this naming: the sinogram reconstruction
      sf_rec is CPU only (sf_rec_t_F), and the re-projection tool rec2rec has
      no _t on its CPU build (rec2rec_F / rec2rec_g_F).

   d. Ring Artifact Removal
      Since version 1.4, a ring removal function based on Vo et al. (2018)
      Algorithm 3 has been implemented. It is executed immediately before
      the CBP computation. This function can be turned ON/OFF by setting
      an environment variable.
      Setting KERNEL_SIZE to 1 turns it OFF. Other positive odd numbers
      change the strength of the effect. The default value is 5 (also used
      when the environment variable is not defined).
      The ring removal processing uses OpenMP-based CPU parallelization,
      with a default of OMP_NUM_THREADS=40. This is independent of
      CBP_THREADS (for back-projection computation) described in section 1c.

   e. Insufficient Transmittance and Missing Angle Handling
      When the sample is thick or dense, the transmitted signal can reach or
      fall below the dark level. I-dark then becomes zero or negative and
      log(I0/(I-dark)) returns +Inf or NaN, which used to abort the GPU ring
      removal with

        CUDA error sort_filter_g.cu:266: an illegal memory access ...
        ring removal image processing failed

      and, on the CPU build, silently corrupted the sort for those columns.

      CT_REC_BLACK_THRESH is a PER-PIXEL floor on the signal, in counts above
      dark (default 2). A pixel at or below it is treated as opaque and the
      floor is used in its place, so the absorbance is capped at
      log(I0/CT_REC_BLACK_THRESH) and can never become infinite. Only the
      affected pixels are touched; everything above the floor is unchanged.
      Values of 0 or below are rejected (the guard cannot be switched off).
      A normal run prints nothing. The values in effect are echoed to stderr
      only when a projection is judged black (both thresholds exceeded), and
      the end-of-run summary (clipped pixels, minimum signal, maximum
      absorbance) appears only when at least one pixel hit the floor.

      For plate-like samples, transmittance can drop drastically at certain
      angles (the missing-angle case). A projection whose clipped fraction
      exceeds CT_REC_BLACK_FRAC (default 0.5) is judged black and replaced by
      the mean of the good projections. Change CT_REC_BLACK_THRESH to 1, 10,
      100, 1000, etc. as before to adjust how readily that happens: for a
      uniformly attenuated line the trigger point is the same as it was.

      Note that the threshold now also caps the absorbance of projections that
      are kept, so a large value visibly flattens a dense sample. Prefer the
      default unless you are handling a plate.

   f. Truncation (Cupping) Correction
      When the sample overfills the field of view (truncated data), the step
      left at both ends of each projection, convolved with the long tail of
      the ramp filter, produces a bright rim (cupping) around the outside of
      the reconstruction. To suppress it, the CBP layer shared by every
      reconstruction program extends each projection on both sides by N/2
      samples before filtering, holding the edge value under a cosine decay.
      The cost of back-projection is unchanged (only the FFT length doubles).
      Suppressing the truncation artifact of an object that extends beyond
      the field of view by extrapolating the projection data follows
      Ohnesorge et al. (2000), Med. Phys. 27(1), 39-46. The implementation
      here is a simplified variant (edge-value hold with a cosine decay), not
      a reproduction of the method described in that paper.
      Whether the pad is applied is decided automatically: it is enabled only
      when the mean amplitude of the sinogram's outermost columns exceeds a
      given ratio of the overall mean. That ratio is the environment variable
      PAD_THRESH, and when it is unset the default is 0.3 (auto-detection ON).
      When the pad becomes active the program prints
        PAD_THRESH: truncation pad enabled (W=...)
      A sample that fits inside the field of view has edge values of nearly
      zero, so the pad is automatically not applied. Setting PAD_THRESH=0
      (or negative) forces it OFF, giving results identical to the previous
      version.
      This applies to ct_rec / hp_tg / his_rec / his_srec / act_tg / p_rec / ofct_rec /
      ofct_srec / ofact_srec / sf_rec / rec2rec alike. Note that the input of rec2rec is already attenuated at
      the edges by the reconstruction circle, so a smaller value of about
      0.1-0.2 is appropriate there.

   g. Other Environment Variables
      Besides the above there are the number of parallel projection-read
      threads (HPTG_READ_THREADS, default 16; effective in hp_tg / ofct_srec),
      the memory-based chunking controls (HPTG_MEM_FRACTION /
      HPTG_MEM_LIMIT_MB / HPTG_CHUNK_ROWS), the CUDA device number to use
      (CUDA_GPU, default 0) and the input file name overrides (RHP_O / RHP_D /
      RHP_Q; RHIS_H / RHIS_C / HIS_SREC_OY for his_rec / his_srec / his_DO).
      See 20260723_CT_env_vars.md for the full list with defaults.
      When several processes are run on the same node, lower
      HPTG_MEM_FRACTION; with the default every process asks for 90% of the
      free memory.

2. 180-degree Scan: Standard Absorption CT Reconstruction

   a. Single Slice Reconstruction
      ct_rec_P_F layer {center} {pixel size} {offsetangle}

      layer: Layer (height) to reconstruct
      center: Rotation axis position (pixels). Auto-estimated if omitted.
      pixel size: Pixel size (um). Defaults to 1.0 if omitted.
      offset angle: Rotation axis origin offset. Defaults to 0.0 if omitted.

      *) Run in the directory containing q????.img or q????.tif files (auto-detected from dark.img / dark.tif).

   b. Continuous Reconstruction
      hp_tg_P_F HiPic Dr RC RA0 rec
      (When the rotation axis is not tilted. All layers.)

      HiPic: Directory containing q????.img or q????.tif files (no trailing /)
      Dr: Pixel size (um)
      RC: Rotation axis position
      RA0: Rotation axis origin offset
      rec: Output directory for reconstructed images (must be created
           before execution)

      hp_tg_P_F HiPic Dr L1 C1 L2 C2 RA0 rec
      (When the rotation axis is tilted, or for partial region computation.)

      HiPic: Directory containing q????.img or q????.tif files (no trailing /)
      Dr: Pixel size (um)
      L1: Start layer
      C1: Rotation axis position at L1
      L2: End layer
      C2: Rotation axis position at L2
      RA0: Rotation axis origin offset
      rec: Output directory for reconstructed images (no trailing /.
           Must be created before execution.)

      *) Run one directory above the directory containing q????.img files.

   c. Continuous Reconstruction from p-images
      p_rec_P_F p rec Dr RC RA0
      (When the rotation axis is not tilted. All layers.)
      (p-files may be numbered from 0 or from 1; every file is used.)

      p: Directory containing p?????.tif files (no trailing /)
      rec: Output directory for reconstructed images (must be created
           before execution)
      Dr: Pixel size (um)
      RC: Rotation axis position
      RA0: Rotation axis origin offset

      p_rec_P_F p rec Dr L1 C1 L2 C2 RA0
      (When the rotation axis is tilted.)

      p: Directory containing p?????.tif files (no trailing /)
      rec: Output directory for reconstructed images (must be created
           before execution)
      Dr: Pixel size (um)
      L1: Start layer
      C1: Rotation axis position at L1
      L2: End layer
      C2: Rotation axis position at L2
      RA0: Rotation axis origin offset

   d. Continuous Reconstruction by Running ct_rec in Parallel
      ct_rec_loop.bat HiPic Dr RC RA0 rec {Njobs}            (Windows)
      ct_rec_loop.bat HiPic Dr L1 C1 L2 C2 RA0 rec {Njobs}
      ct_rec_loop.sh  HiPic Dr RC RA0 rec {Njobs}            (Linux; in bin/)
      ct_rec_loop.sh  HiPic Dr L1 C1 L2 C2 RA0 rec {Njobs}

      The arguments follow hp_tg; only the number of parallel jobs Njobs is
      added at the end (default 1).
      This wrapper produces the same whole-volume reconstruction as hp_tg by
      running the single-slice program ct_rec once per layer, several layers
      at a time. Use it when hp_tg cannot be used, or to load the GPU more
      heavily.
      The input format is auto-detected from dark.img / dark.tif in HiPic.
      In the fixed-center (RC) form, the layer range is taken automatically
      from the height stored in the dark file.
      Results are collected in rec/ and the per-layer logs in rec/log/.
      The executable used is ct_rec_g_c by default; the .sh version can be
      pointed elsewhere with the environment variable CT_EXE.
      Tune Njobs to the available GPU memory.

      *) As with hp_tg, run one directory above the HiPic directory.

   f. Continuous Reconstruction for the Automatic CT System
      (direct 16/8-bit output)
      act_tg_g_r HiPic Dr RC RA0 rec16 LACmin16 LACmax16 rec8 LACmin8 LACmax8
      act_tg_g_r HiPic Dr L1 C1 L2 C2 RA0 rec16 LACmin16 LACmax16 rec8 LACmin8 LACmax8

      HiPic, Dr, RC / L1 C1 L2 C2, RA0: as in hp_tg (2b)
      rec16: Output directory for the 16-bit images rh%05d.tif
      LACmin16 / LACmax16: Normalization range of the 16-bit images
      rec8: Output directory for the 8-bit images ro%05d.tif
      LACmin8 / LACmax8: Normalization range of the 8-bit images

      Identical to hp_tg up to the reconstruction; instead of 32-bit
      rec*.tif, every slice is written at the same time as a 16-bit and an
      8-bit image with the two ranges given on the command line (both
      directories must exist before execution; each LACmin must be smaller
      than its LACmax).  Quantisation, even-size trimming and the 8-field
      ImageDescription are the same as tif_f2i (section 4), so the files
      are identical to what "hp_tg -> tif_f2i" produces.  All environment
      variables of hp_tg apply.  GPU build only; the filter suffix is as
      in hp_tg (act_tg_g_r / act_tg_g_s / act_tg_g_c).

      *) As with hp_tg, run one directory above the HiPic directory.

   g. Reconstruction Directly from a.his (his_rec / his_srec / his_DO / his_sg)
      his_rec_P_F  layer {center} {pixel size} {offsetangle}
      his_srec_P_F HiPic Dr RC RA0 rec
      his_srec_P_F HiPic Dr L1 C1 L2 C2 RA0 rec
      his_DO       HiPic {stride}                 (CPU)
      his_DO_g     HiPic {navg {stride}}          (GPU)
      his_sg       layer {center} {skip}          (CPU)
      his_tb_rec_P_F / his_tb_srec_P_F            (arguments as above)
      his_srec_g_F_n2i HiPic Dr RC RA0 rec K      (GPU; also the L1 C1 L2 C2 form)

      The input is the data as it comes from the beamline: a.his (every
      frame of the scan), conv.bat and output.log in one directory.
      his_rec and his_sg run in that directory (like ct_rec); his_srec and
      his_DO run one level above it with the directory name as HiPic (like hp_tg
      / ofct_DO).  his2img and conv.bat are not run and no .img file is
      written: conv.bat is interpreted in memory.
        his2img <his>                 name of the his file (default a.his)
        img_ave <frames>... dark.img  dark = average of the listed frames
        img_ave <frames>... q<N>.img  q<N> = average of the listed frames
        ren <frame> q<N>.img          q<N> = that frame
        copy q<M>.img q<N>.img        q<N> = same data as q<M>
        hp2do D000 I000 E000 D180 I180 E180 ...   remembered for his_DO
        rem / echo / pause            ignored
      <frame> is a<digits>.img as his2img names it (a0001.img = first frame
      of the file).  Any other command stops the program, as do frame
      numbers outside the file, a q defined twice, or a missing dark.img.
      How many frames are averaged and where the I0 groups sit in time
      are therefore read from conv.bat and output.log, never assumed.

      180-degree or 360-degree offset scan is decided from the angle
      column of output.log: when the sample views reach 350 degrees the
      data are an offset scan (the test ofct_rec has always applied).
        program    180-degree scan   360-degree offset scan
        his_rec    ct_rec (2a)       ofct_rec (3c)
        his_srec   hp_tg (2b)        ofct_srec (3b) with Rc = RC, Oy = 0,
                                     rangeList "-", or its L1 C1 L2 C2 form
        his_DO     hp2DO on the views of the hp2do line of conv.bat
                                     ofct_DO / ofct_DO_g (3a)
        his_sg     sinog (6b)        of_sinog (6b) at center
      When center is omitted, his_rec estimates the rotation axis from the
      layer it reconstructs: on a 180-degree scan with ct_rec's own search,
      on an offset scan with a one-dimensional version of ofct_DO's search
      (the -log transmittances of that layer, smoothed with the same
      Gaussian OFCT_DO_SMOOTH, mirrored opposing views, mean squared
      difference accumulated over all view pairs for shifts in ofct_DO's
      range, parabolic refinement; the smoothing margin is kept out of the
      sum, which otherwise biases the minimum).  The fractional result is
      used as it is.  Only that layer is read, so it takes
      seconds; the result and the ofct_rec command used are printed on
      stderr, with a warning when the minimum sits at the edge of the
      search range (axis outside 0.75*Nx .. Nx-1, which ofct_DO cannot
      handle either).
      The command line of his_srec is always the hp_tg form.  For an
      offset scan RC may be fractional (3b) and Oy is 0 unless the
      environment variable HIS_SREC_OY is set; the L1 C1 L2 C2 form is
      handed to the tilted-axis form of ofct_srec (3b): the axis position
      is interpolated per slice and every slice is written at the size of
      the widest one.  his_DO on a 180-degree
      scan uses hp2DO's estimate on the six views named on the hp2do line
      (dark, I0 and view at 0 deg; dark, I0 and view at 180 deg; without
      such a line q1 q2 and the two views after the last logged one, as in
      ct_rec's centre search), ignores the stride / navg arguments, and
      prints hp2DO's output on stdout; on an offset scan it is ofct_DO.
      In both cases a ready-to-run his_srec command is printed on stderr
      (with 1.0 standing for the pixel size Dr).
      his_sg writes the sinogram of one layer exactly as sinog does
      (s%05d.tif, 32-bit float; a line in ../cmd-hst.log).  On
      an offset scan the two halves are stitched by of_sinog at center,
      which may be fractional; "-" or no center uses the single-layer
      estimate of his_rec described above, and the of_sinog command used
      is printed on stderr.  On a 180-degree scan center is not needed
      and is ignored (write "-" to give skip).
      The averages use the arithmetic of img_ave (integer sum, truncating
      division), so dark / I0 data, the reconstructions and the sinograms
      are bit-identical to the ct_rec / hp_tg / ofct_rec / ofct_srec /
      sinog / of_sinog chains (his_srec only differs in the ARTIST tag of
      the TIFF), and his_DO prints exactly what hp2DO / ofct_DO print.  Only 16-bit his files
      (type 2) are accepted.  Environment variables: RHIS_H (his file
      name, overrides the his2img line), RHIS_C (conv.bat name, default
      conv.bat), RHP_O (output.log path), RHIS_SCAN (Z scans, see
      below), RHIS_TB (mean of the sample views as I0, see below),
      RHIS_VIEWS (view subset k/K of his_srec, see below), HIS_SREC_OY
      (see above) and
      everything that applies to the underlying programs.  Frames are
      read with positional reads on one shared file handle, so the
      parallel projection reading of hp_tg (HPTG_READ_THREADS) and the
      memory-adaptive chunking work unchanged.

      his_tb_rec and his_tb_srec are his_rec and his_srec with the
      incident beam taken from the data themselves: every q that
      output.log marks as an I0 view reads as the integer mean (as
      img_ave) of all sample views listed before the first blank line of
      output.log (the extra view after the blank lines is not included),
      so the result equals a run of ct_rec / hp_tg / ofct_rec / ofct_srec
      after img_ave of those frames has been copied over the I0 images.
      Dark subtraction and everything else are unchanged; the I0 time
      interpolation becomes constant.  his_tb_srec reads the whole file
      once more to form the mean.  The same is obtained from any his_*
      program (his_DO, his_sg, the Z-scan scripts with HIS_REC /
      HIS_SREC) with the environment variable RHIS_TB=1; RHIS_TB=0 turns
      it off in a his_tb_* program.  The start-up line reports
      "I0 = mean of N sample views (RHIS_TB)" and cmd-hst.log carries
      "% RHIS_TB".

      his_srec_g_F_n2i is his_srec for Noise2Inverse: the sample views
      (in angle order) are split into K = 2 or 4 interleaved subsets
      (k, k+K, k+2K, ...) and each subset is reconstructed as his_srec into
      rec/k, with RA0 shifted by the angle of the subset's first view
      (taken from output.log) so that the K volumes share one geometry;
      RC and Dr are the same for every subset, so give the rotation axis
      explicitly.  The number of sample views must be a multiple of K and
      the scan must be a 180-degree scan (an offset scan is refused).
      The subsets run one after another as child processes of the
      program; each child is an ordinary his_srec run with the
      environment variable RHIS_VIEWS=k/K, which can also be set by hand
      with his_srec to reconstruct a single subset.  cmd-hst.log carries
      "% RHIS_VIEWS=k/K".  RHIS_SCAN and RHIS_TB apply to every subset.

      Z scans.  When a.his holds M scans of L frames each (L = the highest
      frame number conv.bat names, M = frames / L; the layout spl /
      his_spl_E split into 001/raw .. 00M/raw with the same conv.bat and
      output.log), the environment variable RHIS_SCAN = 1..M selects the
      scan for every his_* program; such a file is refused without it, so
      the first scan is never used by mistake.  The scripts rc-check-his
      and gen-all-his (exe/*.bat for Windows, bin/*.sh for Linux) replace
      the spl + conv.bat + ct_rec / hp_tg chain of 6m:
        his_info {dir}
            prints "Nx Ny frames L M layer kind" of a.his + conv.bat +
            output.log: frame size, frames in the file, frames per scan L,
            number of scans M, the centre layer Ny/2 and the scan kind 180
            or 360 (readable summary on stderr).
        rc-check-his {layer}
            M, the default layer (centre row) and the kind come from
            his_info; his_rec on that layer of every scan with its axis
            estimate; the axis goes to center.log (column 2) and the
            images to rc-check/00k.tif for a visual check.  center.log
            holds the value to pass to his_srec: for a 360-degree scan 0.5
            is subtracted from the ImageDescription, because ofct_rec
            records the image centre dcnt+0.5 while ofct_srec takes the
            detector axis position Rc = dcnt (ct_rec / hp_tg agree, so a
            180-degree scan needs no correction).
        gen-all-his center.log Dr {LACmin LACmax}
            his_srec of every scan at the axis of center.log into 00k/rec,
            then tif_f2i (8 bit, default -0.5 3.0) into 00k/ro_xy (no
            ro_zx reslice); rec_stk (6f) stacks them as before.
        his_zrec scan layer {center}
            one trial slice: his_rec on layer of that scan (axis estimated
            when center is omitted), kept as rc-check/SSS_LLLLL_C.tif
            (scan %03d, layer %05d, C = the axis used as the his_srec
            value, e.g. 003_01024_1023.5.tif) for comparing layers and
            axis positions.
      The scripts run in the directory of a.his (his_srec takes "." as
      HiPic); HIS_REC / HIS_SREC override the program names (default
      his_rec_g_c / his_srec_g_c).

      *) his_rec, his_sg: run in the directory containing a.his, conv.bat
         and output.log.  his_srec, his_DO: run one directory above it.

3. 360-degree Scan (Offset CT): Standard Absorption CT Reconstruction
   a. Rotation Axis Position Estimation
      ofct_DO   raw {stride}          (CPU)
      ofct_DO_g raw {navg {stride}}   (GPU; same result as ofct_DO)

      raw: Directory containing q????.img or q????.tif files (no trailing /;
           input format auto-detected from dark.img / dark.tif)
      stride: Use only every stride-th view pair. Defaults to 1 (use all).
              This tool is limited by projection reading, so the run time
              becomes roughly 1/stride while the estimated center/Oy hardly
              changes.
      navg: GPU version only. Number of view pairs averaged before taking
            -log. Defaults to 10; the averaging raises the S/N of the search.

      Estimates the rotation-axis position from the offset-CT data and prints
      a ready-to-run ofct_srec command (suggested center and Oy).
      ofct_DO_g computes the per-view-pair MSD on the GPU (host-side
      processing is identical, so the estimate matches the CPU version).
      Both versions apply a gaussian smoothing before the MSD computation.
      Its sigma is set by the environment variable OFCT_DO_SMOOTH (default
      10.0); specifying 0 disables the smoothing.

   b. Reconstruction
      ofct_srec_P_F HiPic Rc Oy rangeList Dr RA0 rec
      ofct_srec_P_F HiPic L1 C1 L2 C2 Oy rangeList Dr RA0 rec
      ofct_srec_P_F HiPic Rc Oy
           (the third form reconstructs nothing: it prints the stitched
            width N, the padded length L, the number of view pairs M and
            the size of the projection buffer 4*N*L*M in GB, as a check
            before a large run)

      HiPic: Directory containing q????.img or q????.tif files (no trailing /)
      Rc: Rotation axis position (pixels from left edge).  May be
           fractional (ver 2.8): the stitching frame follows the integer
           part, the opposing view is resampled 2f pixels to the right
           by linear interpolation and the rotation centre moves by f,
           as for a 180-degree scan.  An integer value gives exactly the
           old result.
      L1 C1 L2 C2: Tilted rotation axis: position C1 at layer L1 and C2 at
           layer L2 (ver 2.8).  Every slice uses the position interpolated
           linearly between them (fractional, see Rc).  The stitched width
           2*floor(Rc(z)) then differs from slice to slice, so every slice is
           written at the size of the widest selected slice; a narrower
           slice is centred in that frame and the margin of its sinogram
           holds the end values of each row (not zero: the truncation
           check of the reconstruction looks at the outermost columns, and
           zero margins would switch the cupping correction off for the
           narrower slices only).  The rotation axis is at the centre of
           every output image, so the slices stack.  With C1 == C2 the
           result is identical to the first form.  Near the image edge the
           values differ slightly from a reconstruction at the slice's own
           width: about 0.1 % of the value range in the centre and up to a
           few per cent at the rim on test data.  L1 and L2 are output
           layer numbers (equal to detector rows for Oy = 0); layers of
           rangeList outside L1..L2 are extrapolated.
      Oy: Vertical shift (always specify 0)
      rangeList: Specify layers to reconstruct.
      Dr: Pixel size
      RA0: Rotation axis origin offset
      rec: Output directory for reconstructed images (must be created
           before execution)

      The second field of the ImageDescription is Rc+0.5 of that slice,
      the centre convention of ofct_rec (dcnt+0.5); for an axis right of
      the middle this equals the stitched-image centre (N+1)/2 written
      before ver 2.8, so existing output is unchanged.  With a tilted or
      fractional axis the field carries the slice's own value.

      *) rangeList examples:
         Single layer 100: 100
         Layers 100-150: 100-150
         All layers: -

   c. Single Slice Reconstruction (offset CT, img/tif)
      ofct_rec_P_F layer center {pixel size} {offsetangle}

      layer: Layer (height) to reconstruct
      center: Rotation axis position (pixels).  May be fractional
           (ver 2.8): the opposing view is resampled by 2f pixels and
           the rotation centre dcnt+0.5 then sits on the mirror axis.
           Before ver 2.8 a fractional value moved the centre but
           stitched the views at the truncated integer, i.e.
           inconsistently; an integer value gives exactly the old result.
      pixel size: Pixel size (um). Defaults to 1.0 if omitted.
      offset angle: Rotation axis origin offset. Defaults to 0.0 if omitted.

      *) Run in the directory containing q????.img or q????.tif files (auto-detected from dark.img / dark.tif).

   d. Offset-CT Continuous Reconstruction for the Automatic CT System
      (direct 16/8-bit output)
      ofact_srec_g_r HiPic Rc Oy rangeList Dr RA0 rec16 LACmin16 LACmax16 rec8 LACmin8 LACmax8

      Identical to ofct_srec (3b) up to the reconstruction; writes
      rec16/rh%05d.tif (16-bit, LACmin16..LACmax16) and rec8/ro%05d.tif
      (8-bit, LACmin8..LACmax8) directly, with the same quantisation and
      ImageDescription as tif_f2i, so the files are identical to
      "ofct_srec -> tif_f2i".  The three-argument form (HiPic Rc Oy) still
      only prints the memory estimate.  GPU build only; the filter suffix
      is as in ofct_srec (ofact_srec_g_r / _s / _c).

4. Normalization of 32-bit TIFF Images
    tif_f2i bit rec out {LACmin LACmax} {x1 y1 x2 y2}

    bit: Bit depth: 0, 8, or 16 (0 only checks the value range)
    rec: Directory containing 32-bit TIFF reconstructed images
    out: Output directory for normalized images (must be created
         before execution)
    LACmin: Minimum value for normalization (optional)
    LACmax: Maximum value for normalization (optional)
    x1: Horizontal crop start (optional)
    y1: Vertical crop start (optional)
    x2: Horizontal crop end (optional)
    y2: Vertical crop end (optional)
    If LAC min/max are omitted, the min/max values found in rec/ are used
    for normalization.
    The min/max values used for normalization are appended to the end of
    the TIFF tags.
    When the slices differ in size (offset CT changes the reconstructed
    size with the rotation-centre position), all slices are aligned to
    the largest one: smaller slices are centred and the outside is
    filled with pixel value 0 (for an odd size difference the left/top
    padding is one pixel smaller).  The number of aligned slices is
    reported at the end.  The crop coordinates x1..y2 refer to this
    largest frame.

5. Applying Gaussian Filter to 32-bit TIFF Images
    rec_gf rec radius out

    rec: Directory containing 32-bit TIFF reconstructed images
    radius: Half-width of the Gaussian filter
    out: Output directory for filtered images (must be created
         before execution)

6. Miscellaneous Utilities
   a. Display Tags Embedded in TIFF Images
      (pixel size, number of projections, etc.)
      pid tiff-file

   b. Sinogram Generation
      sinog    layer {skip}
      of_sinog layer offset {skip}
      Output as 32-bit TIFF. Tags contain only min and max values.
      If skip is specified, the number of projections is reduced by
      a factor of skip.  of_sinog stitches the two halves of a 360-degree
      offset scan (3) at the rotation axis position offset (as the
      center of ofct_rec; fractional allowed from ver 2.8, an integer
      gives the old result); the axis of the stitched sinogram for 6c is
      offset+0.5.  his_sg (2g) does the same directly from a.his.
      The file sino.tmp (width, projections, layer, centre, mean and SD
      of the total absorption, for an old batch chain) is no longer
      written (ver 2.8).

   c. Reconstruction from Sinogram (generated by 6b)
      sf_rec_t_F input output {Dr RC RA0}
      Output as 32-bit TIFF. CPU only; there is no GPU version.
      input: Sinogram as float TIFF
      output: Output as float TIFF
      Dr: Pixel size
      RC: Rotation axis position
      RA0: Rotation axis origin offset

   d. Average img Images
      img_ave file1 file2... output
      Averages img images. Required when running conv.bat.

   e. Split his File
      spl N-shot N-split
      Likely used for Z-scans or energy scans.
      Splits a his file and saves as img files.
      For data acquired as tif, use act_spl2, which also applies filters (6q).

   f. Renumber rec Sequence
      rec_stk num_stack start end
      Used for Z-scans. Renumbers rec files into a single stack.
      rhro_stk num_stack start end {name ...}
      The same for the 16/8-bit output of act_tg / ofact_srec: prints the
      commands that move <scan>/rh/rh%05d.tif and <scan>/ro/ro%05d.tif,
      layers start..end of scans 001..num_stack, into whole/rh and
      whole/ro with one continuous numbering (layer i of scan j becomes
      i + j*(end-start+1), as rec_stk).  Default names rh ro; give names
      to do only some of them.  Example: rhro_stk 8 17 140 > ddd.sh

   g. Convert his File to Sequential img Files
      his2img his-file (x1 x2 y1 y2)
      Cropping is possible by specifying the arguments in parentheses.

   h. Generate Projection Images from 180-degree Scan
      ct_prj_f HiPic prj {pr.par}
      Input format auto-detected: dark.img -> img, otherwise dark.tif -> tif.
      Output is prj/p%05d.tif (32-bit float).

      With a third argument pr.par, Paganin single-distance phase
      retrieval (Paganin et al., J. Microsc. 206, 33-40, 2002, Eq. 10,
      parallel beam) is applied to every projection and the output value
      is mu*t = -ln(T').  Without it the program does plain absorption
      as before.
      pr.par holds five numbers, one per line (# or // starts a comment):
        MU      linear attenuation coefficient  [cm^-1]  > 0
        DELTA   refractive index decrement      [-]      > 0
        R1      source-sample distance          [cm]     > 0 (unused)
        R2      sample-detector distance        [cm]     >= 0
        P_SIZE  effective pixel size            [um]     > 0
      NOTE: P_SIZE is in micrometres.  Do not reuse pr.par files from
      the old mkpms/mkpmbg chain, which took it in cm; a 4-decade unit
      mistake over- or under-smooths drastically.  Extreme values stop
      with an error at start-up, but not every mistake can be caught, so
      check that the rc printed at start-up (smoothing radius in pixels)
      is physically reasonable - normally a few pixels.
      The low-transmission guard (1e) applies to the filtered
      transmission as well.  Thread count: environment variable
      PAGANIN_THREADS (default: the number of logical CPUs).

   i. Histogram from 32-bit TIFF
      tif2hst rec (x1 y1 x2 y2)
      rec: Directory containing 32-bit TIFF reconstructed images
      x1: Horizontal crop start (optional)
      y1: Vertical crop start (optional)
      x2: Horizontal crop end (optional)
      y2: Vertical crop end (optional)

   j. Orthogonal Rotation of 8-bit or 16-bit CT Images
      si_rar orgDir nameFile hAxis vAxis dAxis newDir
      si_rar orgDir nameFile hAxis vAxis dAxis sliceNo newTIFF

      Example: si_rar.exe ro_xy - +y +z +x ro_yz
      Example: si_rar.exe ro_xy - +z +x +y ro_zx

   k. Binning of CT Images
      For 8-bit or 16-bit images:
      si_sir orgDir nameFile Bxyz newDir
      si_sir orgDir nameFile Bx By Bz newDir

      Example: si_sir ro_xy - 2 ro_2x2x2
      Example: si_sir ro_xy - 2 2 1 ro_2x2x1

      For 32-bit TIFF (rec?????.tif):
      rec_sir in out B

      in: Directory containing 32-bit TIFF reconstructed images
      out: Output directory for the binned images (must be created before
           execution)
      B: Binning factor, applied to x, y and z (the slice direction) alike

      Each B x B x B block of voxels is replaced by its average and written
      as 32-bit float (no rounding to an integer). The number of slices and
      both in-plane dimensions become 1/B.
      The first output file takes the number of the first input file and the
      rest follow one by one: rec00003 through rec00011 with B=3 gives
      rec00003, rec00004 and rec00005.
      Columns, rows and slices that do not fill a block are dropped; this is
      the only behavioural difference from si_sir, which keeps them and
      divides by the number of voxels actually present.
      The TIFF tags are updated for the new grid: the pixel size becomes B
      times larger, the rotation axis position becomes (RC-(B-1)/2)/B, the
      number of projections and the rotation angle offset are carried over,
      and the min/max values are recomputed from the result.
      Because the output numbers overlap the input numbers, giving the same
      directory for both input and output stops with an error.

      Example: rec_sir rec rec_b2 2

   l. Gaussian Filter for 8-bit or 16-bit CT Images
      si_gf orgDir nameFile radius {bias} newDir

      Example: si_gf rh - 1.0 rh_gf1
      Example: si_gf rh - 1.0 0.0 rh_gf1

   m. Automatic Z-scan Stacking (Example)
      Number of shots per scan: 7501
      Number of Z-stacks: 8
      Overlap location: 17 140
      Pixel size: 11.31

      spl 7501 8 > aaa.bat
      call aaa.bat
      chk-rc 8 75 > bbb.bat
      call bbb.bat
      set-rc 11.31 > ccc.bat
      call ccc.bat
      rec_stk 8 17 140 > ddd.bat
      call ddd.bat

      Directly from a.his (2g), without spl and conv.bat (the number of
      scans is read from a.his; the layer defaults to the centre row):
      rc-check-his 75
      gen-all-his center.log 11.31
      rec_stk 8 17 140 > ddd.bat
      call ddd.bat

   n. Cropping 32-bit TIFF Images
      rec_crop in out {x1 y1 x2 y2}

      in: Directory containing 32-bit TIFF reconstructed images
      out: Output directory for cropped images (must be created
           before execution)
      x1: Horizontal crop start
      y1: Vertical crop start
      x2: Horizontal crop end
      y2: Vertical crop end

   o. Applying Filters to CT Images of Any Bit Depth
      Available filters: Bilateral filter, Gaussian filter, Median filter.
      Executables ending with _g are GPU versions (check your CUDA Toolkit
      version).
      tif_blf, tif_blf_g, tif_gsf, tif_gsf_g, tif_mdf, tif_mdf_g

      Bilateral filter
      tif_blf[_g] <input_dir> <output_dir> [kernel_size] [spatial_sigma] [intensity_sigma]

      Gaussian filter
      tif_gsf[_g] <input_dir> <output_dir> [sigma]

      Median filter
      tif_md[_g]f <input_dir> <output_dir> [kernel_size]

      Parameters in [] are automatically determined if omitted.

      Example: tif_blf_g rh rh_blf1
      Example: tif_blf rh rh_blf1 5 20 200

      See 20260504_filter_readme.pdf for details on additional filters.

   p. Apply Median Filter Followed by Gaussian Filter to TIFF Images
      tif_mgf[_g] <input_file> <output_file> [median_kernel_size] [gaussian_sigma]
      Used to remove noise such as scattered light in X-ray transmission
      images. _g is the GPU version.

      Example: tif_mgf a000401.tif 001/raw/a0401.tif 3 1
      Example: tif_mgf a000401.tif 001/raw/a0401.tif 0 1
      Example: tif_mgf_g a000401.tif 001/raw/a0401.tif 3 1

   q. Split and Save Continuously Acquired tif Data (with Filtering)
      act_spl2 N-shot N-split {m_kernel_size} {g_kernel_size}
      act_spl  N-shot N-split {kernel_size}

      When a??????.tif files are acquired continuously in a Z-scan or an
      energy scan, these tools distribute them one scan at a time into
      001/raw, 002/raw, ... They are the tif counterpart of spl (6e), which
      splits a his file into img files.
      N-shot: Number of shots per scan
      N-split: Number of splits (number of Z-stacks)
      act_spl2 applies tif_mgf (median + gaussian) while distributing:
        m_kernel_size: Median kernel size (0 if omitted, i.e. no median)
        g_kernel_size: Gaussian sigma (0 if omitted, i.e. no gaussian)
      act_spl applies only a gaussian via gf_sd, and simply moves the files
      without filtering when kernel_size is omitted.
      Both copy conv.bat and output.log from the working directory into every
      raw/, rewriting img to tif inside conv.bat, and append a record to
      cmd-hst.log when finished.
      The filters are run in the background: with start /b on Windows and with
      fork+exec on Linux. On Linux the number of concurrent jobs is limited by
      the environment variable ACT_SPL_JOBS (default 16), and all jobs are
      waited for before the program exits.

      Example: act_spl2 7501 8 3 1
      Example: act_spl  7501 8

   r. Re-projection / Re-reconstruction of Finished CT Slices
      rec2rec_F   rec_in rec_out {Nt}      (CPU)
      rec2rec_g_F rec_in rec_out {Nt}      (GPU)

      rec_in: Directory containing 32-bit float CT slices (rec*.tif, square)
      rec_out: Output directory (same file names as the input)
      Nt: Number of views of the forward projection. If omitted, the number of
          projections recorded in the input TIFF tag is used, falling back to
          the image size N when it cannot be read. Specifying 0 also selects
          the automatic choice.

      Each input slice is forward Radon-transformed into a sinogram and then
      reconstructed again with ring removal + CBP. No scaling is applied, so
      the input values come out as they are. When combining this with the
      truncation correction of 1f, lower PAD_THRESH to about 0.1-0.2.

      Example: rec2rec_g_c rec rec2 1800

   s. Numbered TIFF to OME-Zarr Conversion
      tif2zar tifDir out.zarr {--chunk N} {--pixel P} {--min A --max B}
                              {--levels L} {--codec C} {--clevel N}

      Converts a reconstructed volume (numbered uint16 / 32-bit float
      TIFFs) to OME-Zarr v0.4 (multiscale pyramid + blosc-zstd).  Open it
      in Fiji via Plugins -> BigDataViewer -> N5 Viewer with a
      file:///D:/... URI (typing a bare C:\ path trips a known Fiji bug).
      float32 input is quantised to uint16 using the global range taken
      from the per-slice description tags; for rec data that range is
      dominated by outliers, so passing the same normalisation values as
      tif_f2i via --min/--max is recommended.  Gaps in the series stop
      with an error listing every missing file.  Slices of different
      sizes (offset CT) are centred in the largest frame and zero-padded
      (same rule as tif_f2i, section 4).  Thread count: the
      environment variable TIF2ZAR_THREADS (default: logical CPUs).
      Details and measured performance: 20260909_tif2zar.md (Japanese).

      example: tif2zar rec rec.zarr --min -0.5 --max 3.0
