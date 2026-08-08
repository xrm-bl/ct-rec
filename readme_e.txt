Description of Image Reconstruction Software and Normalization
Based on Nakano's Software

Uesugi

2026.08.06  ver. 2.4
2026.07.30  ver. 2.3
2026.07.02  ver. 2.2
2026.06.30  ver. 2.1
2026.06.30  ver. 2.0
2026.05.04  ver. 1.7

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
      This applies to ct_rec / hp_tg / p_rec / ofct_rec / ofct_srec / sf_rec /
      rec2rec alike. Note that the input of rec2rec is already attenuated at
      the edges by the reconstruction circle, so a smaller value of about
      0.1-0.2 is appropriate there.

   g. Other Environment Variables
      Besides the above there are the number of parallel projection-read
      threads (HPTG_READ_THREADS, default 16; effective in hp_tg / ofct_srec),
      the memory-based chunking controls (HPTG_MEM_FRACTION /
      HPTG_MEM_LIMIT_MB / HPTG_CHUNK_ROWS), the CUDA device number to use
      (CUDA_GPU, default 0) and the input file name overrides (RHP_O / RHP_D /
      RHP_Q). See 20260723_CT_env_vars.md for the full list with defaults.
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

      HiPic: Directory containing q????.img or q????.tif files (no trailing /)
      Rc: Rotation axis position (pixels from left edge)
      Oy: Vertical shift (always specify 0)
      rangeList: Specify layers to reconstruct.
      Dr: Pixel size
      RA0: Rotation axis origin offset
      rec: Output directory for reconstructed images (must be created
           before execution)

      *) rangeList examples:
         Single layer 100: 100
         Layers 100-150: 100-150
         All layers: -

   c. Single Slice Reconstruction (offset CT, img/tif)
      ofct_rec_P_F layer center {pixel size} {offsetangle}

      layer: Layer (height) to reconstruct
      center: Rotation axis position (pixels)
      pixel size: Pixel size (um). Defaults to 1.0 if omitted.
      offset angle: Rotation axis origin offset. Defaults to 0.0 if omitted.

      *) Run in the directory containing q????.img or q????.tif files (auto-detected from dark.img / dark.tif).


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
      sinog layer {skip}
      Output as 32-bit TIFF. Tags contain only min and max values.
      If skip is specified, the number of projections is reduced by
      a factor of skip.

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

   g. Convert his File to Sequential img Files
      his2img his-file (x1 x2 y1 y2)
      Cropping is possible by specifying the arguments in parentheses.

   h. Generate Projection Images from 180-degree Scan
      ct_prj_f HiPic prj
      Input format auto-detected: dark.img -> img, otherwise dark.tif -> tif.

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
      the environment variable ACT_SPL_JOBS (default 8), and all jobs are
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

