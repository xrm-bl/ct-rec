Description of Image Reconstruction Software and Normalization
Based on Nakano's Software

Uesugi

2026.07.02  ver. 2.2
2026.06.30  ver. 2.1
2026.06.30  ver. 2.0
2026.05.04  ver. 1.7

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
      cmd-hst.log upon completion.

   c. Program Suffixes
      The reconstruction software has suffixes such as _t_c.
      These specify the processor and reconstruction filter to use.
      _P: Processor
          _t: Use CPU multi-threading. Controlled by the environment
              variable CBP_THREADS. Default is 8 threads.
          _g: Use GPGPU. The included executables are compiled with
              CUDA Toolkit 10.2.
      _F: Filter
          _c: Chesler filter
          _s: Shepp-Logan filter
          _r: Ramachandran (HAN) filter

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
      CBP_THREADS (for back-projection computation, default 8) described
      in section 1c.

   e. Missing Angle Handling
      For plate-like samples, transmittance can drop drastically at certain
      angles. In extreme cases, some processing is necessary. The "degree"
      of this processing can be controlled.
      Specify using the environment variable CT_REC_BLACK_THRESH (default
      is 1 if not set).
      When missing angles are present, change this value to 1, 10, 100,
      1000, etc. and run ct_rec to adjust the reconstruction behavior.

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

3. 360-degree Scan (Offset CT): Standard Absorption CT Reconstruction
   a. Rotation Axis Position Estimation
      ofct_DO   raw     (CPU)
      ofct_DO_g raw     (GPU; same result as ofct_DO)

      raw: Directory containing q????.img or q????.tif files (no trailing /;
           input format auto-detected from dark.img / dark.tif)

      Estimates the rotation-axis position from the offset-CT data and prints
      a ready-to-run ofct_srec command (suggested center and Oy).
      ofct_DO_g computes the per-view-pair MSD on the GPU (host-side
      processing is identical, so the estimate matches the CPU version).

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
      sf_rec_P_F input output {Dr RC RA0}
      Output as 32-bit TIFF.
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

   k. Binning of 8-bit or 16-bit CT Images
      si_sir orgDir nameFile Bxyz newDir
      si_sir orgDir nameFile Bx By Bz newDir

      Example: si_sir ro_xy - 2 ro_2x2x2
      Example: si_sir ro_xy - 2 2 1 ro_2x2x1

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

      See filter_readme.pdf for details on additional filters.

   p. Apply Median Filter Followed by Gaussian Filter to TIFF Images
      tif_mgf <input_file> <output_file> [median_kernel_size] [gaussian_sigma]
      Used to remove noise such as scattered light in X-ray transmission
      images.

      Example: tif_mgf a000401.tif 001/raw/a0401.tif 3 1
      Example: tif_mgf a000401.tif 001/raw/a0401.tif 0 1
