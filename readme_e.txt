Explanation of image reconstruction software based on Nakano's software and around normalization

K. Uesugi

2025.05.04  ver. 1.4

0. If you find any bugs or requests, please contact the author.

1. Basic rules
   a. Input or raw image files
      img type (ITEX made by Hamamatsu) and tiff monochrome.

   b. Output
      32bit tiff for CT images as rec?????.tif (5 digits).
      Each tiff tag contains the following information
        pixel size, position of the rotation axis, number of projections,
        rotation angle offset, minimum and maximum values
      They are embedded in that order.
      When the image is normalized, the minimum and maximum values at the time
      of normalization are added to the previous tags. For continuous
      reconstruction and normalization, a log is left in cmd-hst.log
      when the execution is completed.

   c. Subscript of a program
      The image reconstruction software has subscripts such as _t_c.
      These specify the processor and reconstruction filter to be used for the operation.
      _P: processor
          _t: Use the multi-threading function of the CPU. Controlled by the environment
              variable CBP_THREADS. Default is 8 threads.
          _g: GPGPU is used. The attached exe has been compiled with CUDA toolkit 11.2.
      _F: Filter function
          _c: Chesler Filter
          _s: Shepp-Logan Filter
          _r: Ramachandran(HAN) Filter

   d. Ring Artifact Removal
      This version 1.4 implements a ring removal function based on Algorithm 3 from Vo et al. (2018).
      The function is executed immediately before CBP calculation. This feature can be turned
      ON/OFF using environment variables. Setting the KERNEL_SIZE environment variable to 1
      turns it OFF. Using other positive odd numbers changes its effect. The default value
      is set to 5 (it also defaults to 5 if the environment variable is not defined).
      Additionally, CPU parallel processing is implemented, with OMP_NUM_THREADS defaulting to 40.
      For information on how to set these environment variables, please refer to the instructions
      at the beginning of the "sort_filter_omp.c".

2. 180deg scan. image reconstruction of standard absorption contrast CT.

   a. Reconstruction one slice
      ct_rec_P_F layer {center} {pixel size} {offsetangle}
      
      layer: Layer to be reconstructed (height)
      center: Position of the rotation axis (pixel). If omitted, it will be estimated automatically.
      pixel: size: Position of the rotation axis (pixel). If omitted, it will be estimated automatically.
      offset angle: Origin offset of the rotation axis. If omitted, the value is 0.0.
      
      *) Run it in the directory where q????.img is located.
   
   b. Reconstruction one slice
      tf_rec_P_F layer {center} {pixel size} {offsetangle}
      
      layer: Layer to be reconstructed (height)
      center: Position of the rotation axis (pixel). If omitted, it will be estimated automatically.
      pixel: size: Position of the rotation axis (pixel). If omitted, it will be estimated automatically.
      offset angle: Origin offset of the rotation axis. If omitted, the value is 0.0.
      
      *) Run it in the directory where q????.tif is located.
   
   c. Continuous reconstruction
      hp_tg_P_F HiPic/ Dr RC RA0 rec/
      tf_tg_P_F HiPic/ Dr RC RA0 rec/
      (When the rotation axis is not tilted. All layers.)
      
      HiPic/: The name of the directory in which the q????.img or q????.tif is stored.(/ is not required)
      Dr: pixel size (um)
      RC: position of rotation axis
      RA0: offset angle of rotation
      rec/: Directory for outputting reconstructed images(To be created before calculation)
      
      hp_tg_P_F HiPic/ Dr L1 C1 L2 C2 RA0 rec/
      tf_tg_P_F HiPic/ Dr L1 C1 L2 C2 RA0 rec/
      (If the rotation axis is tilted. or when calculating only a part of the area)
      
      HiPic/: The name of the directory in which the q????.img or q????.tif is stored.(/ is not required)
      Dr: pixel size (um)
      L1: Calculation start layer
      C1: position of rotation axis at L1
      L2: Calculation end layer
      C2: position of rotation axis at L2
      RA0: offset angle of rotation
      rec/: Directory for outputting reconstructed images(/ is not required。To be created before calculation)
   
      *) Run it on one of the directories where q???? .img is located.

   d. Continuous reconstruction from p-images
      p_rec_P_F p/ rec/ Dr RC RA0
      (When the rotation axis is not tilted. All layers.)
      
      p/: The name of the directory in which the p????.tif is stored.(/ is not required)
      rec/: Directory for outputting reconstructed images(To be created before calculation)
      Dr: pixel size (um)
      RC: position of rotation axis
      RA0: offset angle of rotation
      
      p_rec_P_F p/ rec/ Dr L1 C1 L2 C2 RA0
      (If the rotation axis is tilted.)
      
      p/: The name of the directory in which the p????.tif is stored.(/ is not required)
      rec/: Directory for outputting reconstructed images(To be created before calculation)
      Dr: pixel size (um)
      L1: Calculation start layer
      C1: position of rotation axis at L1
      L2: Calculation end layer
      C2: position of rotation axis at L2
      RA0: offset angle of rotation

3. 360deg scan (offset CT)。Image reconstruction of standard absorption contrast.
   a. Estimation of rotational axis position
      ofct_xy HiPic/ {Ox1 Ox2 Oy1 Oy2} {MSD.tif}
      oftf_xy HiPic/ {Ox1 Ox2 Oy1 Oy2} {MSD.tif}
      
      HiPic/: The name of the directory in which the q????.img or q????.tif is stored.(/ is not required)
      Ox1: Starting point of the horizontal search area. (optional)
      Ox2: End point of the horizontal search area. (optional)
      Oy1: The starting point of the vertical search area. (optional)
      Oy2: End point of vertical search range. (optional)
      MSD.tif: Name of the MSD image output file. (optional)

      *) When outputting the image, the offset of the image is displayed in parentheses.

   b. Check the amount of memory required.
      ofct_srec_P_F HiPic/ Rc Oy
      oftf_srec_P_F HiPic/ Rc Oy

      HiPic/: The name of the directory in which the q????.img or q????.tif is stored.(/ is not required)
      Rc: position of rotation axis (Number of pixels from the left edge)
      Oy: amount of vertical misalignment (Always set to 0).

      This will output the amount of main memory needed to calculate everything at once.
      It can be used as a guide to determine the rangeList.

   c. Reconstruction of offset CT
      ofct_srec_P_F HiPic/ Rc Oy rangeList Dr RA0 rec/
      oftf_srec_P_F HiPic/ Rc Oy rangeList Dr RA0 rec/
      
      HiPic/: The name of the directory in which the q????.img or q????.tif is stored.(/ is not required)
      Rc: position of rotation axis(Number of pixels from the left edge)
      Oy: amount of vertical misalignment (Always set to 0).
      rangeList: Specify the layer to be reconstructed.
      Dr: pixel size
      RA0: offset angle of rotation
      rec/: Directory for outputting reconstructed images(To be created before calculation)

      *) rangeListのEx.。
         For only 100 layers：100
         100-150 layers: 100-150
         All layers: -

4. 32bit tiff image normalization.
    tif_f2i bit rec/ out/ {LACmin LACmax} {x1 y1 x2 y2}
    
    bit: Number of bits. 0 or 8 or 16 (0 is for checking the range.)
    rec/: Directory with 32bit tiff reconstructed images
    out/: Directory for outputting the normalized images (to be created before execution)
    LACmin: The minimum value of normalization. (optional)
    LACmax: The maximum value of normalization. (optional)
    x1: Beginning of lateral cutout. (optional)
    y1: Beginning of the longitudinal cutout. (optional)
    x2: End of the lateral cutout. (optional)
    y2: End of the longitudinal cutout. (optional)
    When the maximum and minimum values of LAC are omitted, the maximum and minimum values
    in rec/ shall be used for normalization. The maximum and minimum values at the time of
    normalization are added to the end of the tiff tag.

5. Applying gaussian filter to 32bit tiff CT image.
    rec_gf rec/ radius out/
    
    rec/: Directory with 32bit tiff reconstructed images
    radius: Radius of gaussian filter.
    out/: Output directory (make beforehand)

6. Others
   a. print image description
      pid tiff-file
      Show tiff tags.
      
   b. make sinogram
      sinog layer {skip}
      32bit tiff output. Maximum and minimum values are described in tags.
      If "skip" is specified, the number of projections is set to 1/skip and output.

   c. Reconstructed from the synogram of b
      sf_rec_P_F input output {Dr RC RA0}
      32bit tiff output.
      input: file name of sinogram image with 32bit float tiff
      output: file name of output CT image
      Dr: pixel size
      RC: position of rotation axis
      RA0: offset angle of rotation

   d. averaging of img images
      img_ave file1 file 2... output

   e. split his file to img files
      spl N-shot N-split
      N-shot: number of shot for a data set
      N-split: number of data set

   f. Re-numbering "rec" files
      rec_stk num_stack start end
      num_stack: number of stack data set
      start: 
      end: 
      Used for Z-scans, to renumber the recs to make a single stack.

   g. Convert his file to a sequentially numbered img.
      his2img his-file (x1 x2 y1 y2)
      crop is available when specified.
      
   h. Create projection images from a 180deg scan
      ct_prj_f HiPic/ prj/
      HiPic/: The name of the directory in which the q????.img is stored.(/ is not required)
      prj/: Directory for projection images (make beforehand)

   i. Histogram from 32bit tiff
      tif2hst rec/ (x1 y1 x2 y2)
      rec/: Directory for 32bit tiff CT images
      x1: Beginning of lateral cutout. (optional)
      y1: Beginning of the longitudinal cutout. (optional)
      x2: End of the lateral cutout. (optional)
      y2: End of the longitudinal cutout. (optional)

   j. Orthogonal rotation of 8bit or 16bit CT image. 
      si_rar orgDir nameFile hAxis vAxis dAxis newDir
      si_rar orgDir nameFile hAxis vAxis dAxis sliceNo newTIFF

      Ex.: si_rar.exe ro_xy - +y +z +x ro_yz
      Ex.: si_rar.exe ro_xy - +z +x +y ro_zx

   k. Binning of 8bit or 16bit CT image.
      si_sir orgDir nameFile Bxyz newDir
      si_sir orgDir nameFile Bx By Bz newDir

      Ex.: si_sir ro_xy - 2 ro_2x2x2
      Ex.: si_sir ro_xy - 2 2 1 ro_2x2x1

   l. Applying Gaussian filter to 8bit or 16bit CT image.
      si_gf orgDir nameFile radius {bias} newDir

      Ex.: si_gf rh - 1.0 rh_gf1
      Ex.: si_gf rh - 1.0 0.0 rh_gf1

   m. Example of auto stack for Z-scan
      Number of images for a scan: 7501
      Number of Z stack: 8
      Overlapping layer: 17 140
      pixel size: 11.31
      
      spl 7501 8 > aaa.bat
      call aaa.bat
      chk-rc 8 75 > bbb.bat
      call bbb.bat
      set-rc 11.31 > ccc.bat
      call ccc.bat
      rec_stk 8 17 140 > ddd.bat
      call ddd.bat

   n. cropping 32bit tiff files
      rec_crop in/ out/ {x1 y1 x2 y2}
    
      in/:  input 32bit tiff files directory
      out/: output directory
      x1: Beginning of lateral cutout.
      y1: Beginning of the longitudinal cutout.
      x2: End of the lateral cutout.
      y2: End of the longitudinal cutout.

