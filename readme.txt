����\�t�g���x�[�X�ɂ����摜�č\���\�t�g�ƋK�i���ӂ�̐���

�㐙

2025.05.04  ver. 1.4

0. �o�O�⃊�N�G�X�g�͍�҂ɘA�����Ă��������B

1. ���ʂ̍l��
   a. ����
      ��{�I��img�`���B

   b. �o��
      CT���� 32bit tiff �ł̏o�͂ƂȂ�Brec?????.tif (���l��5�P�^)
      ���ꂼ��� tiff �^�O�ɂ́A
      ��f�T�C�Y�E��]���̈ʒu�E���e���E��]�p�I�t�Z�b�g�E���̉摜�ł̍ŏ��E�ő�l
      �����̏��Ԃɖ��ߍ��܂�Ă���B
      �K�i�����ɂ́A����܂ł̃^�O�ɉ����A�K�i�����̍ŏ��ő�l���ǋL�����B
      �A���č\���ƋK�i���Ɋւ��ẮA���s�������� cmd-hst.log �Ƀ��O���c��

   c. �v���O�����̓Y����
      �摜�č\���\�t�g�ɂ� _t_c �Ȃǂ̓Y���������Ă���B
      �����͉��Z�Ɏg�p����v���Z�b�T�[�ƍč\���t�B���^�[���w�肷��B
      _P: �v���Z�b�T�[
          _t: CPU �̃}���`�X���b�h�@�\���g�p����B���ϐ�CBP_THREADS�Ő���B�f�t�H���g�͂W�X���B
          _g: GPGPU ���g�p�B�Y�t��exe�� CUDA toolkit 10.2 �ŃR���p�C���ς݁B
      _F: �t�B���^�[
          _c: Chesler �t�B���^�[
          _s: Shepp-Logan �t�B���^�[
          _r: Ramachandran(HAN)�t�B���^�[
      �ƂȂ��Ă���B

   d. �����O�A�[�e�B�t�@�N�g�̏���
      ���̃o�[�W��������Vo et al.(2018)��Argorism 3�^�̃����O�����@�\��݂����B
      CBP�v�Z�̒��O�Ɏ��s���Ă���B���ϐ����w�肷�邱�Ƃł��̋@�\��ON/OFF�ł���B
      ���ϐ��� KERNEL_SIZE ��1�Ɏw�肷���OFF�B����ȊO�̐��̊�Ō��ʂ��ς��B
      �f�t�H���g�l��5�Ƃ��Ă���(���ϐ�����`����Ă��Ȃ��ꍇ��5�ɂȂ�)�B
      �܂��ACPU����v�Z���s���Ă���A�f�t�H���g�l��OMP_NUM_THREADS��40�Ƃ��Ă���B
      ���ϐ��̐ݒ���@�́Asort_filter_omp.c �̖`���ɋL�q������̂ŎQ�Ƃ̂��ƁB

2. 180deg scan�B�W���I�ȋz���̉摜�č\���B

   a. 1�������č\��
      ct_rec_P_F layer {center} {pixel size} {offsetangle}
      
      layer: �č\�����郌�C���[(����)
      center: ��]���̈ʒu(pixel)�B�ȗ������ꍇ�͎������肷��B
      pixel: size: ��f�T�C�Y(um)�B�ȗ������ꍇ��1.0�ɂȂ�B
      offset angle: ��]���̌��_�I�t�Z�b�g�B�ȗ������ꍇ��0.0�ɂȂ�B
      
      *) q????.img ������f�B���N�g���Ŏ��s����B
   
   b. 1�������č\��
      tf_rec_P_F layer {center} {pixel size} {offsetangle}
      
      layer: �č\�����郌�C���[(����)
      center: ��]���̈ʒu(pixel)�B�ȗ������ꍇ�͎������肷��B
      pixel: size: ��f�T�C�Y(um)�B�ȗ������ꍇ��1.0�ɂȂ�B
      offset angle: ��]���̌��_�I�t�Z�b�g�B�ȗ������ꍇ��0.0�ɂȂ�B
      
      *) q????.tif ������f�B���N�g���Ŏ��s����B
   
   c. �A���č\��
      hp_tg_P_F HiPic/ Dr RC RA0 rec/
      (��]�����X���ĂȂ��ꍇ�B�S���C���[)
      
      HiPic/: q????.img ���i�[����Ă���f�B���N�g�����B(/ �͕s�v)
      Dr: ��f�T�C�Y (um)
      RC: ��]���̈ʒu
      RA0: ��]���̌��_�I�t�Z�b�g
      rec/: �č\���摜���o�͂���f�B���N�g��(�v�Z�O�ɍ쐬���邱��)
      
      hp_tg_P_F HiPic/ Dr L1 C1 L2 C2 RA0 rec/
      (��]�����X���Ă���ꍇ�B�������͈ꕔ�̗̈�݂̂̌v�Z��)
      
      HiPic/: q????.img ���i�[����Ă���f�B���N�g�����B(/ �͕s�v)
      Dr: ��f�T�C�Y (um)
      L1: �v�Z�J�n���C���[
      C1: L1�ł̉�]���̈ʒu
      L2: �v�Z�I�����C���[
      C2: L2�ł̉�]���̈ʒu
      RA0: ��]���̌��_�I�t�Z�b�g
      rec/: �č\���摜���o�͂���f�B���N�g��(/ �͕s�v�B�v�Z�O�ɍ쐬���邱��)
   
      *) q????.img ������f�B���N�g���̈��Ŏ��s����B

   d. p�摜����̘A���č\��
      p_rec_P_F p/ rec/ Dr RC RA0
      (��]�����X���ĂȂ��ꍇ�B�S���C���[)
      
      p/: p?????.tif ���i�[����Ă���f�B���N�g�����B(/ �͕s�v)
      rec/: �č\���摜���o�͂���f�B���N�g��(�v�Z�O�ɍ쐬���邱��)
      Dr: ��f�T�C�Y (um)
      RC: ��]���̈ʒu
      RA0: ��]���̌��_�I�t�Z�b�g
      
      p_rec_P_F p/ rec/ Dr L1 C1 L2 C2 RA0
      (��]�����X���Ă���ꍇ)
      
      p/: p?????.tif ���i�[����Ă���f�B���N�g�����B(/ �͕s�v)
      rec/: �č\���摜���o�͂���f�B���N�g��(�v�Z�O�ɍ쐬���邱��)
      Dr: ��f�T�C�Y (um)
      L1: �v�Z�J�n���C���[
      C1: L1�ł̉�]���̈ʒu
      L2: �v�Z�I�����C���[
      C2: L2�ł̉�]���̈ʒu
      RA0: ��]���̌��_�I�t�Z�b�g

3. 360deg scan (offset CT)�B�W���I�ȋz���̉摜�č\��
   a. ��]���ʒu�̐���
      ofct_xy HiPic/ {Ox1 Ox2 Oy1 Oy2} {MSD.tif}
      
      HiPic/: q????.img ���i�[����Ă���f�B���N�g����(/ �͕s�v)
      Ox1: ���̑{���͈͊J�n�_�B(�ȗ���)
      Ox2: ���̑{���͈͏I���_�B(�ȗ���)
      Oy1: �c�̑{���͈͊J�n�_�B(�ȗ���)
      Oy2: �c�̑{���͈͏I���_�B(�ȗ���)
      MSD.tif: MSD���摜�o�̓t�@�C�����B(�ȗ���)

      *) �摜�o�͎��ɂ̓J�b�R���ɉ摜�̃I�t�Z�b�g���\�������B

   b. �K�v�������ʂ̊m�F
      ofct_srec_P_F HiPic/ Rc Oy

      HiPic/: q????.img ���i�[����Ă���f�B���N�g����(/ �͕s�v)
      Rc: ��]���̈ʒu(���[����̉�f��)
      Oy: �c����ʁi��� 0 ���w��j�B

      ����őS������x�Ɍv�Z���邽�߂ɕK�v�ȃ��C���������̗ʂ��o�͂����B
      rangeList�����߂�ۂ̖ڈ��ɂȂ�B

   c. �č\��
      ofct_srec_P_F HiPic/ Rc Oy rangeList Dr RA0 rec/
      
      HiPic/: q????.img ���i�[����Ă���f�B���N�g����(/ �͕s�v)
      Rc: ��]���̈ʒu(���[����̉�f��)
      Oy: �c����ʁi��� 0 ���w��j�B
      rangeList: �č\�����郌�C���[���w��B
      Dr: ��f�T�C�Y
      RA0: ��]���̌��_�I�t�Z�b�g
      rec/: �č\���摜���o�͂���f�B���N�g��(�v�Z�O�ɍ쐬���邱��)

      *) rangeList�̗�B
         100���C���[�����̏ꍇ�F100
         100-150 ���C���[: 100-150
         �S��: -
   d. tiff �f�[�^����1�������č\��
      otf_rec_P_F layer center {pixel size} {offsetangle}
      
      layer: �č\�����郌�C���[(����)
      center: ��]���̈ʒu(pixel)�B
      pixel: size: ��f�T�C�Y(um)�B�ȗ������ꍇ��1.0�ɂȂ�B
      offset angle: ��]���̌��_�I�t�Z�b�g�B�ȗ������ꍇ��0.0�ɂȂ�B
      
      *) q????.tif ������f�B���N�g���Ŏ��s����B
   


4. 32bit tiff �摜�̋K�i��
    tif_f2i bit rec/ out/ {LACmin LACmax} {x1 y1 x2 y2}
    
    bit: �r�b�g�� 0 �� 8 �� 16 (0 �͔͈͂𒲂ׂ邾��)
    rec/: 32bit tiff �č\���摜������f�B���N�g��
    out/: �K�i����̉摜���o�͂���f�B���N�g��(���s�O�ɍ쐬���邱��)
    LACmin: �K�i�����̍ŏ��l�B(�ȗ���)
    LACmax: �K�i�����̍ő�l�B(�ȗ���)
    x1: �������؂�o���̎n�܂�B(�ȗ���)
    x2: �������؂�o���̏I���B(�ȗ���)
    y1: �c�����؂�o���̎n�܂�B(�ȗ���)
    y2: �c�����؂�o���̏I���B(�ȗ���)
    LAC �̍ő�ŏ����ȗ������ꍇ�́Arec/ ���̍ő�ŏ��l��p���ċK�i������B
    �K�i�����̍ő�ŏ��l��tiff�̃^�O�̍Ō�ɕt���������Ă���B

5. 32bit tiff �摜�ւ�gaussian filter�K�p
    rec_gf rec/ radius out/
    
    rec/: 32bit tiff �č\���摜������f�B���N�g��
    radius: gaussian filter �̔��l��
    out/: �t�B���^�[��̉摜���o�͂���f�B���N�g��(���s�O�ɍ쐬���邱��)

6. ���̑����܂�
   a. print image description
      pid tiff-file
      tif�̃^�O��\������
      
   b. �V�m�O�����쐬
      sinog layer {skip}
      32bit tiff �ɂďo�́B�^�O�ɂ͍ő�ƍŏ��l�̂ݓ���B
      skip ���w�肷��ƁA���e����skip����1�ɂ��ďo�͂���B

   c. b �̃V�m�O��������č\��
      sf_rec_P_F input output {Dr RC RA0}
      32bit tiff �ɂďo�́B
      input: float tiff �ŏo�����V�m�O����
      output: float tiff �ŏo��
      Dr: ��f�T�C�Y
      RC: ��]���̈ʒu
      RA0: ��]���̌��_�I�t�Z�b�g

   d. img �摜����
      img_ave file1 file 2... output
      img �摜�𕽋ω�����Bconv.bat ���s���ɕK�{�B

   e. his �����ۑ�
      spl N-shot N-split
      Z�X�L������G�l���M�[�X�L�����������Ƃ��Ɏg���\���������B
      his �t�@�C���𕪊����� img �ɕۑ��B

   f. rec �A�ԐU��Ȃ����B
      rec_stk num_stack start end
      Z�X�L�����������Ƃ��Ɏg���Brec �̘A�Ԃ�U��Ȃ�����1�̃X�^�b�N�ɂ���B

   g. his �t�@�C����A�Ԃ� img �ɕϊ��B
      his2img his-file (x1 x2 y1 y2)
      �J�b�R���̈����w��ɂ�crop���\
      
   h. 180deg �X�L�������瓊�e�����쐬
      ct_prj_f HiPic/ prj/

   i. 32bit tiff ����q�X�g�O����
      tif2hst rec/ (x1 y1 x2 y2)
      rec/: 32bit tiff �č\���摜������f�B���N�g��
      x1: �������؂�o���̎n�܂�B(�ȗ���)
      y1: �c�����؂�o���̎n�܂�B(�ȗ���)
      x2: �������؂�o���̏I���B(�ȗ���)
      y2: �c�����؂�o���̏I���B(�ȗ���)

   j. 8bit or 16bit CT���̒�����]
      si_rar orgDir nameFile hAxis vAxis dAxis newDir
      si_rar orgDir nameFile hAxis vAxis dAxis sliceNo newTIFF

      ��: si_rar.exe ro_xy - +y +z +x ro_yz
      ��: si_rar.exe ro_xy - +z +x +y ro_zx

   k. 8bit or 16bit CT���̃r�j���O
      si_sir orgDir nameFile Bxyz newDir
      si_sir orgDir nameFile Bx By Bz newDir

      ��: si_sir ro_xy - 2 ro_2x2x2
      ��: si_sir ro_xy - 2 2 1 ro_2x2x1

   l. 8bit or 16bit CT����gaussian filter
      si_gf orgDir nameFile radius {bias} newDir

      ��: si_gf rh - 1.0 rh_gf1
      ��: si_gf rh - 1.0 0.0 rh_gf1

   m. Z�X�L�����̎����X�^�b�N(��)
      1��̎B�e����: 7501
      Z�X�^�b�N�̐��F8
      �I�[�o�[���b�v�̏ꏊ�F 17 140
      ��f�T�C�Y 11.31
      
      spl 7501 8 > aaa.bat
      call aaa.bat
      chk-rc 8 75 > bbb.bat
      call bbb.bat
      set-rc 11.31 > ccc.bat
      call ccc.bat
      rec_stk 8 17 140 > ddd.bat
      call ddd.bat

   n. 32bit tiff �摜�̐؂�o��
      rec_crop in/ out/ {x1 y1 x2 y2}
    
      in/: 32bit tiff �č\���摜������f�B���N�g��
      out/: �؂�o����̉摜���o�͂���f�B���N�g��(���s�O�ɍ쐬���邱��)
      x1: �������؂�o���̎n�܂�B
      x2: �������؂�o���̏I���B
      y1: �c�����؂�o���̎n�܂�B
      y2: �c�����؂�o���̏I���B


