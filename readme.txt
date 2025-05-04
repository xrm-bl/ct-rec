中野ソフトをベースにした画像再構成ソフトと規格化辺りの説明

上杉

2025.05.04  ver. 1.4

0. バグやリクエストは作者に連絡してください。

1. 共通の考え
   a. 入力
      基本的にimg形式。

   b. 出力
      CT像は 32bit tiff での出力となる。rec?????.tif (数値は5ケタ)
      それぞれの tiff タグには、
      画素サイズ・回転軸の位置・投影数・回転角オフセット・その画像での最小・最大値
      がその順番に埋め込まれている。
      規格化時には、それまでのタグに加え、規格化時の最小最大値が追記される。
      連続再構成と規格化に関しては、実行完了時に cmd-hst.log にログが残る

   c. プログラムの添え字
      画像再構成ソフトには _t_c などの添え字がついている。
      これらは演算に使用するプロセッサーと再構成フィルターを指定する。
      _P: プロセッサー
          _t: CPU のマルチスレッド機能を使用する。環境変数CBP_THREADSで制御。デフォルトは８スレ。
          _g: GPGPU を使用。添付のexeは CUDA toolkit 10.2 でコンパイル済み。
      _F: フィルター
          _c: Chesler フィルター
          _s: Shepp-Logan フィルター
          _r: Ramachandran(HAN)フィルター
      となっている。

   d. リングアーティファクトの除去
      このバージョンからVo et al.(2018)のArgorism 3型のリング除去機能を設けた。
      CBP計算の直前に実行している。環境変数を指定することでこの機能をON/OFFできる。
      環境変数で KERNEL_SIZE を1に指定するとOFF。それ以外の正の奇数で効果が変わる。
      デフォルト値は5としている(環境変数が定義されていない場合も5になる)。
      また、CPU並列計算も行っており、デフォルト値はOMP_NUM_THREADSを40としている。
      環境変数の設定方法は、sort_filter_omp.c の冒頭に記述があるので参照のこと。

2. 180deg scan。標準的な吸収の画像再構成。

   a. 1枚だけ再構成
      ct_rec_P_F layer {center} {pixel size} {offsetangle}
      
      layer: 再構成するレイヤー(高さ)
      center: 回転軸の位置(pixel)。省略した場合は自動推定する。
      pixel: size: 画素サイズ(um)。省略した場合は1.0になる。
      offset angle: 回転軸の原点オフセット。省略した場合は0.0になる。
      
      *) q????.img があるディレクトリで実行する。
   
   b. 1枚だけ再構成
      tf_rec_P_F layer {center} {pixel size} {offsetangle}
      
      layer: 再構成するレイヤー(高さ)
      center: 回転軸の位置(pixel)。省略した場合は自動推定する。
      pixel: size: 画素サイズ(um)。省略した場合は1.0になる。
      offset angle: 回転軸の原点オフセット。省略した場合は0.0になる。
      
      *) q????.tif があるディレクトリで実行する。
   
   c. 連続再構成
      hp_tg_P_F HiPic/ Dr RC RA0 rec/
      (回転軸が傾いてない場合。全レイヤー)
      
      HiPic/: q????.img が格納されているディレクトリ名。(/ は不要)
      Dr: 画素サイズ (um)
      RC: 回転軸の位置
      RA0: 回転軸の原点オフセット
      rec/: 再構成画像を出力するディレクトリ(計算前に作成すること)
      
      hp_tg_P_F HiPic/ Dr L1 C1 L2 C2 RA0 rec/
      (回転軸が傾いている場合。もしくは一部の領域のみの計算時)
      
      HiPic/: q????.img が格納されているディレクトリ名。(/ は不要)
      Dr: 画素サイズ (um)
      L1: 計算開始レイヤー
      C1: L1での回転軸の位置
      L2: 計算終了レイヤー
      C2: L2での回転軸の位置
      RA0: 回転軸の原点オフセット
      rec/: 再構成画像を出力するディレクトリ(/ は不要。計算前に作成すること)
   
      *) q????.img があるディレクトリの一つ上で実行する。

   d. p画像からの連続再構成
      p_rec_P_F p/ rec/ Dr RC RA0
      (回転軸が傾いてない場合。全レイヤー)
      
      p/: p?????.tif が格納されているディレクトリ名。(/ は不要)
      rec/: 再構成画像を出力するディレクトリ(計算前に作成すること)
      Dr: 画素サイズ (um)
      RC: 回転軸の位置
      RA0: 回転軸の原点オフセット
      
      p_rec_P_F p/ rec/ Dr L1 C1 L2 C2 RA0
      (回転軸が傾いている場合)
      
      p/: p?????.tif が格納されているディレクトリ名。(/ は不要)
      rec/: 再構成画像を出力するディレクトリ(計算前に作成すること)
      Dr: 画素サイズ (um)
      L1: 計算開始レイヤー
      C1: L1での回転軸の位置
      L2: 計算終了レイヤー
      C2: L2での回転軸の位置
      RA0: 回転軸の原点オフセット

3. 360deg scan (offset CT)。標準的な吸収の画像再構成
   a. 回転軸位置の推定
      ofct_xy HiPic/ {Ox1 Ox2 Oy1 Oy2} {MSD.tif}
      
      HiPic/: q????.img が格納されているディレクトリ名(/ は不要)
      Ox1: 横の捜索範囲開始点。(省略可)
      Ox2: 横の捜索範囲終了点。(省略可)
      Oy1: 縦の捜索範囲開始点。(省略可)
      Oy2: 縦の捜索範囲終了点。(省略可)
      MSD.tif: MSDを画像出力ファイル名。(省略可)

      *) 画像出力時にはカッコ内に画像のオフセットが表示される。

   b. 必要メモリ量の確認
      ofct_srec_P_F HiPic/ Rc Oy

      HiPic/: q????.img が格納されているディレクトリ名(/ は不要)
      Rc: 回転軸の位置(左端からの画素数)
      Oy: 縦ずれ量（常に 0 を指定）。

      これで全部を一度に計算するために必要なメインメモリの量が出力される。
      rangeListを決める際の目安になる。

   c. 再構成
      ofct_srec_P_F HiPic/ Rc Oy rangeList Dr RA0 rec/
      
      HiPic/: q????.img が格納されているディレクトリ名(/ は不要)
      Rc: 回転軸の位置(左端からの画素数)
      Oy: 縦ずれ量（常に 0 を指定）。
      rangeList: 再構成するレイヤーを指定。
      Dr: 画素サイズ
      RA0: 回転軸の原点オフセット
      rec/: 再構成画像を出力するディレクトリ(計算前に作成すること)

      *) rangeListの例。
         100レイヤーだけの場合：100
         100-150 レイヤー: 100-150
         全部: -
   d. tiff データから1枚だけ再構成
      otf_rec_P_F layer center {pixel size} {offsetangle}
      
      layer: 再構成するレイヤー(高さ)
      center: 回転軸の位置(pixel)。
      pixel: size: 画素サイズ(um)。省略した場合は1.0になる。
      offset angle: 回転軸の原点オフセット。省略した場合は0.0になる。
      
      *) q????.tif があるディレクトリで実行する。
   


4. 32bit tiff 画像の規格化
    tif_f2i bit rec/ out/ {LACmin LACmax} {x1 y1 x2 y2}
    
    bit: ビット数 0 か 8 か 16 (0 は範囲を調べるだけ)
    rec/: 32bit tiff 再構成画像があるディレクトリ
    out/: 規格化後の画像を出力するディレクトリ(実行前に作成すること)
    LACmin: 規格化時の最小値。(省略可)
    LACmax: 規格化時の最大値。(省略可)
    x1: 横方向切り出しの始まり。(省略可)
    x2: 横方向切り出しの終わり。(省略可)
    y1: 縦方向切り出しの始まり。(省略可)
    y2: 縦方向切り出しの終わり。(省略可)
    LAC の最大最小を省略した場合は、rec/ 中の最大最小値を用いて規格化する。
    規格化時の最大最小値はtiffのタグの最後に付け加えられている。

5. 32bit tiff 画像へのgaussian filter適用
    rec_gf rec/ radius out/
    
    rec/: 32bit tiff 再構成画像があるディレクトリ
    radius: gaussian filter の半値幅
    out/: フィルター後の画像を出力するディレクトリ(実行前に作成すること)

6. その他おまけ
   a. print image description
      pid tiff-file
      tifのタグを表示する
      
   b. シノグラム作成
      sinog layer {skip}
      32bit tiff にて出力。タグには最大と最小値のみ入る。
      skip を指定すると、投影数をskip分の1にして出力する。

   c. b のシノグラムから再構成
      sf_rec_P_F input output {Dr RC RA0}
      32bit tiff にて出力。
      input: float tiff で出来たシノグラム
      output: float tiff で出力
      Dr: 画素サイズ
      RC: 回転軸の位置
      RA0: 回転軸の原点オフセット

   d. img 画像平均
      img_ave file1 file 2... output
      img 画像を平均化する。conv.bat 実行時に必須。

   e. his 分割保存
      spl N-shot N-split
      Zスキャンやエネルギースキャンをしたときに使う可能性が高い。
      his ファイルを分割して img に保存。

   f. rec 連番振りなおし。
      rec_stk num_stack start end
      Zスキャンをしたときに使う。rec の連番を振りなおして1つのスタックにする。

   g. his ファイルを連番の img に変換。
      his2img his-file (x1 x2 y1 y2)
      カッコ内の引数指定にてcropが可能
      
   h. 180deg スキャンから投影像を作成
      ct_prj_f HiPic/ prj/

   i. 32bit tiff からヒストグラム
      tif2hst rec/ (x1 y1 x2 y2)
      rec/: 32bit tiff 再構成画像があるディレクトリ
      x1: 横方向切り出しの始まり。(省略可)
      y1: 縦方向切り出しの始まり。(省略可)
      x2: 横方向切り出しの終わり。(省略可)
      y2: 縦方向切り出しの終わり。(省略可)

   j. 8bit or 16bit CT像の直交回転
      si_rar orgDir nameFile hAxis vAxis dAxis newDir
      si_rar orgDir nameFile hAxis vAxis dAxis sliceNo newTIFF

      例: si_rar.exe ro_xy - +y +z +x ro_yz
      例: si_rar.exe ro_xy - +z +x +y ro_zx

   k. 8bit or 16bit CT像のビニング
      si_sir orgDir nameFile Bxyz newDir
      si_sir orgDir nameFile Bx By Bz newDir

      例: si_sir ro_xy - 2 ro_2x2x2
      例: si_sir ro_xy - 2 2 1 ro_2x2x1

   l. 8bit or 16bit CT像にgaussian filter
      si_gf orgDir nameFile radius {bias} newDir

      例: si_gf rh - 1.0 rh_gf1
      例: si_gf rh - 1.0 0.0 rh_gf1

   m. Zスキャンの自動スタック(例)
      1回の撮影枚数: 7501
      Zスタックの数：8
      オーバーラップの場所： 17 140
      画素サイズ 11.31
      
      spl 7501 8 > aaa.bat
      call aaa.bat
      chk-rc 8 75 > bbb.bat
      call bbb.bat
      set-rc 11.31 > ccc.bat
      call ccc.bat
      rec_stk 8 17 140 > ddd.bat
      call ddd.bat

   n. 32bit tiff 画像の切り出し
      rec_crop in/ out/ {x1 y1 x2 y2}
    
      in/: 32bit tiff 再構成画像があるディレクトリ
      out/: 切り出し後の画像を出力するディレクトリ(実行前に作成すること)
      x1: 横方向切り出しの始まり。
      x2: 横方向切り出しの終わり。
      y1: 縦方向切り出しの始まり。
      y2: 縦方向切り出しの終わり。


