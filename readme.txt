中野ソフトをベースにした画像再構成ソフトと規格化辺りの説明

上杉

2026.08.06  ver. 2.4
2026.07.30  ver. 2.3
2026.07.02  ver. 2.2
2026.06.30  ver. 2.1
2026.06.30  ver. 2.0
2026.05.04  ver. 1.7

【ver 2.4 の変更点】
  ・透過率不足への画素単位ガードを導入。試料が厚い/高密度で透過信号が dark
    レベル以下になると log(I0/(I-dark)) が +Inf/NaN を返し、GPU リング除去が
      CUDA error sort_filter_g.cu:266: an illegal memory access was encountered
    で停止(CPU 版は当該列のソートが黙って破損)していたのを発生源で修正。
  ・環境変数 CT_REC_BLACK_THRESH の意味を「ラインの平均信号による黒投影判定」
    から「画素ごとの信号下限(床)」に変更(既定 1 → 2)。単位(dark を超えた
    カウント数)は不変。これ以下の画素は床止めされ、投影値は
    log(I0/CT_REC_BLACK_THRESH) で頭打ちになるため Inf/NaN は構造的に発生
    しない。床より上の画素は一切変更されない。
  ・黒投影(missing angle)判定は新設の CT_REC_BLACK_FRAC(床止め画素の割合、
    既定 0.5)で行う。板状試料では従来と同じ透過率で発火する。視野のほぼ全面を
    覆う超高密度試料で投影が破棄されすぎる場合は 0.8 程度に上げる。詳細は 1e。
  ・両変数とも起動時に有効値を stderr へ表示し、不正値は既定に戻す。実行終了時に
    床止め画素数・最小信号・最大投影値の要約を出力する。
  ・併せて修正: 空気基準参照の off-by-one(未初期化領域を読んでいたため再構成値が
    わずかに変わる)、全投影が黒判定になったときの 0 除算、回転中心計算の
    無防備な log。技術ノート 20260806_low_transmission_guard.md も参照。

【ver 2.3 の変更点】
  ・打ち切り(カッピング)補正を、全再構成ソフト共通の CBP 層に追加。環境変数
    PAD_THRESH で制御し、未設定時は 0.3(自動判定ON)。PAD_THRESH=0 で強制OFF
    (以前の版と同一の結果)。詳細は 1f。
  ・CPU 版の既定スレッド数を固定8から「実行PCの論理コア数-1」に変更
    (環境変数 CBP_THREADS)。GPU 版の添付 exe は CUDA toolkit 13.2 でビルド
    (CUDA 13 系の最低要件は Turing 世代/sm_75 以降)。
  ・ofct_DO にビュー間引き(stride)、ofct_DO_g にビュー対平均(navg)と
    間引き(stride)を追加。回転軸推定を短時間で回せる。
  ・hp_tg / ofct_srec が投影を並列に読み込むようにした
    (環境変数 HPTG_READ_THREADS、既定16)。
  ・hp_tg 相当の全層再構成を ct_rec の並列実行で行う ct_rec_loop を追加(2d)。
    Windows 用 .bat と Linux 用 .sh の両方があり、dark.tif にも対応。
  ・tif で連続取得したデータの分割保存 act_spl2 / act_spl を追記(6q)。
  ・再構成済み CT 像の再投影・再再構成 rec2rec を追加(6r)。
  ・median + gaussian の tif_mgf に GPU 版 tif_mgf_g を追加(6p)。
  ・同梱の libtiff を 3.6.0 から 4.6.0 に更新。8/16/32bit tiff の入出力は
    従来と完全に互換。
  ・環境変数の一覧(既定値つき)を 20260723_CT_env_vars.md にまとめた。

【ver 2.2 の変更点】
  ・オフセットCT 回転軸推定 ofct_DO の GPU 版 ofct_DO_g（ofct_DO.cu）を追加。
    ホスト側処理は ofct_DO と同一で、推定される中心/Oy は CPU 版と一致する。
    ビュー対ごとの MSD 計算を GPU で実行。CPU 版 ofct_DO は残す。

【ver 2.1 の変更点】
  ・ct_rec と tf_rec を統合。ct_rec が dark ファイル（dark.img / dark.tif）を見て
    .img / .tif を自動判別するようになった。内部で ct_rec.c と tf_rec.c を統合した
    ct_rec_c.c を使用。これに伴い TIFF 専用版 tf_rec を廃止:
      tf_rec_P_F → ct_rec_P_F
  ・オフセットCT 1枚再構成 otf_rec も同様に img/tif 自動判別化し改名（TIFF 専用版 otf_rec は廃止）:
      otf_rec_P_F → otct_rec_P_F
  ・投影像生成 ct_prj_f も img/tif 自動判別化（ct_prj_f.c + tf_prj_f.c → ct_prj_f_c.c）。TIFF 専用版 tf_prj_f は廃止。
  ・オフセットCT 回転軸推定を ofct_xy → ofct_DO に置き換え（img/tif 自動判別、ofct_srec 用コマンドを提案表示）。
  ・オフセットCT 再構成は使用可能メモリを自動判定して自動 chunk するため、手動のメモリ確認手順（旧 3b）を削除。

【ver 2.0 の変更点】
  ・hp_tg と tf_tg を統合。hp_tg が入力ディレクトリの dark ファイル
    （dark.img / dark.tif）を見て .img / .tif を自動判別するようになった。
    内部リーダ rhp.c(.img) と rtf.c(.tif) を統合した rhp_c.c を使用。
  ・これに伴い TIFF 専用版を廃止し、以下に一本化:
      tf_tg_P_F     → hp_tg_P_F
      oftf_srec_P_F → ofct_srec_P_F
      oftf_xy       → ofct_xy

0. バグやリクエストは作者に連絡してください。

1. 共通の考え
   a. 入力
      基本的にimg形式。1枚再構成(ct_rec)・連続再構成(hp_tg)・オフセットCT(ofct_srec/ofct_DO)は、
      dark.img があれば img、無く dark.tif があれば tiff を自動判別して読み込む。

   b. 出力
      CT像は 32bit tiff での出力となる。rec?????.tif (数値は5ケタ)
      それぞれの tiff タグには、
      画素サイズ・回転軸の位置・投影数・回転角オフセット・その画像での最小・最大値
      がその順番に埋め込まれている。
      規格化時には、それまでのタグに加え、規格化時の最小最大値が追記される。
      連続再構成と規格化に関しては、実行完了時に cmd-hst.log にログが残る。
      回転軸推定(ofct_DO)、変換・平均系のユーティリティ、tif_* のフィルタ群も
      同様に cmd-hst.log に記録する。tif_* のフィルタはコマンドとパラメータを
      1行のタブ区切りで残すので、表計算ソフトで開いても崩れない。

   c. プログラムの添え字
      画像再構成ソフトには _t_c などの添え字がついている。
      これらは演算に使用するプロセッサーと再構成フィルターを指定する。
      _P: プロセッサー
          _t: CPU のマルチスレッド機能を使用する。環境変数CBP_THREADSで制御。
              デフォルトは実行PCの論理コア数-1(下限1)。以前は固定8だった。
          _g: GPGPU を使用。添付のexeは CUDA toolkit 13.2 でコンパイル済み
              (CUDA 13 系の最低要件は Turing 世代/sm_75 以降)。
      _F: フィルター
          _c: Chesler フィルター
          _s: Shepp-Logan フィルター
          _r: Ramachandran(HAN)フィルター
      となっている。
      ただし例外があり、シノグラム再構成 sf_rec は CPU 版のみ(sf_rec_t_F)、
      再投影 rec2rec は CPU 版に _t が付かない(rec2rec_F / rec2rec_g_F)。

   d. リングアーティファクトの除去
      バージョン1.4からVo et al.(2018)のAlgorithm 3型のリング除去機能を設けた。
      CBP計算の直前に実行している。環境変数を指定することでこの機能をON/OFFできる。
      環境変数で KERNEL_SIZE を1に指定するとOFF。それ以外の正の奇数で効果が変わる。
      デフォルト値は5としている(環境変数が定義されていない場合も5になる)。
      また、リング除去処理はOpenMPによるCPU並列計算で行っており、デフォルト値は
      OMP_NUM_THREADSを40としている。これは1cで述べたCBP_THREADS（逆投影計算用の
      スレッド数）とは独立した設定である。

   e. 透過率不足と missing angle の処理
      試料が厚い/高密度の場合、透過信号が dark レベルに達するか下回る事がある。
      すると I-dark が 0 または負になり log(I0/(I-dark)) が +Inf / NaN を返す。
      これは GPU リング除去を

        CUDA error sort_filter_g.cu:266: an illegal memory access ...
        ring removal image processing failed

      で停止させ、CPU 版では当該列のソートを黙って壊していた(2026-08 修正)。

      環境変数 CT_REC_BLACK_THRESH は「画素ごと」の信号下限で、単位は dark を
      超えたカウント数(未設定の場合は 2 とする)。これ以下の画素は不透明とみなし
      この値で床止めするので、投影値は log(I0/CT_REC_BLACK_THRESH) で頭打ちになり
      無限大にはならない。床より上の画素は一切変更しない。0 以下の値は不正として
      拒否する(この保護は無効化できない)。有効値は起動時に、床止めした画素数・
      最小信号・最大投影値は終了時に表示される。

      板状試料などの場合、角度によっては透過率が極端に低下する(missing angle)。
      床止めされた画素の割合が CT_REC_BLACK_FRAC(未設定の場合は 0.5)を超えた投影は
      黒とみなし、良い投影の平均で置き換える。従来どおり CT_REC_BLACK_THRESH を
      1, 10, 100, 1000 など変更すればこの起こりやすさを調整できる。一様に減衰した
      ラインでは発火する透過率は従来と同じである。

      ただし、このしきい値は破棄されずに残る投影の投影値も頭打ちにするので、
      大きな値は高密度試料の内部を平坦に潰す。板状試料を扱う場合を除き既定値を推奨する。

   f. 打ち切り(カッピング)補正
      試料が視野からはみ出したデータ(打ち切りデータ)では、投影の両端にできる
      段差とランプフィルタの裾との畳み込みが、再構成像の外周に明るいリム
      (カッピング)を作る。これを抑えるため、フィルタをかける前に投影の両端を
      端の値のままコサインで減衰させながら幅 N/2 だけ外挿する機能を、全再構成
      ソフト共通の CBP 層に設けた。逆投影の計算量は変わらない(FFT長のみ2倍)。
      視野外にはみ出した試料による打ち切りアーティファクトを投影データの外挿で
      抑える手法は Ohnesorge et al.(2000) 型のものである
      (Med. Phys. 27(1), 39-46)。本実装は端値ホールド + コサイン減衰による
      簡易版で、論文の手法をそのまま再現したものではない。
      適用するかどうかは自動判定で、シノグラム両端列の平均振幅が全体平均の
      一定比率を超えたときだけ有効になる。この比率を環境変数 PAD_THRESH で
      指定し、未設定の場合は 0.3(自動判定ON)とする。有効になったときは
        PAD_THRESH: truncation pad enabled (W=...)
      と表示する。視野内に収まっている試料は端の値がほぼ0のため、自動的に
      無適用となる。PAD_THRESH=0(以下)を指定すると強制的にOFFで、以前の版と
      完全に同一の結果になる。
      ct_rec / hp_tg / p_rec / ofct_rec / ofct_srec / sf_rec / rec2rec の
      すべてに効く。ただし rec2rec の入力は再構成円の外で端が減衰しているので、
      この場合は 0.1-0.2 程度の小さい値が適当。

   g. その他の環境変数
      上記のほかに、投影の並列読み込みスレッド数(HPTG_READ_THREADS、既定16。
      hp_tg / ofct_srec で有効)、メモリに合わせたチャンク分割の指定
      (HPTG_MEM_FRACTION / HPTG_MEM_LIMIT_MB / HPTG_CHUNK_ROWS)、使用する
      GPU 番号(CUDA_GPU、既定0)、入力ファイル名の上書き(RHP_O / RHP_D / RHP_Q)
      などがある。既定値を含む一覧は 20260723_CT_env_vars.md を参照のこと。
      同一ノードで複数プロセスを並走させる場合は HPTG_MEM_FRACTION を
      下げること(既定のままだと全プロセスが空きメモリの9割を要求する)。

2. 180deg scan。標準的な吸収の画像再構成。

   a. 1枚だけ再構成
      ct_rec_P_F layer {center} {pixel size} {offsetangle}
      
      layer: 再構成するレイヤー(高さ)
      center: 回転軸の位置(pixel)。省略した場合は自動推定する。
      pixel size: 画素サイズ(um)。省略した場合は1.0になる。
      offset angle: 回転軸の原点オフセット。省略した場合は0.0になる。
      
      *) q????.img もしくは q????.tif があるディレクトリで実行する。(dark.img / dark.tif で自動判別)
   
   b. 連続再構成
      hp_tg_P_F HiPic Dr RC RA0 rec
     (回転軸が傾いてない場合。全レイヤー)
      
      HiPic: q????.img もしくは q????.tif が格納されているディレクトリ名。(/ は不要)
      Dr: 画素サイズ (um)
      RC: 回転軸の位置
      RA0: 回転軸の原点オフセット
      rec: 再構成画像を出力するディレクトリ(計算前に作成すること)
      
      hp_tg_P_F HiPic Dr L1 C1 L2 C2 RA0 rec
      (回転軸が傾いている場合。もしくは一部の領域のみの計算時)
      
      HiPic: q????.img もしくは q????.tif が格納されているディレクトリ名。(/ は不要)
      Dr: 画素サイズ (um)
      L1: 計算開始レイヤー
      C1: L1での回転軸の位置
      L2: 計算終了レイヤー
      C2: L2での回転軸の位置
      RA0: 回転軸の原点オフセット
      rec: 再構成画像を出力するディレクトリ(/ は不要。計算前に作成すること)
   
      *) q????.img があるディレクトリの一つ上で実行する。

   c. p画像からの連続再構成
      p_rec_P_F p rec Dr RC RA0
      (回転軸が傾いてない場合。全レイヤー)
      
      p: p?????.tif が格納されているディレクトリ名。(/ は不要)
      rec: 再構成画像を出力するディレクトリ(計算前に作成すること)
      Dr: 画素サイズ (um)
      RC: 回転軸の位置
      RA0: 回転軸の原点オフセット
      
      p_rec_P_F p rec Dr L1 C1 L2 C2 RA0
      (回転軸が傾いている場合)
      
      p: p?????.tif が格納されているディレクトリ名。(/ は不要)
      rec: 再構成画像を出力するディレクトリ(計算前に作成すること)
      Dr: 画素サイズ (um)
      L1: 計算開始レイヤー
      C1: L1での回転軸の位置
      L2: 計算終了レイヤー
      C2: L2での回転軸の位置
      RA0: 回転軸の原点オフセット

   d. ct_rec の並列実行による連続再構成
      ct_rec_loop.bat HiPic Dr RC RA0 rec {Njobs}            (Windows)
      ct_rec_loop.bat HiPic Dr L1 C1 L2 C2 RA0 rec {Njobs}
      ct_rec_loop.sh  HiPic Dr RC RA0 rec {Njobs}            (Linux。bin/ にある)
      ct_rec_loop.sh  HiPic Dr L1 C1 L2 C2 RA0 rec {Njobs}

      引数は hp_tg と同じで、最後に並列ジョブ数 Njobs を足すだけ(省略時は1)。
      hp_tg と同じ全層再構成を、1枚再構成の ct_rec をレイヤーごとに並列実行して
      行うラッパー。hp_tg が使えない場合や、GPU をもっと使いたい場合に用いる。
      HiPic の dark.img / dark.tif を見て img/tif を自動判別する。
      RC 指定(中心固定)の形では dark ファイルの高さから全レイヤーを自動で決める。
      結果は rec/ に集め、レイヤーごとのログは rec/log/ に置く。
      使用する実行ファイルは既定で ct_rec_g_c。.sh 版は環境変数 CT_EXE で変更可能。
      Njobs は GPU メモリに合わせて調整すること。

      *) hp_tg と同じく、HiPic ディレクトリの一つ上で実行する。

3. 360deg scan (offset CT)。標準的な吸収の画像再構成
   a. 回転軸位置の推定
       ofct_DO   raw {stride}          (CPU)
       ofct_DO_g raw {navg {stride}}   (GPU; 結果は ofct_DO と同一)

       raw: q????.img もしくは q????.tif が格納されているディレクトリ名(/ は不要)
            dark.img / dark.tif の有無で img/tif を自動判別。
       stride: ビュー対を stride 個に1つだけ使って間引く。省略した場合は 1(全部使用)。
               この処理は投影の読み込み律速なので、実行時間はほぼ 1/stride になる。
               推定される中心/Oy はほとんど変わらない。
       navg: GPU 版のみの引数。-log を取る前に平均するビュー対の数。省略した場合は 10。
             平均によって軸探索の S/N が上がる。

       オフセットCTデータから回転軸位置を推定し、そのまま実行できる
       ofct_srec コマンド（中心と Oy）を提案表示する。

       ofct_DO_g はビュー対ごとの MSD を GPU で計算する（ホスト側処理は同一で結果も一致）。

       いずれも MSD の計算前にガウシアン平滑化をかけている。環境変数
       OFCT_DO_SMOOTH で σ を指定し、未設定の場合は 1.0。0 を指定すると平滑化しない。

   b. 再構成
      ofct_srec_P_F HiPic Rc Oy rangeList Dr RA0 rec
      
      HiPic: q????.img もしくは q????.tif が格納されているディレクトリ名(/ は不要)
      Rc: 回転軸の位置(左端からの画素数)
      Oy: 縦ずれ量（常に 0 を指定）。
      rangeList: 再構成するレイヤーを指定。
      Dr: 画素サイズ
      RA0: 回転軸の原点オフセット
      rec: 再構成画像を出力するディレクトリ(計算前に作成すること)

      *) rangeListの例。
         100レイヤーだけの場合：100
         100-150 レイヤー: 100-150
         全部: -

   c. img/tif データから1枚だけ再構成（オフセットCT）
      ofct_rec_P_F layer center {pixel size} {offsetangle}
      
      layer: 再構成するレイヤー(高さ)
      center: 回転軸の位置(pixel)。
      pixel size: 画素サイズ(um)。省略した場合は1.0になる。
      offset angle: 回転軸の原点オフセット。省略した場合は0.0になる。
      
      *) q????.img もしくは q????.tif があるディレクトリで実行する。(dark.img / dark.tif で自動判別)
   


4. 32bit tiff 画像の規格化
    tif_f2i bit rec out {LACmin LACmax} {x1 y1 x2 y2}
    
    bit: ビット数 0 か 8 か 16 (0 は範囲を調べるだけ)
    rec: 32bit tiff 再構成画像があるディレクトリ
    out: 規格化後の画像を出力するディレクトリ(実行前に作成すること)
    LACmin: 規格化時の最小値。(省略可)
    LACmax: 規格化時の最大値。(省略可)
    x1: 横方向切り出しの始まり。(省略可)
    y1: 縦方向切り出しの始まり。(省略可)
    x2: 横方向切り出しの終わり。(省略可)
    y2: 縦方向切り出しの終わり。(省略可)
    LAC の最大最小を省略した場合は、rec/ 中の最大最小値を用いて規格化する。
    規格化時の最大最小値はtiffのタグの最後に付け加えられている。

5. 32bit tiff 画像へのgaussian filter適用
    rec_gf rec radius out
    
    rec: 32bit tiff 再構成画像があるディレクトリ
    radius: gaussian filter の半値幅
    out: フィルター後の画像を出力するディレクトリ(実行前に作成すること)

6. その他おまけ
   a. tif 画像に埋め込まれたタグ(画素サイズ、投影数など)を表示する
      pid tiff-file
      
   b. シノグラム作成
      sinog layer {skip}
      32bit tiff にて出力。タグには最大と最小値のみ入る。
      skip を指定すると、投影数をskip分の1にして出力する。

   c. b のシノグラムから再構成
      sf_rec_t_F input output {Dr RC RA0}
      32bit tiff にて出力。CPU 版のみで GPU 版はない。
      input: float tiff で出来たシノグラム
      output: float tiff で出力
      Dr: 画素サイズ
      RC: 回転軸の位置
      RA0: 回転軸の原点オフセット

   d. img 画像平均
      img_ave file1 file2... output
      img 画像を平均化する。conv.bat 実行時に必須。

   e. his 分割保存
      spl N-shot N-split
      Zスキャンやエネルギースキャンをしたときに使う可能性が高い。
      his ファイルを分割して img に保存。
      tif で連続取得した場合はフィルタ処理も兼ねる act_spl2 を使う(6q)。

   f. rec 連番振りなおし。
      rec_stk num_stack start end
      Zスキャンをしたときに使う。rec の連番を振りなおして1つのスタックにする。

   g. his ファイルを連番の img に変換。
      his2img his-file (x1 x2 y1 y2)
      カッコ内の引数指定にてcropが可能
      
   h. 180deg スキャンから投影像を作成
      ct_prj_f HiPic prj
      入力は dark.img があれば img、無く dark.tif があれば tiff を自動判別。

   i. 32bit tiff からヒストグラム
      tif2hst rec (x1 y1 x2 y2)
      rec: 32bit tiff 再構成画像があるディレクトリ
      x1: 横方向切り出しの始まり。(省略可)
      y1: 縦方向切り出しの始まり。(省略可)
      x2: 横方向切り出しの終わり。(省略可)
      y2: 縦方向切り出しの終わり。(省略可)

   j. 8bit or 16bit CT像の直交回転
      si_rar orgDir nameFile hAxis vAxis dAxis newDir
      si_rar orgDir nameFile hAxis vAxis dAxis sliceNo newTIFF

      例: si_rar.exe ro_xy - +y +z +x ro_yz
      例: si_rar.exe ro_xy - +z +x +y ro_zx

   k. CT像のビニング
      8bit or 16bit の場合:
      si_sir orgDir nameFile Bxyz newDir
      si_sir orgDir nameFile Bx By Bz newDir

      例: si_sir ro_xy - 2 ro_2x2x2
      例: si_sir ro_xy - 2 2 1 ro_2x2x1

      32bit tiff (rec?????.tif) の場合:
      rec_sir in out B

      in: 32bit tiff 再構成画像があるディレクトリ
      out: ビニング後の画像を出力するディレクトリ(実行前に作成すること)
      B: ビニング数。x, y と z(スライス方向)すべてに同じ値を適用する。

      B×B×B ボクセルの平均を 32bit float のまま出力する(整数への丸めはしない)。
      出力の枚数と縦横は 1/B になる。ファイル名は先頭が入力の先頭番号と一致し、
      以降 1 ずつ増える。例えば rec00003 から rec00011 を B=3 で処理すると
      rec00003, rec00004, rec00005 の 3 枚になる。
      ビニング数で割り切れない端の列・行・スライスは切り捨てる(si_sir は端も
      残して実画素数で平均するので、この点だけ挙動が異なる)。
      tiff タグは新しい格子に合わせて更新する。画素サイズは B 倍、回転軸の
      位置は (RC-(B-1)/2)/B、投影数と回転角オフセットはそのまま、最小最大値は
      結果から再計算する。
      出力の番号が入力の番号と重なるため、入力と同じディレクトリを出力に
      指定した場合はエラーで停止する。

      例: rec_sir rec rec_b2 2

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
      rec_crop in out {x1 y1 x2 y2}
    
      in: 32bit tiff 再構成画像があるディレクトリ
      out: 切り出し後の画像を出力するディレクトリ(実行前に作成すること)
      x1: 横方向切り出しの始まり。
      y1: 縦方向切り出しの始まり。
      x2: 横方向切り出しの終わり。
      y2: 縦方向切り出しの終わり。

   o. 任意bit数のCT像に filterかける
      Bilateral filter, Gaussian Filter, Median filter などがある。
      各ソフト _g で終わるものはGPU使用版(Cuda toolkit のバージョンを確認すること)。
      tif_blf, tif_blf_g, tif_gsf, tif_gsf_g, tif_mdf, tif_mdf_g
      
      Bilateral filter
      tif_blf[_g] <input_dir> <output_dir>  [kernel_size] [spatial_sigma] [intensity_sigma]
      
      Gaussian filter
      tif_gsf[_g] <input_dir> <output_dir> [sigma]
      
      Median filter
      tif_md[_g]f <input_dir> <output_dir> [kernel_size]
      
      []の変数は省略すると、自動的にパラメータを決める。

      例: tif_blf_g rh rh_blf1
      例: tif_blf rh rh_blf1 5 20 200
      
      その他のフィルタ等の詳細は20260504_filter_readme.pdfを参照のこと。
      

   p. tif画像にmedian filterかけてからgaussian filterをかける。
      tif_mgf[_g] <input_file> <output_file> [median_kernel_size] [gaussian_sigma]
      X線透過像に散乱光のノイズなどがある場合に、除去するようなときに用いる。
      _g は GPU 使用版。

      例: tif_mgf a000401.tif 001/raw/a0401.tif 3 1
      例: tif_mgf a000401.tif 001/raw/a0401.tif 0 1
      例: tif_mgf_g a000401.tif 001/raw/a0401.tif 3 1

   q. tif で連続取得したデータの分割保存(フィルタ処理つき)
      act_spl2 N-shot N-split {m_kernel_size} {g_kernel_size}
      act_spl  N-shot N-split {kernel_size}

      Zスキャンやエネルギースキャンで a??????.tif を連続取得した場合に、
      1回の撮影分ずつ 001/raw, 002/raw, ... へ振り分けて保存する。
      6e の spl (his を img に分割) の tif 版にあたる。
      N-shot: 1回の撮影枚数
      N-split: 分割数(Zスタックの数)
      act_spl2 は振り分けと同時に tif_mgf (median + gaussian) をかける。
        m_kernel_size: median のカーネルサイズ(省略時は 0 = かけない)
        g_kernel_size: gaussian の σ(省略時は 0 = かけない)
      act_spl は gf_sd による gaussian のみで、kernel_size を省略した場合は
      フィルタなしの移動だけを行う。
      どちらも実行ディレクトリの conv.bat と output.log を各 raw/ にコピーし、
      conv.bat 中の img を tif に書き換える。実行後 cmd-hst.log に記録が残る。
      フィルタは Windows では start /b、Linux では fork+exec で背景実行する。
      Linux での同時実行数は環境変数 ACT_SPL_JOBS(既定8)で制限し、終了前に
      全ジョブの完了を待ち合わせる。

      例: act_spl2 7501 8 3 1
      例: act_spl  7501 8

   r. 再構成済み CT 像の再投影・再再構成
      rec2rec_F   rec_in rec_out {Nt}      (CPU)
      rec2rec_g_F rec_in rec_out {Nt}      (GPU)

      rec_in: 32bit float の CT 断面(rec*.tif、正方形)があるディレクトリ
      rec_out: 出力ディレクトリ(入力と同じファイル名で出力)
      Nt: 順投影のビュー数。省略した場合は入力 tiff のタグにある投影数を使い、
          読めない場合は画像サイズ N とする。0 を指定しても自動になる。

      入力像を順投影(Radon 変換)してシノグラムを作り、リング除去 + CBP で
      もう一度再構成する。値のスケール変換はしないので、入力の値がそのまま
      出てくる。1f の打ち切り補正を併用する場合は、PAD_THRESH を 0.1-0.2
      程度に下げるとよい。

      例: rec2rec_g_c rec rec2 1800

      

