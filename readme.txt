中野ソフトをベースにした画像再構成ソフトと規格化辺りの説明

上杉

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

【ver 2.8 の変更点】
  ・his_rec_P_F、his_srec_P_F、his_DO、his_sg(2g)を新設。a.his から
    直接、再構成・回転軸推定・シノグラム作成を行う。ビームラインの conv.bat(フレームを
    dark.img / q????.img にする手順: img_ave / ren / copy)をメモリ上で解釈し、
    output.log は従来どおり読むので、his2img と conv.bat の実行は不要で .img
    ファイルは作られない。output.log の角度欄から 180° スキャンと 360° オフ
    セットスキャン(試料像が 350° に達する。ofct_rec が従来から用いている判定)
    を自動判別し、対応する本体を呼ぶ: his_rec = ct_rec / ofct_rec(1 枚、a.his
    と同じ階層で実行。center 省略時はそのレイヤーから回転軸を推定、オフセット
    では ofct_DO の探索の 1 レイヤー版)、his_srec = hp_tg / ofct_srec(全スライス、1 つ上で実行。
    引数はどちらも hp_tg 形式)、his_DO = hp2DO / ofct_DO(回転軸推定)、
    his_sg = sinog / of_sinog(1 レイヤーのシノグラム、a.his と同じ階層で実行)。
    Z スキャン(1 つの a.his に複数スキャン)は環境変数 RHIS_SCAN=k でスキャンを選び、
    rc-check-his / gen-all-his(exe/*.bat、bin/*.sh)が spl と conv.bat なしに全スキャン
    の軸チェックと全体再構成を行う(スキャン数と中央レイヤーは his_info が報告)。
    his_tb_rec_P_F / his_tb_srec_P_F(または his_* 全般で RHIS_TB=1)は、入射光画像の
    代わりに全試料像の整数平均を使う。his_srec_g_F_n2i は Noise2Inverse 用に
    ビューを K = 2 または 4 組に交互分割して rec/0～rec/K-1 に再構成する(GPU のみ)。
    平均は
    img_ave と同じ整数演算(切り捨て)なので、結果は従来の手順とビット一致する
    (CPU 版・GPU 版で検証済み。his_srec は TIFF の ARTIST タグのみ異なる)。
    16bit の his のみ対応。暗電流・入射光の枚数や入射光のタイミングは
    conv.bat と output.log から読み取り、固定値は持たない。既存プログラムの
    動作は不変(2 本を 1 つの his_* に同居させるため、ファイル内の名前を
    static にしただけ)。
  ・ofct_srec / ofct_rec / of_sinog: オフセットスキャンの回転軸位置に 180° スキャン
    と同様に小数を指定できるようにした。対向像を小数部の 2 倍だけ線形補間でずらして
    重ね、回転中心を小数部だけ動かす。整数指定の結果は従来と同一。ofct_rec は
    従来も小数を受け付けたが貼り合わせは切り捨てた整数で行っており、これを修正。
    his_rec の単層の回転軸推定も小数のまま使う(his_sg も同じ推定を使う)。
  ・ofct_srec: 傾いた回転軸に対応(L1 C1 L2 C2 形式、3b)。軸位置をスライスごとに
    補間して使い、全スライスを最大幅で出力(狭いスライスは中央配置、余白は行の
    端の値)、ImageDescription の 2 列目にそのスライスの Rc+0.5 を記録する
    (ofct_rec と同じ規約)。his_srec の L1 C1 L2 C2 形式はオフセット時にこれを使う。
    C1 = C2 なら出力は従来と同一。
  ・CPU 版(cbp_thread_int.c): 逆投影の検出器座標の固定小数点が 32bit で、再構成幅が
    5792 画素を超えると桁あふれしていた。これに達するのはオフセット系(ofct_srec /
    ofct_rec / his_* の 360° スキャン)で 4096 画素の検出器を使った場合だけ(例 Rc 3809
    → 幅 7618)で、CPU 版の画像の幾何が崩れていた(GPU 版は正常)。座標を 64bit にし、
    それ以下の幅では従来とビット一致。GPU 版はもともと影響なし。

【ver 2.7 の変更点】
  ・act_tg_g_F(2f)と ofact_srec_g_F(3d)を新設。16bit + 8bit を直接出力する
    連続再構成で、再構成までは hp_tg / ofct_srec と同一。32bit の rec*.tif の
    代わりに、各スライスを rec16/rh%05d.tif と rec8/ro%05d.tif として同時に
    書き出し、規格化範囲は引数で 2 組与える(rec16 LACmin16 LACmax16 rec8
    LACmin8 LACmax8)。量子化・偶数サイズへの切り詰め・8 欄の ImageDescription
    は tif_f2i と同じ規則で、「hp_tg / ofct_srec → tif_f2i」で作ったファイルと
    同一(CPU 版・GPU 版で検証済み)。自動CT装置用に作成したが、通常の CT 計測
    にも使用できる。ただし 16bit と 8bit の規格化定数は実行時に必ず指定する
    (自動設定はない)。ビルドは GPU 版のみ(_r / _s / _c)。
  ・act_rec_g_r / ofact_rec_g_r(ver 2.6 で新設した 1 枚再構成、8/16bit 直接
    出力)を廃止。ビルド・配布の対象外。act_tg / ofact_srec、または ct_rec /
    ofct_rec と tif_f2i を用いる。
  ・hp_tg: パイプライン版(ver 2.6 更新分)で、傾き軸形式(L1 C1 L2 C2)のとき
    ImageDescription の回転中心欄が次スライスの値になることがあった(画素値は
    無関係)。修正。
  ・Windows ビルド: 日本語コメントを含むソースに UTF-8 BOM を付け、cl に /utf-8
    を指定するようにして、コードページ関連の警告 C4819 をなくした。Linux と
    生成コードには影響なし。
  ・ct_rec / ofct_rec(2026-10-13 更新): スライスの全投影が黒と判定された
    場合(X線が全く入っていない。1e 参照)、何も出力しない代わりに、通常と同じ
    名前・サイズのプレースホルダ画像を書き出す。ファイル名、「fail!」、
    「all projections black」を大きな文字で描いたもの(背景 0、文字 1.0)で、
    ImageDescription の末尾に「FAILED: all projections black」が付く。終了
    コードは従来どおり非ゼロ(1)なので、バッチからの失敗検知は変わらない。
    ofct_rec はこの場合に全ゼロの画像を終了コード 0 で出力していたが、ct_rec
    と同じ挙動にした。良い投影が 1 枚でもあるスライスの出力は不変。

【ver 2.6 の変更点】
  ・自動CT装置用の1枚再構成 act_rec_g_r / ofact_rec_g_r(8/16bit 直接出力)を
    新設(ver 2.7 で廃止)。
  ・GPU 版全ツール: GPUメモリ不足時のエラーメッセージを改良。空き/総容量と
    対処ヒントを表示し、cuFFT はエラーコード名も表示する。再構成(cbp 層)は
    確保前に必要量を見積もり、不足見込みなら警告する(Windows のドライバは
    VRAM 不足時にエラーを出さずシステムRAMへ退避して大幅に遅くなるだけの
    ため、この警告が唯一の手掛かりになる)。
  ・act_spl2 / act_spl(Linux): 背景フィルタジョブの同時実行数上限(環境変数
    ACT_SPL_JOBS)の既定値を 8 から 16 に変更。上限の動作自体は検証済み
    (既定 16、ACT_SPL_JOBS=3 で 3)。効いていないように見える場合は、バイナリ
    が 2026-07-23 より古いか、tif_mgf / gf_sd が PATH にないことが原因。
    Windows は従来どおり start /b で無制限。
  ・GPU 版 ct_rec のメモリ使用量と並列実行数の見積もりメモ
    20260919_ct_rec_gpu_memory.md を追加。
  ・hp_tg / ofct_srec: スライスごとのオーバーヘッドを削減(結果は不変)。従来は
    スライスごとに約 1.4 GB の作業バッファ(リング除去の入出力、出力画像)を
    malloc/free しており、1 回の実行で約 1.85 億回のページフォルトが発生して
    sys 時間の大半を占め、さらにシノグラムと再構成像の 1 スレッドコピーが
    5 本あった。現在はリング除去を CBP の投影バッファ上でインプレース実行
    (GPU 版は pinned なので転送もフル速度)、hp_tg は再構成バッファから直接
    TIFF を書き出し、ofct_srec は float 変換バッファを保持、残るコピー(投影
    の取り込み、黒投影判定、再構成像の受け渡し)は OpenMP 並列にした。
    sort_filter_g.cu はデバイスバッファを呼び出し間で保持する(従来は毎回
    cudaMalloc/cudaFree × 4。解放用 sort_filter_gpu_release() を追加。
    result_data と image_data は同一バッファ可)。各演算の順序は従来と同一
    なので出力はバイト単位で一致する(黒投影を含む合成データで CPU 版を
    旧版と比較して検証済み)。
    さらに hp_tg / ofct_srec は 1 スライス先行のパイプライン(cbp_pipe.h)で
    動く。投影バッファを 2 面持ち、GPU がスライス z を再構成している間に
    ホストはスライス z+1 の準備(投影コピー、黒投影判定、専用 CUDA ストリーム
    でのリング除去)を済ませ、z の TIFF コピー・書き出しは z+1 の GPU 計算と
    重なる。このために cbp.cu に AllocCBPProjection() / SelectCBPProjection()
    を追加(既存の CBP() / BeginCBP() / EndCBP() は不変)。角度表の転送を FFT
    より前に移し、BeginCBP() がフィルタ処理の完了を待たずに戻るようにした。
    ct_rec_g_* を多数同時実行する場合の注意:
    リング除去のデバイスバッファ(3 × Nx × Nt × 4 バイト)がプロセス終了
    まで確保されたままになる。
  ・p_rec: 2 件の修正と hp_tg と同じ高速化。(1) 読み込みループが最初の
    p ファイルを読み飛ばし、投影 0 が常に空(Nt 枚中 1 枚欠け)になっていた。
    全ファイルを使うよう修正し、結果はその 1 投影分だけ変わる。番号は 0 始まり
    でも 1 始まりでもよい。(2) 異常値除外が float に int 版 abs() を適用して
    いたため NaN/Inf を除外できず、前スライスの値が残ることがあった。fabsf()
    で判定し除外画素は 0 にする(しきい値 100 は不変、通常データは影響なし)。
    読み込みは並列(HPTG_READ_THREADS)でバンドに必要な行のストリップだけを
    読み、リング除去はインプレース、出力バッファは使い回し、単位変換と
    min/max は OpenMP、TIFF 書き出しは別スレッド、hp_tg と同じ 1 スライス
    先行パイプライン(HPTG_PIPELINE=0 で逐次順)。(1) を除いて旧版と CPU 版・
    GPU 版ともビット単位で一致することを確認済み。

【ver 2.5 の変更点】
  ・投影像生成 ct_prj_f に Paganin 位相回復(単一距離・平行ビーム)をオプション
    追加。第3引数にパラメータファイル pr.par を与えたときだけ有効(6h 参照)。
    省略時は従来どおり吸収のみで、結果も従来と完全に同じ。
  ・pr.par の P_SIZE(実効画素サイズ)はマイクロメートル単位。旧 mkpms/mkpmbg
    系は cm 単位だったため、旧ファイルの流用は 4 桁の取り違えになる。1 nm
    未満の値は cm 単位の疑いとしてエラー停止する。
  ・境界はミラーパディング(幅 = 核の減衰長 rc = sqrt(R2・δ/μ) の 10 倍、
    最低 32 画素)。旧 mkpms + mkpmbg + pm_sino 連鎖と背景減算後で一致する
    ことを合成データで検証済み(差 ~4e-7、32bit float の量子化と同水準)。
    詳細は 20260909_paganin_phase_retrieval.md。
  ・FFT は自前の radix-2(単精度・実数パック・OpenMP 並列)で外部ライブラリ
    不要。2048x2048 の投影 1 枚あたり約 0.1 秒(16 スレッド)。スレッド数は
    環境変数 PAGANIN_THREADS(既定は論理 CPU 数、OMP_NUM_THREADS とは独立)。
  ・TIFF連番 → OME-Zarr 変換 tif2zar を新設(6s)。再構成済みボリュームを
    Fiji + BigDataViewer/MoBIE 向けの多重解像度・blosc-zstd 圧縮形式へ
    1 パスで変換する。詳細は 20260909_tif2zar.md。

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
  ・両変数とも不正値は既定に戻す。正常な実行では何も表示しない。有効値の表示は
    黒投影が発生した(=両しきい値を超えた)場合のみ、床止め画素数・最小信号・
    最大投影値の要約は床止めが1画素以上あった場合のみ出力する。
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
      拒否する(この保護は無効化できない)。正常な実行では何も表示しない。
      有効値は黒投影が発生した場合のみ、床止め画素数・最小信号・最大投影値の
      要約は床止めが1画素以上あった場合のみ表示される。

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
      ct_rec / hp_tg / act_tg / p_rec / ofct_rec / ofct_srec / ofact_srec /
      sf_rec / rec2rec の
      すべてに効く。ただし rec2rec の入力は再構成円の外で端が減衰しているので、
      この場合は 0.1-0.2 程度の小さい値が適当。

   g. その他の環境変数
      上記のほかに、投影の並列読み込みスレッド数(HPTG_READ_THREADS、既定16。
      hp_tg / ofct_srec で有効)、メモリに合わせたチャンク分割の指定
      (HPTG_MEM_FRACTION / HPTG_MEM_LIMIT_MB / HPTG_CHUNK_ROWS)、使用する
      GPU 番号(CUDA_GPU、既定0)、入力ファイル名の上書き(RHP_O / RHP_D / RHP_Q、
      his_rec / his_srec / his_DO では RHIS_H / RHIS_C / HIS_SREC_OY)
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
      (p ファイルの番号は 0 始まりでも 1 始まりでもよく、全ファイルを使う)
      
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

   f. 自動CT用の連続再構成(16/8bit 直接出力)
      act_tg_g_r HiPic Dr RC RA0 rec16 LACmin16 LACmax16 rec8 LACmin8 LACmax8
      act_tg_g_r HiPic Dr L1 C1 L2 C2 RA0 rec16 LACmin16 LACmax16 rec8 LACmin8 LACmax8

      HiPic, Dr, RC / L1 C1 L2 C2, RA0: hp_tg(2b)と同じ
      rec16: 16bit 画像 rh%05d.tif の出力ディレクトリ
      LACmin16 / LACmax16: 16bit 画像の規格化範囲
      rec8: 8bit 画像 ro%05d.tif の出力ディレクトリ
      LACmin8 / LACmax8: 8bit 画像の規格化範囲

      再構成までは hp_tg と同一。32bit の rec*.tif の代わりに、各スライスを
      引数の 2 組の範囲で 16bit と 8bit に規格化して同時に書き出す(両ディレク
      トリは計算前に作成すること。各 LACmin は LACmax より小さいこと)。量子化・
      偶数サイズへの切り詰め・8 欄の ImageDescription は tif_f2i(4章)と同じ
      規則で、「hp_tg → tif_f2i」で作ったファイルと同一になる。hp_tg の環境変数
      はすべて有効。ビルドは GPU 版のみで、フィルタの接尾辞は hp_tg と同じ
      (act_tg_g_r / act_tg_g_s / act_tg_g_c)。

      *) hp_tg と同じく、HiPic ディレクトリの一つ上で実行する。

   g. a.his からの直接再構成(his_rec / his_srec / his_DO / his_sg)
      his_rec_P_F  layer {center} {pixel size} {offsetangle}
      his_srec_P_F HiPic Dr RC RA0 rec
      his_srec_P_F HiPic Dr L1 C1 L2 C2 RA0 rec
      his_DO       HiPic {stride}                 (CPU)
      his_DO_g     HiPic {navg {stride}}          (GPU)
      his_sg       layer {center} {skip}          (CPU)
      his_tb_rec_P_F / his_tb_srec_P_F            (引数は上と同じ)
      his_srec_g_F_n2i HiPic Dr RC RA0 rec K      (GPU。L1 C1 L2 C2 形式も可)

      入力はビームラインから持ち帰ったままの a.his(全撮影フレーム)・conv.bat・
      output.log が置かれた 1 つのディレクトリ。his_rec と his_sg はそのディレ
      クトリ内で(ct_rec と同様)、his_srec と his_DO は 1 つ上の階層でディレクトリ名を
      HiPic に指定して(hp_tg / ofct_DO と同様)実行する。his2img と conv.bat は
      実行せず、.img ファイルも作られない。conv.bat はメモリ上で次のように
      解釈する。
        his2img <his>                 his ファイル名(既定 a.his)
        img_ave <frames>... dark.img  暗電流 = 列挙フレームの平均
        img_ave <frames>... q<N>.img  q<N> = 列挙フレームの平均
        ren <frame> q<N>.img          q<N> = そのフレーム
        copy q<M>.img q<N>.img        q<N> = q<M> と同じデータ
        hp2do D000 I000 E000 D180 I180 E180 ...   his_DO 用に記憶
        rem / echo / pause            無視
      <frame> は his2img の命名 a<数字>.img(a0001.img がファイルの先頭フレー
      ム)。これ以外の命令、ファイル外のフレーム番号、同じ q の二重定義、
      dark.img の未定義はエラーで停止する。したがって平均枚数や入射光群の
      時間位置は conv.bat と output.log から読み取られ、固定値は持たない。

      180° スキャンか 360° オフセットスキャンかは output.log の角度欄で判定
      する(試料像が 350° に達すればオフセット。ofct_rec が従来から用いている
      判定)。
        プログラム  180° スキャン      360° オフセットスキャン
        his_rec     ct_rec(2a)         ofct_rec(3c)
        his_srec    hp_tg(2b)          ofct_srec(3b)。Rc = RC、Oy = 0、
                                       rangeList "-"、または L1 C1 L2 C2 形式
        his_DO      conv.bat の hp2do 行の画像で hp2DO
                                       ofct_DO / ofct_DO_g(3a)
        his_sg      sinog(6b)          of_sinog(6b)。center で貼り合わせ
      center を省略すると his_rec は再構成するレイヤーから回転軸を推定する。
      180° スキャンでは ct_rec 自身の推定、オフセットスキャンでは ofct_DO の
      探索を 1 次元にしたもの(そのレイヤーの -log 透過率を同じガウス平滑化
      OFCT_DO_SMOOTH で平滑化し、対向像を鏡像反転して、ofct_DO と同じ範囲の
      ずらし量について全ビュー対の平均二乗差を累積、放物線で補間。平滑化の
      幅だけ重なりの両端を除外する。含めると最小値が偏る)を用い、小数の推定値を
      そのまま使う。読むのはそのレイヤーの行だけなので数秒で済む。推定値と実行した
      ofct_rec のコマンドを標準エラーに表示し、最小値が探索範囲の端にある
      (回転軸が 0.75*Nx～Nx-1 の外。ofct_DO でも扱えない)ときは警告する。
      his_srec の引数は常に hp_tg 形式。オフセット時の RC は小数も可(3b)、
      Oy は環境変数 HIS_SREC_OY がなければ 0。L1 C1 L2 C2 形式は ofct_srec の
      傾き軸形式(3b)に渡され、軸位置をスライスごとに補間し、全スライスを
      最大幅で出力する。his_DO は
      180° スキャンでは hp2do 行の 6 画像(暗電流・入射光・0° の試料像、暗電流・
      入射光・180° の試料像。hp2do 行がなければ ct_rec の中心推定と同じ q1 q2
      と最終記録の次の 2 枚)で hp2DO の推定を行い(stride / navg 引数は無視)、
      hp2DO と同じ出力を標準出力に書く。オフセットスキャンでは ofct_DO。
      いずれも実行可能な his_srec のコマンド例を標準エラーに表示する(1.0 は
      画素サイズ Dr に置き換える)。
      his_sg は 1 レイヤーのシノグラムを sinog とまったく同じ形で書く(s%05d.tif、
      32bit float、../cmd-hst.log への 1 行)。オフセットスキャンでは
      of_sinog が center で両半分を貼り合わせる。center は小数も可で、"-" または
      省略なら上記の his_rec と同じ単層推定を使い、実行した of_sinog のコマンドを
      標準エラーに表示する。180° スキャンでは center は不要で無視される(skip を
      与えるときは "-" と書く)。
      平均は img_ave と同じ演算(整数和を枚数で切り捨て除算)なので、暗電流・
      入射光データも再構成結果もシノグラムも ct_rec / hp_tg / ofct_rec /
      ofct_srec / sinog / of_sinog の手順とビット一致し(his_srec は TIFF の
      ARTIST タグのみ異なる)、his_DO の
      出力は hp2DO / ofct_DO と同一になる。his は 16bit(type 2)のみ対応。
      環境変数は RHIS_H(his ファイル名。his2img 行より優先)、RHIS_C(conv.bat
      名、既定 conv.bat)、RHP_O(output.log のパス)、RHIS_SCAN(Z スキャン、
      後述)、RHIS_TB(試料像の平均を入射光にする。後述)、RHIS_VIEWS(his_srec の
      ビューの組 k/K。後述)、HIS_SREC_OY(上記)のほか、
      元のプログラムのものがすべて有効。フレームは 1 つのファイルハンドルへの
      位置指定読み出しで読むので、hp_tg の並列読み込み(HPTG_READ_THREADS)と
      メモリ量に応じたチャンク分割はそのまま働く。

      his_tb_rec と his_tb_srec は、入射光をデータ自身から取る his_rec / his_srec。
      output.log で入射光(it = 0)とされた q はすべて、output.log の最初の空行より
      前に並ぶ試料像(it = 1。空行の後の末尾の 1 枚は含まない)の整数平均(img_ave と
      同じ)として読まれる。したがって結果は、それらのフレームを img_ave で平均して
      入射光画像に上書きしてから ct_rec / hp_tg / ofct_rec / ofct_srec を実行した
      場合と一致する。暗電流の差し引きなど他は不変で、入射光の時間内挿は定数になる。
      his_tb_srec は平均を作るために a.his をもう 1 回全読みする。his_* のどの
      プログラム(his_DO、his_sg、HIS_REC / HIS_SREC で指定した Z スキャンのスクリ
      プト)でも環境変数 RHIS_TB=1 で同じ動作になり、his_tb_* では RHIS_TB=0 で
      切れる。起動時の要約行に「I0 = mean of N sample views (RHIS_TB)」、
      cmd-hst.log に「% RHIS_TB」が付く。

      his_srec_g_F_n2i は Noise2Inverse 用の his_srec。試料像を角度順に K = 2 または
      4 組に交互分割し(k, k+K, k+2K, …)、各組を his_srec として rec/k に再構成する。
      RA0 には組の最初のビューの角度(output.log から取る)を加えるので、K 枚の像は
      同じ幾何になる。RC と Dr は全組共通なので回転軸は明示的に与えること。試料像の
      枚数は K の倍数でなければならず、180° スキャンのみ(オフセットスキャンは
      エラー)。各組はこのプログラムの子プロセスとして順に実行され、子は環境変数
      RHIS_VIEWS=k/K を持つ通常の his_srec なので、手で RHIS_VIEWS を設定して
      his_srec を実行すれば 1 組だけ再構成できる。cmd-hst.log に「% RHIS_VIEWS=k/K」
      が付く。RHIS_SCAN と RHIS_TB は全組に効く。

      Z スキャン。a.his に L 枚のスキャンが M 本入っている場合(L は conv.bat が
      参照する最大フレーム番号、M = 総フレーム数 / L。spl / his_spl_E が同じ
      conv.bat・output.log とともに 001/raw～00M/raw に展開していた配置)、環境
      変数 RHIS_SCAN = 1～M で his_* 全プログラムのスキャンを選ぶ。RHIS_SCAN なし
      ではエラー停止し、誤って 1 本目を使うことはない。スクリプト rc-check-his
      と gen-all-his(Windows は exe/*.bat、Linux は bin/*.sh)が 6m の spl +
      conv.bat + ct_rec / hp_tg の手順を置き換える。
        his_info {dir}
            a.his + conv.bat + output.log の「Nx Ny frames L M layer kind」を 1 行で
            表示する: 画像サイズ、総フレーム数、1 スキャンのフレーム数 L、スキャン
            数 M、中央レイヤー Ny/2、スキャン種別 180 / 360(標準エラーには読みやすい
            要約)。
        rc-check-his {layer}
            スキャン数 M、layer の既定値(中央行)、スキャン種別は his_info から取る。
            全スキャンについてその layer を his_rec で軸推定つきで再構成し、軸を
            center.log(2 列目)に、画像を目視確認用に rc-check/00k.tif に書く。
            center.log の値は his_srec に渡す値で、360° スキャンでは
            ImageDescription から 0.5 を引く(ofct_rec は画像中心 dcnt+0.5 を記録し、
            ofct_srec は検出器上の軸位置 Rc = dcnt を取るため。ct_rec / hp_tg は
            規約が同じなので 180° では補正なし)。
        gen-all-his center.log Dr {LACmin LACmax}
            全スキャンを center.log の軸で his_srec により 00k/rec に再構成し、
            tif_f2i(8bit、既定 -0.5 3.0)で 00k/ro_xy を作る(ro_zx は作らない)。
            積層は従来どおり rec_stk(6f)。
        his_zrec scan layer {center}
            試し再構成 1 枚: そのスキャンの layer を his_rec で再構成し(center 省略
            時は推定)、rc-check/SSS_LLLLL_C.tif(撮影番号 %03d、レイヤー %05d、C は
            使った軸位置を his_srec に渡す値で表したもの。例 003_01024_1023.5.tif)
            として保存する。レイヤーや軸位置を
            変えた結果を並べて比較するためのもの。
      スクリプトはいずれも a.his のあるディレクトリで実行する(his_srec は HiPic
      に "." を取れる)。プログラム名は HIS_REC / HIS_SREC で変更できる(既定
      his_rec_g_c / his_srec_g_c)。

      *) his_rec と his_sg は a.his・conv.bat・output.log のあるディレクトリで、
         his_srec と his_DO はその 1 つ上で実行する。

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
       OFCT_DO_SMOOTH で σ を指定し、未設定の場合は 10.0。0 を指定すると平滑化しない。

   b. 再構成
      ofct_srec_P_F HiPic Rc Oy rangeList Dr RA0 rec
      ofct_srec_P_F HiPic L1 C1 L2 C2 Oy rangeList Dr RA0 rec
      ofct_srec_P_F HiPic Rc Oy
           (第 3 の形式は再構成せず、貼り合わせ幅 N、パディング後の長さ L、
            対向ビュー対の数 M と投影バッファの大きさ 4*N*L*M を GB で表示する。
            大きな計算の前の確認用)
      
      HiPic: q????.img もしくは q????.tif が格納されているディレクトリ名(/ は不要)
      Rc: 回転軸の位置(左端からの画素数)。小数も可(ver 2.8)。貼り合わせの枠は
           整数部で決め、対向像を 2f 画素だけ線形補間でずらして重ね、回転中心を
           f だけ動かす(180° スキャンと同様)。整数を与えれば従来と同一の結果。
      L1 C1 L2 C2: 傾いた回転軸(ver 2.8)。レイヤー L1 で位置 C1、L2 で C2 とし、
           各スライスは線形補間した位置を使う(小数。Rc の項を参照)。
           貼り合わせ幅 2*floor(Rc(z)) がスライスごとに変わるので、選択
           スライスのうち最大の幅で全スライスを出力し、狭いスライスは中央に
           置いてシノグラムの余白を各行の端の値で埋める(0 にすると、両端の列を
           見る打ち切り補正の判定が狭いスライスだけ常に無効になるため)。回転軸
           は常に画像中央にあり、そのまま積層できる。C1 = C2 なら第 1 の形式と
           同一の結果。画像の縁付近では、そのスライス固有の幅で再構成した場合と
           わずかに異なる(試験データで中心部は値域の 0.1 % 程度、縁で数 %)。
           L1、L2 は出力レイヤー番号(Oy = 0 なら検出器の行番号と同じ)。rangeList
           で L1～L2 の外を指定すると外挿になる。
      Oy: 縦ずれ量（常に 0 を指定）。
      rangeList: 再構成するレイヤーを指定。
      Dr: 画素サイズ
      RA0: 回転軸の原点オフセット
      rec: 再構成画像を出力するディレクトリ(計算前に作成すること)

      *) rangeListの例。
         100レイヤーだけの場合：100
         100-150 レイヤー: 100-150
         全部: -

      ImageDescription の 2 列目はそのスライスの Rc+0.5(ofct_rec の dcnt+0.5 と
      同じ規約)。軸が中央より右なら ver 2.8 より前の貼り合わせ画像の中心 (N+1)/2
      と同値で、既存の出力は変わらない。傾きや小数の軸ではスライス固有の値になる。

   c. img/tif データから1枚だけ再構成（オフセットCT）
      ofct_rec_P_F layer center {pixel size} {offsetangle}
      
      layer: 再構成するレイヤー(高さ)
      center: 回転軸の位置(pixel)。小数も可(ver 2.8)。対向像を 2f 画素ずらして
           重ね、回転中心 dcnt+0.5 が鏡映の軸に一致する。ver 2.8 より前は小数を
           与えると中心だけ動き、貼り合わせは切り捨てた整数で行われていた
           (不整合)。整数なら従来と同一の結果。
      pixel size: 画素サイズ(um)。省略した場合は1.0になる。
      offset angle: 回転軸の原点オフセット。省略した場合は0.0になる。
      
      *) q????.img もしくは q????.tif があるディレクトリで実行する。(dark.img / dark.tif で自動判別)
   


   d. 自動CT用のオフセットCT連続再構成(16/8bit 直接出力)
      ofact_srec_g_r HiPic Rc Oy rangeList Dr RA0 rec16 LACmin16 LACmax16 rec8 LACmin8 LACmax8

      再構成までは ofct_srec(3b)と同一。rec16/rh%05d.tif(16bit、LACmin16..
      LACmax16)と rec8/ro%05d.tif(8bit、LACmin8..LACmax8)を直接書き出す。
      量子化と ImageDescription は tif_f2i と同じ規則で、「ofct_srec → tif_f2i」
      と同一になる。引数 3 個の形(HiPic Rc Oy)は従来どおりメモリ見積もりのみ。
      ビルドは GPU 版のみで、フィルタの接尾辞は ofct_srec と同じ(ofact_srec_g_r / _s / _c)。

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
    スライスの画素数が揃っていない場合(オフセットCTでは回転中心位置に
    より再構成後の画素数が変わる)は、全スライスの最大サイズに揃える:
    小さいスライスは中心合わせで配置し、外側は画素値 0 で埋める(サイズ差
    が奇数のときは左・上側のパディングを 1 画素少なくする)。揃え処理を
    行った枚数は終了時に表示される。切り出し座標 x1..y2 はこの最大サイズ
    の座標系で指定する。

5. 32bit tiff 画像へのgaussian filter適用
    rec_gf rec radius out
    
    rec: 32bit tiff 再構成画像があるディレクトリ
    radius: gaussian filter の半値幅
    out: フィルター後の画像を出力するディレクトリ(実行前に作成すること)

6. その他おまけ
   a. tif 画像に埋め込まれたタグ(画素サイズ、投影数など)を表示する
      pid tiff-file
      
   b. シノグラム作成
      sinog    layer {skip}
      of_sinog layer offset {skip}
      32bit tiff にて出力。タグには最大と最小値のみ入る。
      skip を指定すると、投影数をskip分の1にして出力する。of_sinog は 360°
      オフセットスキャン(3)の両半分を回転軸位置 offset(ofct_rec の center と
      同じ。ver 2.8 から小数も可、整数なら従来と同じ結果)で貼り合わせる。
      貼り合わせたシノグラムの回転軸位置(6c 用)は offset+0.5。his_sg(2g)は
      同じことを a.his から直接行う。sino.tmp(幅・投影数・レイヤー・中心・
      全吸収の平均と SD。旧バッチ処理用)は ver 2.8 から書かない。

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
      rhro_stk num_stack start end {name ...}
      act_tg / ofact_srec の 16/8bit 出力用。各スキャンの rh/rh%05d.tif と
      ro/ro%05d.tif(レイヤー start～end)を whole/rh、whole/ro に 1 本の連番
      (スキャン j のレイヤー i → i + j*(end-start+1)、rec_stk と同じ)で移動する
      コマンドを出力する。name を省略すると rh と ro、指定すればその組だけ。
      例: rhro_stk 8 17 140 > ddd.bat

   g. his ファイルを連番の img に変換。
      his2img his-file (x1 x2 y1 y2)
      カッコ内の引数指定にてcropが可能
      
   h. 180deg スキャンから投影像を作成
      ct_prj_f HiPic prj {pr.par}
      入力は dark.img があれば img、無く dark.tif があれば tiff を自動判別。
      出力は prj/p%05d.tif(32bit float)。

      第3引数 pr.par を与えると、Paganin の単一距離位相回復(Paganin et al.,
      J. Microsc. 206, 33-40, 2002 の式 10、平行ビーム)を投影ごとに適用し、
      値は μt = -ln(T') になる。省略すると従来どおり吸収のみ。
      pr.par は次の 5 行(1 行 1 数値。# または // 以降はコメント可):
        MU      線吸収係数            [cm^-1]  > 0
        DELTA   屈折率デクリメント δ  [-]      > 0
        R1      光源-試料距離         [cm]     > 0(平行ビームでは未使用)
        R2      試料-検出器距離       [cm]     >= 0
        P_SIZE  実効画素サイズ        [um]     > 0
      【注意】P_SIZE はマイクロメートル単位。旧 mkpms/mkpmbg 系(cm 単位)の
      pr.par を流用しないこと。単位を 4 桁間違えると平滑化が過多/過少になる。
      極端な値は起動時にエラーで止まるが、止められない場合もあるので、起動時
      に表示される rc(平滑化半径、画素単位)が物理的に妥当(通常は数画素)で
      あることを確認すること。
      低透過率ガード(1e)は位相回復後の透過率にも同様に効く。スレッド数は
      環境変数 PAGANIN_THREADS で指定(既定は論理 CPU 数)。

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

      a.his から直接(2g)、spl と conv.bat なしで行う場合(スキャン数は a.his
      から読む。layer を省けば中央行):
      rc-check-his 75
      gen-all-his center.log 11.31
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
      Linux での同時実行数は環境変数 ACT_SPL_JOBS(既定16)で制限し、終了前に
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

      

   s. TIFF連番 → OME-Zarr 変換
      tif2zar tifDir out.zarr {--chunk N} {--pixel P} {--min A --max B}
                              {--levels L} {--codec C} {--clevel N}

      再構成済みボリューム(uint16 / 32bit float の連番 tiff)を OME-Zarr
      v0.4(多重解像度ピラミッド + blosc-zstd 圧縮)へ変換する。Fiji の
      Plugins → BigDataViewer → N5 Viewer 等で file:///D:/... 形式の URI
      を与えて開く(C:\ の直接入力は Fiji 側の既知の不具合で落ちる)。
      float32 入力はタグの min/max から全体レンジを決めて uint16 に量子化
      する。rec 直変換では外れ値にレンジが引きずられるため、tif_f2i と
      同じ規格化値を --min/--max で与えることを推奨。連番の欠番はエラー
      停止(欠損ファイル名を列挙)。スライスの画素数が揃っていない場合
      (オフセットCT)は全スライスの最大サイズに中心合わせで揃え、外側を
      画素値 0 で埋める(tif_f2i の 4 章と同じ規則)。スレッド数は環境変数 TIF2ZAR_THREADS
      (既定は論理 CPU 数)。詳細・実測性能・出力仕様は 20260909_tif2zar.md
      を参照。

      例: tif2zar rec rec.zarr --min -0.5 --max 3.0
