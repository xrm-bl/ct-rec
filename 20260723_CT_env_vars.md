# CT再構成ソフトウェア 環境変数一覧

xrm-bl/ct-rec `src/` の getenv 横断調査に基づく / **2026-07-27 版**（初版 2026-07-23、随時更新）
（`20260620_CT_env_vars_EN/J.pdf` を置き換える更新版。差分は末尾参照）

---

## スレッド / 並列数の制御

| 変数 | 用途 | 既定 | 使用ファイル (src/) | 備考 |
|---|---|---|---|---|
| `CBP_THREADS` | CPU版CBP(畳み込み+逆投影)のワーカースレッド数 | **実行PCの論理コア数−1 (下限1)** | cbp_thread.c cbp_thread_int.c cbp_thread_nai.c cbp_thread_avx.c | 2026-07 変更(旧: 固定8)。`-DCBP_THREADS=n` でコンパイル時固定も可。1〜M(投影数)の範囲外は `bad number of CBP_THREADS` で停止 |
| `HPTG_READ_THREADS` | hp_tg / ofct_srec の投影**並列読み込み**スレッド数 | 16 | hp_tg_ku.c ofct_srec.c | 2026-07 新設。InfiniBand等のディスクアレイは並列読みで帯域が出る。投影数超は自動クランプ |
| `THREADS` | 回転軸探索の比較スレッド数 | 40 | oct_xy.c otf_xy.c ofct_DO.c | 3本とも既定40。0以下は `bad number of THREADS` |
| `OMP_NUM_THREADS` | リング除去ソートフィルタの並列数 | 40 | sort_filter_omp.c sort_filter_g.cu | `get_num_threads_from_env()` が読む。**未設定時はOpenMP既定ではなく40を返す**。0以下は無効→40 |
| `ACT_SPL_JOBS` | act_spl / act_spl2 (Linux) の背景ジョブ同時実行数 | 8 | act_spl.c act_spl2.c | 2026-07 新設。Linuxで gf_sd / tif_mgf を fork+exec 背景実行する際の上限(Windows は従来どおり `start /b`)。終了前に全ジョブを待ち合わせる |

※ GPUソートは内部 `SORT_THREADS=256`、GPUフィルタは `BLOCK_SIZE_X/Y/Z=8/8/4` の固定値(環境変数ではない)。

## リング除去

| 変数 | 用途 | 既定 | 使用ファイル (src/) | 備考 |
|---|---|---|---|---|
| `KERNEL_SIZE` | リング除去のソート/メディアンフィルタのカーネル幅 | 5 | sort_filter_omp.c sort_filter_g.cu | 判定は `size > -1`(0以上を許容) |

## 再構成品質 / 前処理

| 変数 | 用途 | 既定 | 使用ファイル (src/) | 備考 |
|---|---|---|---|---|
| `CT_REC_BLACK_THRESH` | **画素ごと**の信号下限 [dark を超えたカウント数]。これ以下の画素は不透明とみなし、この値で床止めして対数変換する | **2.0** | blacklim.h ct_rec.c ct_rec_c.c ofct_rec.c otf_rec.c tf_rec.c sinog.c of_sinog.c ct_prj_f.c ct_prj_f_c.c ict_prj_fc.c rhp.c rhp_c.c rtf.c | **2026-08 に意味を変更**(旧: ラインの平均信号がこれ未満の投影を黒とみなす)。単位は従来と同じなので readme の 1/10/100/1000 のはしごはそのまま使える。床止めにより `log(I0/(I−dark))` が +Inf / NaN になり得なくなる(→GPUリング除去のクラッシュ根絶)。投影値の頭打ちは `log(I0/しきい値)`。**0以下は不正値として拒否**し既定に戻す(無効化はできない)。有効値は起動時に stderr へエコー。rhp_c.c は ReadHiPic / ReadHiPicBand の両方で参照(共通リーダ経由で hp_tg / ofct_srec / ofct_DO にも効く) |
| `CT_REC_BLACK_FRAC` | 黒(missing angle)投影の判定。`CT_REC_BLACK_THRESH` で床止めされた画素の割合がこれを超えた投影を黒とみなし、良い投影の平均で置き換える | 0.5 | 同上 | **2026-08 新設**。旧来の「平均信号 < しきい値」を置き換える。一様に減衰したライン(板状試料の edge-on)では両判定が同じ透過率で発火するため板状試料の挙動は不変。視野に素抜けの余白が残る場合、平均には (1−被覆率)×信号 のハードフロアがあり旧判定は原理的に発火できなかった。有効域 0 < v ≤ 1、範囲外は既定維持 |
| `OFCT_DO_SMOOTH` | ofct_DO(_g) の SSD 計算前ガウシアン平滑化 σ | **10.0**（**既定で有効**） | ofct_DO.c ofct_DO.cu | **0 で無効化**。分離型2次元等方ガウシアン(水平→垂直2パス)。2026-08-06 に既定を 1.0→10.0 へ変更(それまで CPU 版 1.0 / GPU 版 2.0 と食い違っていたのを統一) |
| `PAD_THRESH` | **全再構成ソフト共通**(CBP層)の打ち切り(カッピング)補正しきい値(**比率**) | **0.3(既定ON・自動判定)** ※未設定時。既定値は `cbp.h` の `PAD_THRESH_DEFAULT` | cbp_thread.c cbp_thread_int.c cbp_thread_nai.c cbp_thread_avx.c cbp.cu | シノグラム端列の平均振幅が全体平均の この比率を超えたら(=試料が視野をはみ出していたら)フィルタ入力を端値ホールド+コサイン減衰で幅0.5N外挿(逆投影コスト不変、FFT長のみ2倍)。視野内試料は端≈0で自動的に無適用。ct_rec/hp_tg/ofct_rec/ofct_srec/p_rec/sf_rec/tf_rec/rec2rec すべてに有効。**`PAD_THRESH=0`(以下)で強制OFF**(旧既定と完全一致)。rec2rec入力(再構成円で端が減衰)は 0.1〜0.2 推奨。既定値変更は `cbp.h` の1定義、または `-DPAD_THRESH_DEFAULT=…` |

## メモリ対応チャンク分割 (3プログラム共通)

| 変数 | 用途 | 既定 | 使用ファイル (src/) | 備考 |
|---|---|---|---|---|
| `HPTG_MEM_FRACTION` | 使用を許す空きメモリの割合(バンドサイズ算定) | **hp_tg: 0.95 / ofct_srec, p_rec: 0.9** | hp_tg_ku.c ofct_srec.c p_rec.c | 有効域 0.05 < v ≤ 0.95、範囲外は既定維持。**同一ノードで複数プロセスを並走させる場合は 0.8/プロセス数 程度に下げること**(全員が9割を要求するとスワップする) |
| `HPTG_MEM_LIMIT_MB` | メモリ上限の絶対値(MB) | 未設定=割合ベース | hp_tg_ku.c ofct_srec.c p_rec.c | mb>0 のとき budget を mb×1024×1024 に置換(HPTG_MEM_FRACTION より優先) |
| `HPTG_CHUNK_ROWS` | 1バンドあたりのスライス行数を手動指定 | 未設定=自動 | hp_tg_ku.c ofct_srec.c p_rec.c | atoi>0 でその行数に強制(1〜total_rows にクランプ)。3変数の優先順位: CHUNK_ROWS > MEM_LIMIT_MB > MEM_FRACTION |

## 入出力 (HiPic/TIFF 共通リーダ)

| 変数 | 用途 | 既定 | 使用ファイル (src/) | 備考 |
|---|---|---|---|---|
| `RHP_O` | 読み込む output.log のパス上書き | output.log | rhp.c rhp_c.c rtf.c | `-` で stdin |
| `RHP_D` | ダーク画像ファイルのパス上書き | dark.img / dark.tif | rhp.c rhp_c.c rtf.c | 拡張子で img/TIFF を判別 |
| `RHP_Q` | 投影ファイルの接頭文字の上書き | q | rhp.c rhp_c.c rtf.c | 先頭1文字(a〜z)。例 `q[0-9]*.tif` の "q" を変更 |

## GPU

| 変数 | 用途 | 既定 | 使用ファイル (src/) | 備考 |
|---|---|---|---|---|
| `CUDA_GPU` | 使用する CUDA デバイス番号 | 0 | cu.h (SETUP_CUDA_GPU) | GPU各プログラム(cbp.cu 経由)で使用。負値/不正は "no device assigned" |

---

## 2026-06-20 版からの差分

| 変数 | 変更内容 |
|---|---|
| `CBP_THREADS` | 既定 8 → **論理コア数−1(下限1)**。対象に cbp_thread_avx.c(新設・選択制)追加 |
| `HPTG_READ_THREADS` | **新設**(既定16)。hp_tg / ofct_srec の投影並列読み込み |
| `ACT_SPL_JOBS` | **新設**(既定8)。act_spl / act_spl2 の Linux 背景ジョブ上限 |
| `OFCT_DO_SMOOTH` | **新設**(既定1.0、0で無効)。ofct_DO(_g) の前処理平滑化 |
| `PAD_THRESH` | **新設**(既定0.3=ON自動判定、`cbp.h` の `PAD_THRESH_DEFAULT`)。CBP層共通の打ち切り(カッピング)補正。全再構成ソフトに有効。`PAD_THRESH=0` で強制OFF |
| `HPTG_MEM_FRACTION` | 既定 0.8 → **0.9**(ofct_srec, p_rec)。hp_tg のみ **0.95** |
| `THREADS` | ofct_DO の既定が 8 → **40**(3本とも40に統一) |
| `OMP_NUM_THREADS` | 未設定時の実挙動は「OpenMP既定」ではなく **40** と判明(記載修正) |
| `CT_REC_BLACK_THRESH` | 参照ファイルが拡大: ct_rec_c.c / ofct_rec.c / rhp_c.c / rtf.c を追記 |

## 2026-08-06 の変更 (低透過率ガード)

| 変数 | 変更内容 |
|---|---|
| `CT_REC_BLACK_THRESH` | **意味を「ラインの平均」から「画素ごとの床」へ変更、既定 1.0 → 2.0**。単位(dark を超えたカウント数)は不変。新しい共通ヘッダ `src/blacklim.h` に実装を集約し、参照ファイルを 13 本に統一 |
| `CT_REC_BLACK_FRAC` | **新設**(既定 0.5)。黒(missing angle)投影の判定を、床止めされた画素の割合で行う |

背景: 試料が厚い/高密度で `I − dark` が 0 以下になると `log(I0/(I−dark))` が +Inf / NaN を返し、
GPU リング除去 (`sort_filter_g.cu`) の列ソートが詰め物の番兵 (index −1) を実ランクへ押し込むため
`perm[]` に 4294967295 が入り、`median_scatter` が確保領域の遥か外へ書いて
`CUDA error sort_filter_g.cu:266: an illegal memory access was encountered` /
`ring removal image processing failed` で停止していた。CPU/OpenMP 版は停止しない代わりに
当該列のソート結果が黙って壊れていた。画素ごとの床止めにより両方が構造的に解消する。

旧判定が機能しなかった理由: 判定量がラインの**平均**で、危険を決めるのは**最小値**(1画素でも 0 以下なら
クラッシュ)だったため。さらに視野に素抜けの余白が残ると平均には (1−被覆率)×信号 のハードフロアが生じ、
被覆率 60% 程度では中心の透過率を 0 にしても既定 1.0 では発火しない。両端 10 画素を空気基準に使う
処理(`p_ave`)がそもそも余白の存在を前提としているため、このフロアは設計上必ず存在する。
| `RHP_O/D/Q` | rhp_c.c(現行リーダ)を追記 |

調査方法: `grep -rnoE 'getenv\("[A-Z_0-9]+"\)' src/*.c src/*.cu src/*.h` で横断し、各既定値をソースで確認。
