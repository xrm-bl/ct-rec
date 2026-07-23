# CT再構成ソフトウェア 環境変数一覧

xrm-bl/ct-rec `src/` の getenv 横断調査に基づく / **2026-07-23 版**
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
| `CT_REC_BLACK_THRESH` | 黒(ダーク)投影判定しきい値。平均信号 (T−dark) がこれ未満の投影を黒とみなし補正 | 1.0 | ct_rec.c ct_rec_c.c ofct_rec.c otf_rec.c tf_rec.c rhp.c rhp_c.c rtf.c | rhp_c.c は ReadHiPic / ReadHiPicBand の両方で参照(共通リーダ経由で hp_tg / ofct_srec / ofct_DO にも効く) |
| `OFCT_DO_SMOOTH` | ofct_DO(_g) の SSD 計算前ガウシアン平滑化 σ | 1.0（**既定で有効**） | ofct_DO.c ofct_DO.cu | **0 で無効化**。分離型2次元等方ガウシアン(水平→垂直2パス) |

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
| `HPTG_MEM_FRACTION` | 既定 0.8 → **0.9**(ofct_srec, p_rec)。hp_tg のみ **0.95** |
| `THREADS` | ofct_DO の既定が 8 → **40**(3本とも40に統一) |
| `OMP_NUM_THREADS` | 未設定時の実挙動は「OpenMP既定」ではなく **40** と判明(記載修正) |
| `CT_REC_BLACK_THRESH` | 参照ファイルが拡大: ct_rec_c.c / ofct_rec.c / rhp_c.c / rtf.c を追記 |
| `RHP_O/D/Q` | rhp_c.c(現行リーダ)を追記 |

調査方法: `grep -rnoE 'getenv\("[A-Z_0-9]+"\)' src/*.c src/*.cu src/*.h` で横断し、各既定値をソースで確認。
