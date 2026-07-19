# ct-rec 高速化・libtiff更新 作業ログ (2026-07)

作業期間: 2026-07-17 〜 2026-07-19
対象コミット: `ecde4be` 〜 `384f738`
検証環境: Windows 10 / MSVC (VS2022) / CUDA 13.2 / RTX 3070 (8GB, sm_86) / 8コア16スレッド

このセッションで行った libtiff の更新と一連の高速化、および「検討したが採用しなかった」
事項をまとめる。数値はすべて合成データでの実測値（best-of-N、ローカルSSD）。

---

## 0. 全体サマリ

| 項目 | 効果 | コミット |
|---|---|---|
| libtiff 3.6.0 → 4.6.0 移行 | 互換維持・最新化 | `ecde4be` |
| 保存高速化 (raw strip 複数行化) | **保存 ×3.6** | `1607fd8` |
| rif_f 高速読み込みパス | 読み込み ×18〜23（現状未使用経路） | `09289db` |
| CBP GPU 転送最適化 | **再構成/スライス ×8** | `1c2aed7` |
| CBP CPU 2フェーズ化 | **CPU再構成 ×2**・決定的出力・メモリ1/K | `07c1670` |
| CBP AVX2 変種 (選択制) | 精度200倍・速度はint版に劣る→既定はint | `8c317aa` |
| tif_mgf GPU版 | median+gaussian の GPU化 | `337eb0b` |
| ofct_DO_g ペア平均 | 軸探索のS/N向上 | `8d3ff2b` |
| HPTG_MEM_FRACTION 0.8→0.9 | チャンク既定調整 | `4a982e8` |
| 起動時サイズ検証の間引き | **起動 10〜20秒 → 大幅短縮** | `1eb6872` |

---

## 1. libtiff 3.6.0 → 4.6.0 移行 (`ecde4be`)

詳細は `20260717_LIBTIFF-MIGRATION.md` 参照。要点:

- 4.4.0 で削除された旧 typedef を標準の `<stdint.h>` 型へ (`uint32`→`uint32_t` 等)。
- Windows は vcpkg で `x64-windows-static`(/MT) をビルドし、`src/` に配置:
  `libtiff.lib` + 依存 `jpeg.lib`/`lzma.lib`/`zs.lib`(zlib) + ヘッダ4点。
- **CRTモデルは /MT**（既定 cl も nvcc も /MT）。バッチのコンパイルフラグ変更は不要。
- 旧3.6.0一式は `src/libtiff-3.6.0/` に退避。
- クラシックTIFF(8/16/32bit)の入出力は完全に互換。

---

## 2. TIFF 保存の高速化 — raw strip 複数行化 (`1607fd8`)

### 背景
再構成/変換系の writer は `TIFFWriteRawStrip` を「1行=1ストリップ」でループしていた
(`RowsPerStrip=1`)。1行ごとに strip とwrite syscall が発生し遅い。

### 実装
`src/tifwrite.h` に `ct_write_raw_strips()` を新設。約1MB/ストリップになる行数を
scanline幅から算出し、連続する複数行を1回の `TIFFWriteRawStrip` で書く。
書き込むバイト列は同一（連続行）なので**出力画素はビット不変**。24 writer を置換。

### ベンチ (4096×3072 float, 50MB, ローカルSSD)
ストリップ戦略の谷は **1〜4MB**。全ライン一括はむしろ遅い（下表）。

| 方式 | rps | 時間 | 対1行 |
|---|---|---|---|
| Scanline rps=1 (旧) | 1 | 49.4 ms | 1.0× |
| Scanline rps=64 (~1MB) | 64 | 15.9 ms | 3.1× |
| **Scanline rps=WHOLE** | 3072 | 42.1 ms | 1.2×(遅) |
| EncStrip rps=256 (~4MB) | 256 | 13.7 ms | 3.7× |

→ ~1MB 目標 (`TIFFWriteScanline`相当) を採用。実測 **49.4 → 13.7 ms (×3.6)**。

### 重要な注意 (再発防止)
`TIFFWriteRawStrip` は **strip番号とバイト数が RowsPerStrip と一致していないと壊れる**。
「RowsPerStripだけ変えて書き込みループを直さない」と出力が破損する（実際に一度この
状態を作り、rps=1へ戻してから両者を同時に直した）。scanline系 writer
(`tif_mgf`/`hp2DO`/`tf2DO`) は `TIFFWriteScanline` なので rps>1 でも安全、
`TIFFDefaultStripSize` を使用。

---

## 3. rif_f 高速読み込みパス (`09289db`)

`ReadImageFile_Float` は1サンプルを関数ポインタ多段(最終 `fgetc`)で読み、
4096×3072 float で約0.9秒かかっていた。無圧縮・通常fill order・バイト境界サンプル
の場合に「ストリップ丸ごと `fread`→一括変換」のファストパスを追加（LZW/PackBits/
予測子等は従来経路へフォールバック）。**×18〜23**、旧版とビット一致。

注意: 現行ツール (`tif_f2i`/`rec_crop`/`rec_gf`) は本リーダをヘッダ取得のみ
(`cell=NULL`) で使い、画素は各自の libtiff scanline で読むため、**この高速化の
恩恵を受ける実行経路は現状ない**（将来 `cell!=NULL` で使う時のための健全な実装）。
旧実装は `rif_f.c.org` に保管。

---

## 4. CBP (convolution back projection) 高速化

CBP は re-構成の律速。GPU/CPU 両実装を最適化した。

### 4.1 GPU: cbp.cu 転送最適化 (`1c2aed7`)
フェーズ別計測 (N=2048,M=1800) で **転送が全体の89%**、カーネルは既に軽いと判明:
prepare 65ms(投影ごとM回memcpy) + execute 17ms + end 70ms(行ごとN回memcpy)。

- **PrepareCBP**: `cudaMemset`1回 + `cudaMemcpy2D`1回で全投影転送。ホストバッファを pinned化。
- **EndCBP**: 連続領域を1回の D2H で直接受け取り (Float==float)。
- **InitCBP/ExecuteCBP**: `L2*M>CUFFT_LIMIT` の退化フォールバック(投影ごと単発FFT+
  D2Dコピー)を廃止し、常時チャンクバッチFFTに。`CUFFT_LIMIT` は1チャンクの上限要素数に意味変更。
- `Float!=double` ビルド用に旧経路をフォールバック温存。

実測: 152.2 → **19.0 ms/スライス (×8.0)** (N=2048)、64.4 → 3.0 ms (N=1024)。
出力は旧実装と**ビット一致**(max abs diff=0)、実データでも一致確認済み。

### 4.2 CPU: cbp_thread_int.c 2フェーズ化 (`07c1670`, C1+C2)
- **C1**: 既定スレッド数を「実行PCの論理コア数−1(下限1)」に (旧: 定数8)。
  環境変数 `CBP_THREADS` / `-DCBP_THREADS=n` で上書き可。3変種すべてに適用。
- **C2**: スレッド分割を「投影分割の1パス」から「A:投影分割でFFT畳み込み → B:画像
  行帯分割で逆投影」の2フェーズに再構成。
  - スレッドごとのN²画像複製(K倍メモリ)と直列合算を廃止 → **メモリ1/K**
    (K=15,N=4096で約1GB→67MB)。
  - 量子化スケールを全投影のグローバル最大に統一 → **出力がスレッド数に依存しない
    (決定的)**。
  - 実測 1862 → **約900〜1130 ms (×2)**。旧版(K=1)と**ビット一致**、K=4/15間も一致。

### 4.3 CPU: cbp_thread_avx.c 追加 (選択制) (`8c317aa`)
AVX2/FMA・float の逆投影変種。GPU版と同じ直接評価 (r=Xcosθ+Ysinθ-r0) を8画素同時、
内部ブロックは連続load+permute、端はマスク付きgather。
- 精度: double基準比 max abs diff 5e-9（**int版の約200倍良い**、実質倍精度一致）。
- 速度: このデスクトップ(8コア)では **int版より3〜5割遅い**（gather/permuteコスト）。
  → **既定は int版のまま**。将来のXeon等 gather高速CPUで再評価する余地あり。
- 注意: GCCビルドの `#pragma GCC target("avx2,fma")` は **Linux未検証**。

---

## 5. その他の機能追加・調整

- **tif_mgf_g.cu** (`337eb0b`): `tif_mgf` (2D median+gaussian) の GPU版。単一ファイル・
  同一CLI。median は整数演算でCPU版とビット一致、gaussian は単精度(≤1 LSB差)。
  `cuda13_compat.h` 準拠、メタデータはガード付きコピー(NULL-ImageDescription回避)。
  `makefileGPU.bat` / `MakefileGPU` の `filters_g` に登録。
- **ofct_DO_g** (`8d3ff2b`): 第2引数 `group`(既定10)。K個の投影ペアを透過率ドメインで
  平均してから -log することで低S/N時の noise×noise を抑え、回転中心推定を安定化。
- **HPTG_MEM_FRACTION 0.8→0.9** (`4a982e8`): hp_tg/ofct_srec/p_rec のチャンク行数決定に
  使う空きメモリ使用率の既定値。

---

## 6. 起動時間の短縮 — サイズ検証の間引き (`1eb6872`)

`InitReadHiPic` (rhp_c.c) は output.log に載る全画像を1枚ずつ開いて寸法検証しており、
数千枚の入射+透過フレームで **起動〜初表示に10〜20秒** (純粋なファイルオープンI/O)。
GPU/cuFFT初期化は<1秒で無関係と実測確認。

対策: サイズ検証を **100枚ごと(先頭は必ず)** に間引き、オープン回数を約1/100に。
基準サイズ(darkフレーム)取得とNi/Nt計数・ログ記録は従来どおり全数実施。
`rhp_c.c` を使う全プログラム(hp_tg/ofct_srec/ofct_DO等)に効く。
トレードオフ: 均一撮影中の稀な混在サイズ異常が間をすり抜け得る(通常データでは問題なし)。
※ 生データがネットワーク/OneDrive上だと更に遅いので、ローカルSSD配置も推奨。

---

## 7. 検討したが「現状維持」とした事項 (将来の再検討用)

### 7.1 cbp.cu の更なる超高速化 → 実装改良では不可能
BPカーネルだけ無効化して execute を分離計測:

| N, M | execute | FFT等のみ | **逆投影(BP)** | BP占有 |
|---|---|---|---|---|
| 2048,1800 | 16.5 ms | 1.2 ms | 15.3 ms | 93% |
| 4096,1800 | 68.3 ms | 2.4 ms | 65.9 ms | 97% |
| 8000,900 | 130.5 ms | 3.7 ms | 126.8 ms | 97% |

- BPスループット ≈ 360〜390 Gサンプル/s (N非依存=キャッシュ律速)。
- 命令発行ピークの**約50%**、Q読みはL1/L2で捌けDRAM帯域は未飽和。
- **同一アルゴリズム・同一GPUでの実装改良の上限は約2倍**
  (レジスタブロッキング×1.3〜1.6 + ストリーム重畳×1.2 程度)。
- **テクスチャ線形補間の復活は逆効果**(fp32テクスチャ~160G/s < 現行389G/s、精度も8bit)。
  誰かが思いついても戻さないこと。
- 桁を変えるには **NUFFT/gridding 等のアルゴリズム変更**（O(N²M)→O(N²logN)、理論10倍超）
  しかないが、定量CTの画質・CT値忠実度の検証を伴う研究開発になる。
- 手軽な確実解は **GPU更新**(Ada世代で約3倍・変更ゼロ) か **マルチGPUのスライス並列**。
- 結論: **現状維持**。BP_GMFは既にハード限界の約半分で動作している。

### 7.2 読み込みA (libtiff scanline系, 約30本) → 改善余地なし
`TIFFReadScanline`+memcpy を `TIFFReadEncodedStrip` 直読みに置換しても
×0.9〜1.28 (16bit投影ではむしろ悪化)。既に ~2.5GB/s でメモリコピー速度近傍。触らない。

### 7.3 ofct_srec の「N≈8000で異常値」→ ソフトのバグではなかった
ユーザのパラメータ誤入力が原因と判明。調査の過程で以下を**全画素数で健全と確認**:
書き出し(8000²=256MB, 10000²=400MB含む6000+通り)、CBP(CPU/GPU, 全N, 両r0 regime)、
GPUリング除去(sort_filter_g)。共有基盤は大画像でも正しく動作する。

---

## 8. 検証方法（再現用メモ）

- すべて合成ファントム(一様円盤+視野依存リップル)で **best-of-N** 計測。
- 数値検証は `cbp_thread_nai.c`(double) を基準に RMS/max abs diff、および
  「旧実装との bit 一致」「スレッド数/チャンク数非依存性」を要求。
- 書き出しは「libtiffで読み戻して元バッファとバイト一致」を多数の画素数で総当たり。
- 計測ハーネスはリポジトリ外(scratchpad)で実施しており、本リポジトリには含めない。

---

## 9. 未完了・申し送り

- **AVX版CBP のLinux/GCCビルド未検証**。Linux利用時に一度確認すること。
- **配布バイナリ `exe/`** はリポジトリで追跡している。ソース変更時は最終ビルド後に
  `exe/` を再コミットする運用(このセッションでもそうした)。
- 総再構成時間を更に縮めたい場合の次の一手は 7.1 の GPU更新 or NUFFT検討。
