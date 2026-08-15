# 純粋な CBP 以外に適用している 3 種類の補正

2026-08-15

本ソフトの再構成は、基本的には平行ビームの畳み込み逆投影(CBP / filtered back
projection)である。すなわち dark と I0 による規格化 → 対数変換 → ランプ系フィルタ
の畳み込み → 逆投影、という手順そのものは教科書どおりである。

これに加えて、標準で 3 種類の補正が入っている。以下にそれぞれの実装箇所、原理、
参考文献をまとめる。

---

## 1. リング(ストライプ)アーチファクト除去 — ソートベースフィルタ(Vo の Algorithm 3)

**実装**: [src/sort_filter_omp.c](src/sort_filter_omp.c)(CPU/OpenMP)、
[src/sort_filter_g.cu](src/sort_filter_g.cu)(GPU)。全再構成プログラムでシノグラム
構築直後に適用。環境変数 `KERNEL_SIZE`(既定 5、0 で無効)。

**原理**: シノグラムを列(検出器チャンネル)ごとに投影角方向で値の大きさ順にソート
し、ソート後の像に対して横方向(チャンネル方向)へメディアンフィルタをかけ、記録して
おいた元位置インデックスで並びを戻す。ソートによって試料由来の構造は角度方向に
単調に並び替えられる一方、特定チャンネルに固定されたゲイン誤差は列全体のオフセット
として残るため、チャンネル方向のメディアンだけがそれを削り、試料信号はほとんど触られ
ない。特に partial stripe(部分ストライプ)に強く、かつ余計なストライプを新たに
生まないのが利点とされている。

**文献**: Vo, Atwood, Drakopoulos, "Superior techniques for eliminating ring artifacts
in X-ray micro-tomography", *Opt. Express* 26(22), 28396-28412 (2018),
doi:10.1364/OE.26.028396 — 本実装は
同論文の Algorithm 3(sorting-based technique)に相当する。従来法との比較には Münch et al.,
*Opt. Express* 17(10), 8567-8591 (2009)(wavelet-FFT 法)が参照になる。

---

## 2. 打ち切り(カッピング)補正 — 投影データの外挿パッド

**実装**: CBP 層に共通 — [src/cbp_thread.c](src/cbp_thread.c) の `DetectPad()` /
`FillPad()`(および `cbp_thread_int.c` / `cbp_thread_nai.c` / `cbp_thread_avx.c`)、
GPU は [src/cbp.cu](src/cbp.cu)。環境変数 `PAD_THRESH`(既定 0.3 = 自動判定 ON、
`PAD_THRESH=0` で強制 OFF)。

**原理**: 試料が視野をはみ出すと投影の両端がゼロに落ちず、ランプフィルタの長距離裾
-1/(2 pi^2 r^2) がその段差を拾って再構成像にカッピング(周辺が持ち上がり中心が沈む)が
出る。本実装は、シノグラム端列の平均振幅が全体平均の `PAD_THRESH` 倍を超えたときだけ
打ち切りありと判定し、フィルタ入力の左右に N/2 画素分(N = 1 投影あたりの横方向画素数
= シノグラムの幅)を付け足してから畳み込む。付け足す値は端の値 p[0] / p[N-1] に
レイズドコサイン窓 w(n) = (1 + cos(pi n / (N/2))) / 2, n = 1…N/2 を掛けたもので、端の
高さから滑らかに 0 へ落とす(ゼロ埋めでは段差が残り、それ自体がカッピングの原因になる)。
逆 FFT 後は先頭 N サンプルだけを取り出すので再構成領域は元の N×N のままであり、視野内に
収まる試料は端が約 0 のため自動的に無適用となる。

**文献**: Ohnesorge, Flohr, Schwarz, Heiken, Bae, "Efficient correction for CT image
artifacts caused by objects extending outside the scan field of view", *Med. Phys.*
27(1), 39-46 (2000)。
また、視野拡張まで踏み込んだものとして Hsieh et al., "A novel reconstruction algorithm
to extend the CT scan field-of-view", *Med. Phys.* 31(9), 2385-2391 (2004)。
放射光ローカル CT 関連では Kyrieleis et al., *Nucl. Instrum. Methods A* 607,
677-684 (2009)。

---

## 3. 低透過率ガード — 対数変換のフロア処理

**実装**: [src/blacklim.h](src/blacklim.h)(`BlackLog()` を各シノグラム構築側が呼ぶ)。
環境変数 `CT_REC_BLACK_THRESH`(既定 2.0 カウント、dark 上)と `CT_REC_BLACK_FRAC`
(既定 0.5)。技術ノート:
[20260806_low_transmission_guard.md](20260806_low_transmission_guard.md)。

**原理**: p = log(I0/(I - dark)) は、厚い/重い試料で I - dark が 0 以下になると +Inf や
NaN を生む。本ガードは画素ごとに信号の下限(フロア)を設け、吸収を log(I0/threshold) で
頭打ちにする。さらに 1 投影のうちクリップ画素の割合が `CT_REC_BLACK_FRAC` を超えた場合は
その投影を「黒」と判定し、良好投影の平均で置換する(欠測角の扱い)。

**文献**: 物理的背景は photon starvation である。Hsieh, "Adaptive streak artifact
reduction in computed tomography resulting from excessive x-ray photon noise",
*Med. Phys.* 25(11), 2139-2147 (1998) が古典的出典。低計数域での対数変換のバイアス
については Whiting et al., "Properties of preprocessed sinogram data in x-ray computed
tomography", *Med. Phys.* 33(9), 3290-3303 (2006)。なお、投影置換の判定則
(`FRAC` によるライン単位判定)は本ソフト固有の実装で、対応する文献はない。
