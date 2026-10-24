# n2i — Noise2Inverse による CT 像ノイズ除去

Self-supervised denoising of tomographic reconstructions with Noise2Inverse, implemented in PyTorch (2D / 2.5D / 3D U-Net). No clean reference data is needed: the network is trained on reconstructions from disjoint subsets of the projections. Documentation is in Japanese.

Noise2Inverse 法を PyTorch で独自に実装した、CT 再構成像のノイズ除去ツールです。投影データを分割して別々に再構成した像同士から学習するため、ノイズのない正解像は必要ありません。

## 特長

- **2D／2.5D／3D の U-Net**：2.5D は隣接スライスを入力に使い、スライス方向のちらつきを抑えます。
- **アーチファクト対策**：resize-conv によるアップサンプリング、重なり付きタイル推論と cos² ブレンドで、格子模様やタイル境界の継ぎ目を抑えます。
- **大容量データ対応**：学習に使う中央付近のスライスだけを memmap キャッシュにまとめ（RAM・ディスク容量に合わせて自動調整）、推論は TIFF を順に読みます。4k × 4k × 2.3k 級のデータも処理できます。
- **32bit float TIFF の入出力**：出力は元と同じファイル名で、TIFF タグも引き継ぎます。[ct-rec](https://github.com/xrm-bl/ct-rec) 形式の ImageDescription は最小値・最大値を更新するので、ct-rec の後処理（8bit 化など）をそのまま使えます。
- **GPU メモリ不足への自動対処**：学習はバッチを分割し、推論はタイルを縮小して続行します。
- **実行ログ**：コマンド、環境、出力、異常終了時の状況を `n2i.log` に追記します。
- **セットアップスクリプト**：Windows と Ubuntu で、CUDA 対応 PyTorch を含む環境を自動で構築します。

## 動作確認環境

| OS | GPU | 備考 |
| --- | --- | --- |
| Windows 11 | NVIDIA Quadro RTX 8000（48 GB） | 2.5D で学習・推論を確認 |
| Ubuntu <バージョン> | <GPU> | <確認状況> |

NVIDIA ドライバが必要です。CUDA Toolkit は不要です（PyTorch が CUDA ランタイムを同梱しています）。GPU なしでも動きますが、実用的な速度は出ません。

## ファイル構成

| ファイル | 内容 |
| --- | --- |
| [n2i_v2.py](n2i_v2.py) | 本体（学習 `train`／推論 `infer`） |
| [setup_n2i.bat](setup_n2i.bat), [setup_n2i.ps1](setup_n2i.ps1) | Windows 用セットアップ |
| [setup_n2i.sh](setup_n2i.sh) | Ubuntu 用セットアップ |
| [n2i_v2_usage.md](n2i_v2_usage.md) | 詳しい使い方、引数一覧、インストール方法、手法と実装の詳細、引用文献 |
| [n2i_large_volume.md](n2i_large_volume.md) | 大容量データ（4k × 4k × 2.3k 級）の処理 |
| [N2I_install.md](N2I_install.md) | a.his から分割再構成して N2I を行う手順（BL47XU の再構成ソフト `his_rec_g_c`／`his_srec_g_c_n2i` を使う場合） |

## 入力データの準備

1. 投影データを角度インデックスでインターリーブ分割します（i mod K、K は 2 または 4 を推奨）。
2. 各組を同じ条件（回転中心、フィルタ、画素サイズ、スライス範囲）で再構成し、組ごとに別のディレクトリ（例：`rec\0`〜`rec\3`）へ 32bit float TIFF の連番として保存します。
3. ファイル名はゼロ埋めしてください（例：`rec_0001.tif`）。ファイル名順がスライス順として扱われます。

連続した角度ブロックで分割すると、組ごとに限定角度のアーチファクトが異なり、うまく学習できません。

BL47XU の a.his から処理する場合は、分割再構成から推論までの手順を [N2I_install.md](N2I_install.md) にまとめています。

## 注意事項

- 本ツールは Noise2Inverse の公式実装（[ahendriksen/noise2inverse](https://github.com/ahendriksen/noise2inverse)）ではなく、論文に基づく独立した実装です。
- 検出器起因のリングなど、分割間で共通のノイズは除去されません。
- 本ソフトウェアは無保証で提供されます。結果の妥当性は利用者の責任で確認してください。
- 開発にあたり、生成 AI（Anthropic Claude）の支援を受けています。