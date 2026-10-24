# n2i_v2.py 使い方

Noise2Inverse（X:1 戦略）による CT 像ノイズ除去スクリプト。2D／2.5D／3D U-Net に対応。

## 前提

- 入力は K 個（K≥2、通常は 2 または 4）のサブ再構成フォルダ。各フォルダに 32bit float tiff のスライスを置く。
- 各サブ再構成は、投影を角度インデックスでインターリーブ分割（i mod K）した互いに重ならない組から、同一条件で再構成したもの。
- 全フォルダでスライス枚数・並び順が同じであること。
- ファイル名は文字列としてソートされるため、`rec_0001.tif` のようにゼロ埋めしておくこと。ゼロ埋めされていないと `rec_1, rec_10, rec_100, rec_2…` の順になり、2.5D／3D では隣接スライスの関係が崩れる。
- BL47XU の a.his から分割再構成する手順は [N2I_install.md](N2I_install.md) を参照。

## 必要なリソースの目安

| 項目 | 目安 |
|---|---|
| ディスク（キャッシュ） | K × 学習スライス数 × H × W × 4 バイト。学習スライスは既定で中央付近の最大 500 枚（例：K=2、2048×2048 画像 500 枚で約 17 GB）。ローカル SSD に置く |
| ディスク（出力） | スライス数 × H × W × 4 バイト（入力の 1 分割分。上の例で約 19 GB） |
| GPU メモリ | 2.5D の既定設定で、学習時 7 GB 程度（2048×2048 画像の例）。3D はさらに多く必要。不足時は「トラブル対応」を参照 |
| RAM | 2D／2.5D は数 GB 程度。3D 推論は K × (slab + 2×ov_z) × H × W × 4 バイト（K=4、2048×2048 画像、既定値で約 6.4 GB）。学習のキャッシュは RAM の 70% 以内に自動で制限される |

## 基本の流れ

データのあるディレクトリ（例：`D:\data\sample01`）に移動し、その下の `rec\0`〜`rec\3` にサブ再構成が入っている場合の例。

**1. 学習**

```
C:\n2i\n2i.bat train --mode 25d --splits .\rec\0 .\rec\1 .\rec\2 .\rec\3 --cache .\n2i.npy --out .\n2i_25d.pt --workers 4
```

学習してモデル（`.pt`）を作る。`--cache` を付けると、学習に使うスライス（既定では中央付近の最大 500 枚）だけをまとめたキャッシュを作成し、読み込みが高速になる（GPU の稼働率が上がる）。学習開始時に次のような行が表示されるので、`val=[...]` と `std=...` の値を控えておく（手順 2 で使う）。

```
training slices 326-825 of 1152 (central 500 slices (--max-train-slices))
mode=25d K=4 N=500 2048x2048 train-starts=... val=[246,254) mean=1.589 std=3.944
  (val is indexed within the training range / cache: slices 572-579 of the stack)
```

`val=[...]` は学習範囲（＝キャッシュ）の中での番号。スタック全体での番号は次の行に表示される。

**2. 統計情報の確認（ノイズ分散の下限）**

```
C:\n2i\venv\Scripts\python.exe -c "import numpy as np; z0,z1,s=246,254,3.944; m=np.load(r'.\n2i.npy',mmap_mode='r'); H,W=m.shape[2:]; y,x=(H-512)//2,(W-512)//2; d=(m[0,z0:z1,y:y+512,x:x+512]-m[1,z0:z1,y:y+512,x:x+512])/s; print(d.var()/2)"
```

`z0,z1,s=246,254,3.944` の部分を、手順 1 で控えた値（`val=[z0,z1)` の 2 つの数と `std`）に置き換えて実行する。

- **何を計算しているか**：検証用スライスの中央 512×512 領域（学習時の検証と同じ領域）で、分割 0 と分割 1 の差の分散の 1/2 を求めている。2 つの分割は信号が同じでノイズだけが独立なので、この値は「1 つの分割に含まれるノイズの分散」（正規化した単位）の見積もりになる。
- **なぜ必要か**：学習ログの val 損失は「ネットワークの誤差」と「ターゲット分割自体のノイズ分散」の和なので、0 にはならず、値だけでは良し悪しが判断できない。この値が val 損失の下限にあたる。
- **判断の目安**：学習ログの最良の val からこの値を引いた差が、ネットワークの誤差の目安。差が下限の数 % 程度まで小さければ、学習は十分で、回数を増やしても改善はほとんど見込めない。差が大きければ、学習スライスを増やす（`--train-slices` または `--start`／`--stop`）か、`--iters` を増やす。
- **前提**：キャッシュ（`--cache`）を使っていること、画像が 512×512 以上で `--patch` が既定の 256 であること（検証領域が 512×512 になるため）。分割数 K によらず使える。

**3. 推論**

```
C:\n2i\n2i.bat infer --splits .\rec\0 .\rec\1 .\rec\2 .\rec\3 --model .\n2i_25d.pt --outdir .\denoised --fp16
```

モデルを使ってノイズ除去した tiff を `--outdir` に書き出す。モード（2d／25d／3d）はモデルに記録されているので `--mode` は不要。推論はスライスを先頭から順に TIFF から読むので、キャッシュは不要（学習用のキャッシュは学習範囲しか含まないため、推論では使われない）。

出力は 32bit float tiff で、ファイル名は入力（`--splits` の最初のフォルダ）と同じ。値は学習時の正規化を戻してあり、入力の再構成像と同じ単位（線吸収係数など）になっている。結果は「全投影で再構成した像（K 組の単純平均）」のノイズを除去したものに相当する。

各スライスは、最初のフォルダの同じファイル名のスライスから TIFF タグ（ImageDescription、解像度、Software、独自タグなど）を引き継ぐ。ImageDescription が ct-rec 形式（`画素サイズ\t回転中心\t投影数\t面内回転角\t最小値\t最大値`）の場合は、最小値・最大値だけをノイズ除去後の値に置き換える。そのため、ct-rec の 8bit 化などの後処理をそのまま使える。

## Ubuntu での実行

Ubuntu では次のように読み替える。ランチャーは `~/n2i/n2i`（セットアップ時に `--link` を付けた場合は `n2i` だけでよい）、パスの区切りは `/`。`--workers` は既定の 8 のままでよい。

```bash
# 1. 学習
~/n2i/n2i train --mode 25d --splits ./rec/0 ./rec/1 ./rec/2 ./rec/3 --cache ./n2i.npy --out ./n2i_25d.pt

# 2. 統計情報の確認（z0,z1,s を学習ログの値に置き換える）
~/n2i/venv/bin/python -c "import numpy as np; z0,z1,s=246,254,3.944; m=np.load('./n2i.npy',mmap_mode='r'); H,W=m.shape[2:]; y,x=(H-512)//2,(W-512)//2; d=(m[0,z0:z1,y:y+512,x:x+512]-m[1,z0:z1,y:y+512,x:x+512])/s; print(d.var()/2)"

# 3. 推論
~/n2i/n2i infer --splits ./rec/0 ./rec/1 ./rec/2 ./rec/3 --model ./n2i_25d.pt --outdir ./denoised --fp16
```

以降の節のコマンド例は Windows 形式で書いてあるので、同じように読み替える（`rem` で始まる行はコメント。bash では `#`）。

## 開始と終了

**開始**：学習のコマンド（`train`）を実行する。最初に学習範囲（`training slices ...`）と設定確認の行（`mode=... K=... N=... val=[...]`）が出るので、スライス数・画像サイズが想定どおりか確認する。キャッシュを作る場合は、その進み具合が表示される。その後、`--log` 間隔（既定 500 回）ごとに train／val の損失と、`[gpu 6.8GB  ram_free 210GB  35.2min]`（GPU メモリの最大使用量、空き RAM、経過時間）が表示される。

**終了**：`--iters` 回に達すると自動で終了する。途中で止めるときは Ctrl+C。

- val 損失が最小値を更新するたびに `--out` へ保存される（表示に `*saved`）。中断しても最後に保存されたモデルは使える。
- `_last.pt`（最終状態）は最後まで走った場合のみ作られる。
- 途中から再開する機能はない。中断した場合は最初からやり直し。
- val が下がらなくなったあたりが打ち切りの目安。
- 学習回数を増やしたいときは、`--iters` を大きくして最初から実行し直す（学習率は `--iters` に合わせて減衰するため、続きから学習することはできない）。キャッシュは再利用される。`--out` の名前を変えておくと、前のモデルと比べられる。

## ログ

実行のたびに、カレントディレクトリの `n2i.log` に追記される（`--log-file` で別のファイル、`--no-log` で無効）。記録される内容は次のとおり。

- コマンドライン、作業フォルダ、実行日時
- Python／PyTorch／CUDA のバージョン、GPU 名とメモリ量、RAM 容量
- 既定値を含むすべての引数
- 画面に出た内容すべて（進み具合、エラー表示）
- 終了状態（`=== finished normally`、`=== interrupted by user (Ctrl+C)`、`=== stopped: ...`、`=== terminated by an error`）
- PyTorch・CUDA・ドライバ内部で異常終了した場合は、その時点の状況（`Windows fatal exception` / `Fatal Python error`）

実行ごとに `====…` の区切り線が入る。

## 共通の引数（train／infer）

| 引数 | 内容 |
|---|---|
| `--splits` | サブ再構成のフォルダ（必須、K 個） |
| `--cache` | キャッシュ（.npy）のパス。train では学習範囲のスライスだけを作成・再利用する（3D では強く推奨）。infer では全スライス分の一致するキャッシュがある場合だけ使い、なければ TIFF を読む |
| `--start` `--stop` | 使うスライス範囲。ファイル名順の 0 始まりの通し番号で、`--start` は含み `--stop` は含まない（`--start 100 --stop 200` なら 100〜199 番）。train では学習範囲、infer では推論範囲 |
| `--log-file` | ログファイル（既定：カレントディレクトリの `n2i.log`、追記） |
| `--no-log` | ログを書かない |

## キャッシュの扱い

- **作成範囲**：キャッシュは学習範囲のスライスだけで作られる。学習範囲は、`--start`／`--stop` を指定すればその範囲、`--train-slices` を指定すれば中央付近のその枚数、どちらも指定しなければ中央付近の最大 500 枚（`--max-train-slices`）。ただし、キャッシュが RAM の 70%（`--mem-frac`）とディスクの空きの 90% に収まるように、枚数が自動で減らされる。
- **自動検証**：キャッシュの横に `n2i.npy.json` が作られ、元のファイル（パス・サイズ・更新時刻）と学習範囲が記録される。再構成し直した場合や学習範囲を変えた場合は、古いキャッシュを自動で検出して作り直す。
- **作成時の表示**：TIFF を読み込む進み具合のあと、`writing cache to disk ...` と表示されてから `cache written` が出るまでは、ディスクへの書き出しで時間がかかる（大きいキャッシュでは数分〜数十分）。中断しても作りかけのファイル（`.tmp`）が残るだけで、壊れたキャッシュが使われることはない。
- **推論との関係**：推論は TIFF を順に読むのでキャッシュは不要。学習用のキャッシュは学習範囲しか含まないため、推論では使われない。
- **旧版のキャッシュ**：`.json` がない旧版のキャッシュは、検証できない旨の警告を出したうえで、形が合えば使う。データを変えた場合は削除する。
- **削除**：推論まで終わったら、`n2i.npy` と `n2i.npy.json` を削除してよい。

## train の引数

| 引数 | 既定値 | 内容 |
|---|---|---|
| `--mode` | 2d | `2d`／`25d`／`3d` |
| `--out` | n2i.pt | 保存するモデルのファイル名 |
| `--iters` | 2000 | 学習回数 |
| `--s` | 2 | 2.5D で使う前後のスライス数（2s+1 枚を入力） |
| `--pz` | 32 | 3D パッチの z 方向の厚さ（2^depth の倍数） |
| `--patch` | 256（3D は 128） | 学習パッチの xy サイズ |
| `--batch` | 16（3D は 4） | バッチサイズ |
| `--ppr` | 4（3D は 2） | 1 回の読み込みから切り出すパッチ数 |
| `--depth` | 4（3D は 3） | U-Net の段数 |
| `--ch` | 32 | 第 1 段のチャネル数 |
| `--lr` | 5e-4 | 学習率（cosine で減衰） |
| `--up` | resize | アップサンプリング方式（`transpose` で旧方式） |
| `--nval` | 8 | 検証用に除外するスライス数（範囲の中央から取る） |
| `--workers` | 8 | データ読み込みの並列数（Windows では 4 程度を推奨） |
| `--log` | 500 | 損失表示・保存判定の間隔 |
| `--amp` | なし | bf16 で学習（RTX 30xx 以降） |
| `--train-slices` | なし | 学習に使う中央付近のスライス数（`--start`／`--stop` を指定しない場合） |
| `--max-train-slices` | 500 | 自動で選ぶ学習スライス数の上限。`0` で上限なし |
| `--mem-frac` | 0.7 | キャッシュの大きさを RAM のこの割合以内に制限する |

## infer の引数

| 引数 | 既定値 | 内容 |
|---|---|---|
| `--model` | （必須） | 学習したモデル |
| `--outdir` | （必須） | 出力フォルダ。元と同じファイル名で書き出す |
| `--tile` | 2D: `1024 1024`／3D: `64 256 256` | 推論タイルの大きさ（2D は y x、3D は z y x） |
| `--ov` | 2D: `64 64`／3D: `16 32 32` | タイルの重なり幅 |
| `--slab` | 64 | 3D で一度に処理するスライス数 |
| `--fp16` | なし | 半精度で推論（高速・省メモリ） |

## よく使う例

まず一部のスライスで短く試し、結果を確認してから全体を処理する。

```
rem 短い試験学習（1000 回、スライス 200〜599 だけ使用。--train-slices 400 なら中央 400 枚）
C:\n2i\n2i.bat train --mode 25d --splits .\rec\0 .\rec\1 .\rec\2 .\rec\3 --cache .\n2i.npy --out .\test.pt --iters 1000 --start 200 --stop 600 --workers 4

rem 一部スライスだけ推論して確認
C:\n2i\n2i.bat infer --splits .\rec\0 .\rec\1 .\rec\2 .\rec\3 --model .\test.pt --outdir .\test_out --start 400 --stop 410 --fp16
```

3D の場合：

```
C:\n2i\n2i.bat train --mode 3d --splits .\rec\0 .\rec\1 .\rec\2 .\rec\3 --cache D:\cache\n2i.npy --out .\n2i_3d.pt --amp --workers 4
C:\n2i\n2i.bat infer --splits .\rec\0 .\rec\1 .\rec\2 .\rec\3 --model .\n2i_3d.pt --outdir .\denoised_3d --fp16
```

## トラブル対応

**GPU メモリ不足（`CUDA out of memory`）**

自動で次の対処が行われ、画面とログにその旨が表示される。

- 学習時：バッチを分割して計算する（勾配を足し合わせるので、学習としては同じバッチサイズのまま）。パッチ 1 枚でも足りない場合は停止するので、`--patch` を小さくする（3D なら `--pz 16 --patch 96 --ch 24` など）。
- 推論時：タイルを半分にして続ける。

手動で設定する場合は、学習は `--batch`／`--patch`、推論は `--tile` を小さくする。いずれも `--ov` は残す（重なりを減らすとタイル境界が目立ちやすくなる）。

**エラー表示なしで止まる**

1. `n2i.log` の最後を見る。`=== finished normally` なら正常終了。Traceback と `=== terminated by an error` があれば Python のエラー。`Windows fatal exception` / `Fatal Python error` があれば PyTorch・CUDA・ドライバ内部での異常終了。進み具合の行で途切れていれば外部からの強制終了。
2. 止まった直後に `echo %ERRORLEVEL%` を実行する。`-1073741819` はアクセス違反（ドライバや CUDA 内部の異常）、`-1073741510` は Ctrl+C またはウィンドウを閉じた、`1` は Python のエラー。
3. Windows のイベントビューアー（Windows ログ）を見る。「アプリケーション」に `python.exe` の Application Error があれば、異常終了したモジュールがわかる。「システム」にソース Display、イベント ID 4101（ディスプレイ ドライバー nvlddmkm の応答停止）があれば TDR が原因。

**TDR（GPU のタイムアウト）**：画面表示も担当している GPU（WDDM モード）では、1 回の GPU 処理が約 2 秒を超えると Windows がドライバをリセットし、プログラムが表示なしで終了することがある。推論の `--tile` を小さくする（`2048 2048` や既定の 1024）。それでも起きる場合は、レジストリ `HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\GraphicsDrivers` に DWORD 値 `TdrDelay`（例：60 秒）を追加して再起動する（管理者権限が必要。所内の管理ルールに従うこと）。

**3D 推論の RAM 使用量**

おおよそ K × (slab + 2×ov_z) × H × W × 4 バイト。2048² スライスで既定値なら約 6.4 GB。足りなければ `--slab` を小さくする。

**読み込みが遅い**

- 学習では `--cache` を使い、ローカル SSD に置く（ネットワークドライブ不可）。
- Windows Defender のリアルタイム保護からデータフォルダとキャッシュを除外する。

**GPU の稼働率が低い**

- 稼働率は別のコマンドプロンプト（端末）で `nvidia-smi -l 1` を実行し、GPU-Util で確認する。Windows のタスクマネージャーの GPU グラフは、既定では CUDA の計算を表示しない。
- 低い原因の多くはデータの読み込み待ち。`--cache` を使う。キャッシュを使えない場合は `--ppr` を増やす（例：`--batch 32 --ppr 16`）。
- 読み込みが解消してから `--batch` を大きくする。

**Windows で起動直後に止まる・エラーになる**

- `--workers 0` で実行して切り分ける。

## 結果の確認

- K 組の単純平均（≒フル再構成）とノイズ除去結果を並べて比較する。
- 差分（平均 − ノイズ除去後）に構造が見える場合は細部が削られている可能性がある。学習回数を減らす、学習スライスを増やす、`--ch` を下げる、などを試す。
- 検出器起因のリングなど、分割間で共通のノイズは除去されない。
- 2.5D／3D では、ボリュームの最上端・最下端の数枚は隣接スライスが足りないため、端のスライスを繰り返して代用している。この数枚はわずかに質が落ちる可能性がある。
- 学習範囲から離れたスライス（上端・中央・下端など）も、差分画像で確認しておくとよい。

### 検証損失の評価（ノイズ分散の下限）

「基本の流れ」の手順 2 を参照。

## モデルの使い回し

- **同じデータの別のスライス**：試料の構造がスライス方向でほぼ同じなら、一部のスライス（既定の中央 500 枚、`--train-slices`、`--start`／`--stop`）で学習したモデルで全スライスを推論してよい。学習範囲から離れたスライスで差分画像を確認すること。構造が大きく変わる場合（試料の端、容器の壁際など）は、全範囲で学習し直す。
- **別の測定データ**：撮影条件（エネルギー、露光時間、投影数）、再構成の条件、分割数が同じで、試料も似ていれば使える可能性がある。ただしモデルには学習データの平均・標準偏差が記録されていて推論時もその値で正規化するので、値の範囲が大きく違うデータには使えない。ノイズの強さが違う場合も学習し直す。Noise2Inverse は対象データ自体から学習する手法なので、データごとに学習するのが基本。短い学習（スライス範囲を絞って 1000〜2000 回）でも十分な場合が多い。

---

# 付帯情報

## ファイル構成

| ファイル | 内容 |
|---|---|
| n2i_v2.py | 本体。学習（train）と推論（infer） |
| setup_n2i.bat / setup_n2i.ps1 | Windows 用セットアップ |
| setup_n2i.sh | Ubuntu 用セットアップ |
| n2i_v2_usage.md | 本書 |
| n2i_large_volume.md | 大容量データ（4k × 4k × 2.3k 級）の処理 |
| N2I_install.md | a.his から分割再構成して N2I を行う手順（BL47XU の再構成ソフト his_rec_g_c／his_srec_g_c_n2i を使う場合） |
| README.md | 概要 |

必要なライブラリは PyTorch、NumPy、tifffile のみ。Noise2Inverse の手法を独自に実装したもので、公式実装（[ahendriksen/noise2inverse](https://github.com/ahendriksen/noise2inverse)）や、それが依存する msd_pytorch は使っていない。論文・報告書では「Noise2Inverse 法を独自に実装して適用した」と記述し、手法の原著論文を引用する。

## インストール

`setup_n2i.*` と `n2i_v2.py` を同じフォルダに置いてセットアップスクリプトを実行する。スクリプトは次の処理を自動で行う。

1. `nvidia-smi` で NVIDIA ドライバと、対応する CUDA の上限バージョンを確認する。
2. 64bit の Python（3.10 以上）を探す。
3. インストール先に仮想環境（`venv`）を作る。
4. ドライバが対応する範囲で最も新しい CUDA ビルドの PyTorch を入れ、GPU で畳み込みを実際に実行して確認する。失敗した場合は一つ古いビルドに切り替える。
5. NumPy と tifffile を入れ、`n2i_v2.py` をコピーし、起動用ランチャーを作る。

CUDA Toolkit のインストールは不要（PyTorch が CUDA ランタイムを同梱しているため、NVIDIA ドライバだけあればよい）。PyTorch 本体は 2〜3 GB あるので初回は時間がかかる。再実行した場合、正常に動く PyTorch が既に入っていればそのまま残す。

### Windows

`setup_n2i.bat` をダブルクリックするか、コマンドプロンプトで実行する。既定のインストール先は `C:\n2i`。

```
setup_n2i.bat
setup_n2i.bat -InstallDir D:\n2i -InstallPython
```

| オプション | 内容 |
|---|---|
| `-InstallDir <パス>` | インストール先（既定 `C:\n2i`。PyTorch はパスが深いので短いパスにする） |
| `-CudaTag cu126` | PyTorch の CUDA ビルドを固定 |
| `-InstallPython` | Python が無ければ winget で Python 3.12 を入れる |
| `-Cpu` | GPU なしの PC 用（CPU 版 PyTorch） |
| `-Reinstall` | PyTorch を入れ直す |

実行は `C:\n2i\n2i.bat train ...` のようにランチャー経由で行う（仮想環境のアクティベートや実行ポリシーの変更は不要）。

### Ubuntu

```bash
bash setup_n2i.sh
bash setup_n2i.sh --install-deps --link
```

既定のインストール先は `~/n2i`。

| オプション | 内容 |
|---|---|
| `--dir <パス>` | インストール先（既定 `~/n2i`） |
| `--cuda-tag cu126` | PyTorch の CUDA ビルドを固定 |
| `--install-deps` | `python3-venv` などを apt で入れる（sudo を使用） |
| `--cpu` | GPU なしのマシン用 |
| `--reinstall` | PyTorch を入れ直す |
| `--link` | `~/.local/bin/n2i` にリンクを作り、`n2i train ...` で実行できるようにする |
| `--python <コマンド>` | 使う Python を指定（例 `python3.12`） |

実行は `~/n2i/n2i train ...` のようにランチャー経由で行う。

### インストール時の注意

- **ファイル名**：`n2i_v2.py` という名前でセットアップスクリプトと同じフォルダに置く。ダウンロード時に名前が変わっている（`n2iv2.py`、`n2i_v2 (1).py`、`n2i_v2.py.txt` など）とコピーされない。Windows ではエクスプローラーで拡張子を表示する設定にしておくとよい。
- **ドライバ**：`nvidia-smi` が動かない場合は NVIDIA ドライバを入れる（Ubuntu では `sudo ubuntu-drivers install` の後に再起動）。
- **Python（Ubuntu）**：`python3-venv` が入っていないと仮想環境の作成で止まる。`--install-deps` を付けるか `sudo apt install python3-venv` を実行する。
- **プロキシ**：プロキシ経由の接続が必要な環境では、実行前に設定する（Windows：`set HTTPS_PROXY=http://proxy:port`、Ubuntu：`export https_proxy=http://proxy:port`）。
- **改行コード（Ubuntu）**：`setup_n2i.sh` を Windows で編集・保存すると CRLF になり、`$'\r': command not found` で止まる。`sed -i 's/\r$//' setup_n2i.sh` で直す。
- **CUDA ビルドの一覧**：試す CUDA ビルドはスクリプト冒頭の `$TorchTags`（Windows）／`TORCH_TAGS`（Ubuntu）に書いてある。PyTorch 側で新しいビルドが追加・削除された場合はこの行を編集する（存在しないビルドは自動でスキップされる）。
- **bf16（`--amp`）**：セットアップの最後に、GPU が bf16 に対応しているかが表示される。対応しない GPU（Turing 世代以前）では `--amp` を使わない。

## 手法

投影データを角度インデックスでインターリーブ分割（i mod K）し、それぞれを独立に再構成した K 個のサブ再構成を入力とする。

- **学習（X:1 戦略）**：分割 k をターゲット、残り K−1 個の平均を入力とし、k をランダムに選んで L2 損失で学習する。入力とターゲットのノイズが独立なので、ネットワークはノイズのない像の条件付き期待値を推定するように学習される。
- **推論**：K 通りの「入力 → 出力」を計算し、その平均を最終結果とする。

## ネットワーク

| 項目 | 内容 |
|---|---|
| 構造 | U-Net（2D／3D 畳み込みを切り替え） |
| 段数・チャネル | 2D／2.5D：4 段、3D：3 段。第 1 段 32ch で段ごとに倍増 |
| 畳み込みブロック | 3×3（3D は 3×3×3）畳み込み ×2、各々 BatchNorm と LeakyReLU（傾き 0.1）を付加。パディングは replicate |
| ダウンサンプリング | MaxPool（2 倍） |
| アップサンプリング | 双線形／三線形補間で 2 倍にした後 3×3 畳み込み（resize-conv）。転置畳み込みによる格子状アーチファクトを避けるため。`--up transpose` で転置畳み込みも選択可 |
| 出力 | 1×1 畳み込みに入力中央スライスを加える残差構成 |
| 2.5D | 隣接 2s+1 枚をチャネル入力し、中央 1 枚を出力 |

## 学習（既定値）

| 項目 | 2D／2.5D | 3D |
|---|---|---|
| 損失関数 | 平均二乗誤差（L2） | 同左 |
| 最適化 | Adam、学習率 5×10⁻⁴、cosine annealing で減衰 | 同左 |
| 学習回数 | 2000 回 | 同左 |
| パッチ | 256×256、バッチ 16（1 回の読み込みから 4 パッチ） | 32×128×128、バッチ 4（1 回の読み込みから 2 パッチ） |
| データ拡張 | 90° 回転、左右反転、スライス方向の反転 | 同左 |
| 正規化 | 学習データから求めた平均・標準偏差で標準化（モデルに保存され、推論でも使用） | 同左 |
| 学習範囲 | `--cache` 使用時、`--start`／`--stop` を指定しなければ中央付近の最大 500 枚（RAM の 70%・ディスク空きの範囲内） | 同左 |
| 検証 | 学習範囲の中央 `--nval` 枚（既定 8 枚）を学習から除外 | 同左 |
| 精度 | float32（`--amp` で bf16） | 同左 |

## 推論

画像をタイル（2D：1024×1024、重なり 64 画素／3D：64×256×256、重なり 16×32×32）に分け、重なり部分を cos² 重みでブレンドしてつなぐ。3D は z 方向にマージンを付けたスラブ単位で処理する。BatchNorm は推論時に固定統計を使うため、出力はタイルの切り方に依存しない。`--fp16` で半精度推論。GPU メモリが不足するとタイルを自動で半分にする。出力の TIFF タグは入力スライスから引き継ぐ（ct-rec 形式の ImageDescription は最小値・最大値を更新）。

## 引用文献

**必ず引用するもの**

1. A. A. Hendriksen, D. M. Pelt, K. J. Batenburg, "Noise2Inverse: Self-Supervised Deep Convolutional Denoising for Tomography," *IEEE Transactions on Computational Imaging* **6**, 1320–1335 (2020). doi:10.1109/TCI.2020.3019647（arXiv:2001.11801）
   - 手法の原著論文。
2. A. A. Hendriksen, M. Bührer, L. Leone, M. Merlini, N. Viganò, D. M. Pelt, F. Marone, M. di Michiel, K. J. Batenburg, "Deep denoising for multi-dimensional synchrotron X-ray tomography without high-quality reference data," *Scientific Reports* **11**, 11895 (2021). doi:10.1038/s41598-021-91084-8
   - 放射光 CT への適用と、空間・時間・スペクトル方向への拡張。2.5D／3D を使う場合に特に関係する。

**手法の背景・構成要素として引用を推奨するもの**

3. J. Lehtinen et al., "Noise2Noise: Learning Image Restoration without Clean Data," *Proc. ICML 2018*, PMLR **80**, 2965–2974 (2018).
4. O. Ronneberger, P. Fischer, T. Brox, "U-Net: Convolutional Networks for Biomedical Image Segmentation," *MICCAI 2015*, LNCS **9351**, 234–241 (2015). doi:10.1007/978-3-319-24574-4_28
5. Ö. Çiçek et al., "3D U-Net: Learning Dense Volumetric Segmentation from Sparse Annotation," *MICCAI 2016*, LNCS **9901**, 424–432 (2016). doi:10.1007/978-3-319-46723-8_49（3D モード使用時）
6. A. Paszke et al., "PyTorch: An Imperative Style, High-Performance Deep Learning Library," *Advances in Neural Information Processing Systems* **32** (NeurIPS 2019).

**実装の細部まで記述する場合に引用するもの**

7. A. Odena, V. Dumoulin, C. Olah, "Deconvolution and Checkerboard Artifacts," *Distill* (2016). doi:10.23915/distill.00003（resize-conv の根拠）
8. K. Zhang et al., "Beyond a Gaussian Denoiser: Residual Learning of Deep CNN for Image Denoising," *IEEE Trans. Image Process.* **26**, 3142–3155 (2017).（残差学習）
9. D. P. Kingma, J. Ba, "Adam: A Method for Stochastic Optimization," *ICLR 2015* (arXiv:1412.6980).
10. S. Ioffe, C. Szegedy, "Batch Normalization: Accelerating Deep Network Training by Reducing Internal Covariate Shift," *Proc. ICML 2015*, PMLR **37**, 448–456 (2015).
11. I. Loshchilov, F. Hutter, "SGDR: Stochastic Gradient Descent with Warm Restarts," *ICLR 2017*.（cosine annealing）
12. C. R. Harris et al., "Array programming with NumPy," *Nature* **585**, 357–362 (2020).
13. C. Gohlke, tifffile（ソフトウェア。Zenodo の DOI と使用バージョンを記載）

引用の際は、各出版元で最新の書誌情報を確認すること。

## 方法の記述例（英文テンプレート）

`< >` 内を実際の条件に置き換えて使う。

> Reconstructed volumes were denoised using Noise2Inverse [1,2], implemented in-house in PyTorch [6]. The <N> projections were split into <K> interleaved subsets by angle index, which were reconstructed independently. A <2.5D> U-Net [4] (<five adjacent slices as input channels>, bilinear upsampling followed by convolution [7], residual output) was trained with the X:1 strategy using an L2 loss and the Adam optimizer [9] (learning rate 5×10⁻⁴, cosine annealing [11], <iterations> iterations, <patch> patches, batch size <batch>). Training used slices <range>, with <n> slices held out for validation. The final image was obtained by averaging the network outputs over all <K> splitting strategies.