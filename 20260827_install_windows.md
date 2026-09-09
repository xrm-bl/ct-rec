# Windows PC への導入と最初の再構成まで

2026-08-27

リリース版の zip をダウンロードしてから、実際に1枚再構成するまでの手順。
同梱の Windows x64 バイナリ(`exe/`)をそのまま使う場合を対象とする。ソースから
ビルドする場合は README.md の "Binaries and building" を参照のこと。

---

## 1. 前提条件

OS / ハード: Windows 10 / 11 の 64bit。`exe/` の実行ファイルはすべて x64
ネイティブである。

CPU 版(`_t` で終わるプログラム): Microsoft Visual C++ 2015-2022 再頒布可能
パッケージ(x64)。C ランタイムは静的リンクされているため通常は何も要らないが、
OpenMP を使う 33 本(`ct_rec_t_*`、`hp_tg_t_*` など)だけが `vcomp140.dll` を
必要とし、これが再頒布可能パッケージに含まれる。

GPU 版(`_g` で終わるプログラム): 次の3点が必要。

| 必要なもの | 由来 | 対象 |
|---|---|---|
| NVIDIA ドライバ (`nvcuda.dll`) | ドライバのインストール | GPU 版 28 本すべて |
| `cufft64_12.dll` | CUDA Toolkit (12.x / 13.x) | `ct_rec_g_*`, `hp_tg_g_*`, `ofct_rec_g_*`, `ofct_srec_g_*`, `p_rec_g_*`, `rec2rec_g_*` の 18 本 |
| Turing (sm_75) 以降の GPU | - | 同梱 exe は CUDA Toolkit 13.2 でビルド |

ドライバだけでは足りず、cuFFT を使う 18 本は CUDA Toolkit の導入が必要である。
GPU が無い PC でも CPU 版(`_t`)だけで一通りの処理はできる。

> `ct_rec_tif_g_c` / `_r` / `_s` の3本だけは古い CUDA 11 世代のビルドが残って
> おり、`cudart64_110.dll` と `cufft64_10.dll` を要求する。ver 2.1 で `ct_rec`
> に統合済みの旧プログラムなので、通常は `ct_rec_g_*` を使うこと。

画像表示(任意): 出力は 32bit float TIFF なので、ImageJ または Fiji があると
便利である。同梱の `imagej-plugins/SP8CT_Plugins.jar` を `plugins/` に置くと
HIS/IMG も直接開ける。詳細は imagej-plugins/readme.txt を参照。

## 2. ダウンロードと展開

1. リリースページ https://github.com/xrm-bl/ct-rec/releases/latest から
   `ct-rec_v2.4.zip` を取得する。
2. 展開前に zip のブロックを解除する。zip を右クリック → プロパティ → 全般
   タブ下部の「セキュリティ: 許可する」にチェック → OK。これを忘れると展開した
   exe に「Web からダウンロードされたファイル」の印が残り、SmartScreen に
   止められることがある。
3. 展開する。中身は `ct-rec_v2.4\` という1つのフォルダにまとまっている。

## 3. 設置場所

日本語やスペースを含まないパスに置くのが無難である。バッチファイルが pushd /
popd でディレクトリを渡り歩くため、深いパスや同期フォルダ(OneDrive 等)の下は
避けること。

```
C:\ct-rec\ct-rec_v2.4\exe\      ← 実行ファイルとバッチファイル
C:\ct-rec\ct-rec_v2.4\src\      ← ソース(使わなければ触らなくてよい)
C:\ct-rec\ct-rec_v2.4\readme.txt
```

更新時は新しい zip を別フォルダに展開し、PATH の向き先を差し替えるのが安全で
ある(バージョンを併存させられる)。

## 4. PATH を通す

`exe` フォルダを PATH に追加する。GUI で行うのが確実である。

1. `Win + R` → `rundll32 sysdm.cpl,EditEnvironmentVariables` → Enter
2. 上段の「ユーザー環境変数」の `Path` を選んで「編集」
3. 「新規」で `C:\ct-rec\ct-rec_v2.4\exe` を追加 → OK
4. 開いているコマンドプロンプトは閉じて開き直す(既存プロセスには反映されない)

`setx PATH "%PATH%;C:\ct-rec\ct-rec_v2.4\exe"` でも設定できるが、setx は 1024
文字で切り詰めるため既存の PATH が壊れることがある。GUI を推奨する。

PATH を汚したくない場合は、作業フォルダに次の `env.bat` を置き、作業開始時に
叩く方法もある。

```bat
@echo off
set PATH=C:\ct-rec\ct-rec_v2.4\exe;%PATH%
set KERNEL_SIZE=11
set OMP_NUM_THREADS=16
cmd /k
```

## 5. 動作確認

コマンドプロンプトを開き、引数なしで実行する。使い方が表示されれば導入成功で
ある。

```
> ct_rec_t_c
parameter was wrong!!!
usage : ...\ct_rec_t_c.exe layer (center) (pixel size) (offsetangle)
default pixel size 1.0um
```

続いて GPU 版も確認する。

```
> ct_rec_g_c
```

同じ usage が出れば GPU 版の依存 DLL も揃っている。ここでウィンドウが一瞬で
消える、または「cufft64_12.dll が見つかりません」と出る場合は CUDA Toolkit が
未導入である。`nvidia-smi` でドライバと GPU 世代も確認しておくとよい。

## 6. 環境変数の設定

よく触るものだけ挙げる。全一覧は 20260723_CT_env_vars.md にある。

| 変数 | 既定 | 推奨・注意 |
|---|---|---|
| `KERNEL_SIZE` | 5 | リング除去の強さ。0 または 1 で無効。際立った構造がない試料では 11-17 程度 |
| `OMP_NUM_THREADS` | 40 | リング除去の並列数。未設定だと 40 固定なので、コア数の少ない PC では論理コア数程度に下げる |
| `CBP_THREADS` | 論理コア数-1 | CPU 版逆投影のスレッド数。通常は既定のまま |
| `CUDA_GPU` | 0 | 複数 GPU がある場合の使用番号 |

その場限りなら `set KERNEL_SIZE=11`、恒久設定なら手順4と同じ GUI で「ユーザー
環境変数」に追加する。`set` はそのコマンドプロンプトを閉じると消える。

## 7. 最初の再構成

データの準備: 投影像 `q0001.img …`(または `q0001.tif …`)と暗電流 `dark.img`
(または `dark.tif`)が同じディレクトリにある状態にする。入力形式は dark ファイル
の拡張子で自動判別される。ビームラインから `a.his` + `conv.bat` + `output.log`
の3点で持ち帰った場合は、そのディレクトリで `conv.bat` を実行すると分割・変換
される。

1枚だけ再構成:

```
> cd D:\data\001\raw
> ct_rec_g_c 120
```

`120` は再構成する高さ(レイヤー)である。回転軸位置を省略すると自動推定される。
完了すると `rec00120.tif`(32bit float TIFF)ができる。画素サイズと回転軸位置を
明示する場合は `ct_rec_g_c 120 1024.5 5.64` のように続ける。

確認: `pid rec00120.tif` で TIFF タグに埋め込まれた画素サイズ・回転軸位置・
投影数・最小最大値が表示される。ImageJ で開いてリングやカッピングが無いかも
見ること。

全レイヤーの再構成:

```
> mkdir rec
> hp_tg_g_c raw 5.64 1024.5 0 rec
```

引数は順に、投影像のあるディレクトリ、画素サイズ(um)、回転軸位置、回転軸原点
オフセット、出力ディレクトリである。出力ディレクトリは事前に作成しておく必要が
ある。

8/16bit への変換:

```
> mkdir ro_xy
> tif_f2i 8 rec ro_xy -0.5 3.0
```

## 8. 連続撮影データの一括処理

複数測定をまとめて処理する手順は rec-all_memo.md にある。`rc-check.bat 2101 4`
で試し再構成 → `rc-check\` の画像を目視確認 → `gen-all center.log` で本再構成、
という流れである(自動 CT の場合は `act_rc-check.bat` / `act_gen-all.bat`)。

`gen-all.bat` は画素サイズなどが直書きされているので、exe フォルダのものを直接
編集せず、作業フォルダにコピーしてから編集すること(更新時に上書きされない)。

## 9. うまくいかないとき

| 症状 | 原因と対処 |
|---|---|
| `'ct_rec_g_c' は、内部コマンドまたは外部コマンド…ではありません` | PATH が通っていない、またはコマンドプロンプトを開き直していない |
| GPU 版だけ何も表示されず終了する | `cufft64_12.dll` が無い。CUDA Toolkit を導入する |
| `VCOMP140.DLL が見つかりません` | VC++ 2015-2022 再頒布可能パッケージ(x64)を導入する |
| `CUDA error sort_filter_g.cu:266: an illegal memory access` | ver 2.4 で対策済み。古い exe を使っていないか確認する |
| 低透過率の警告が出る | 試料が厚すぎるか露光不足。20260806_low_transmission_guard.md を参照 |
| 出力ディレクトリが無いと言われる | hp_tg 系は出力先を自動作成しない。mkdir してから実行する |
