# libtiff 3.6.0 → 4.6.0 移行メモ

CT 再構成スイート（ct-rec）が使う libtiff を **3.6.0 → 4.6.0** に更新した際の記録。
クラシック TIFF（8/16/32bit）の読み書き互換は完全に保たれる。

## 何を変えたか

### 1. ソース（型の標準化）
libtiff 4.4.0 で旧整数 typedef（`uint8`/`uint16`/`uint32`）が削除され、
標準の `<stdint.h>` 型（`uint8_t`/`uint16_t`/`uint32_t`）へ統一された。
これに未対応だった 4 ファイルを修正（幅は同一なので挙動は不変）:

- `tif_blf.c`, `tif_gsf.c`, `tif_mdf.c`, `tif_mgf.c`
  - `#include <stdint.h>` を追加
  - `uint32`→`uint32_t`, `uint16`→`uint16_t`（コメント含む）

他のファイルは既に `*_t` + `<stdint.h>` へ移行済みだった。CUDA 版は `MIGRATION.md` 参照。

### 2. ヘッダ / ライブラリ（Windows）
`src/` 直下の以下を 4.6.0 のものへ差し替え:

| ファイル | 内容 |
|---|---|
| `tiff.h`, `tiffio.h`, `tiffvers.h`, `tiffconf.h` | 4.6.0 の公開ヘッダ（4.x は `tiffconf.h` も必要） |
| `libtiff.lib` | 4.6.0 静的ライブラリ（vcpkg の `tiff.lib` をリネーム） |
| `jpeg.lib`, `lzma.lib`, `zs.lib` | 依存の静的ライブラリ（libjpeg-turbo / liblzma / zlib） |

旧 3.6.0 一式は `src/libtiff-3.6.0/` に退避（ロールバック用）。

### 3. ビルドスクリプト（リンク行）
静的 tiff は依存（jpeg/lzma/zlib）も一緒にリンクする必要がある。
各 Windows バッチに `TIFFLIB` マクロを追加し、`libtiff.lib` の記述を置換:

```bat
set TIFFLIB=libtiff.lib jpeg.lib lzma.lib zs.lib
```

対象: `makefileCPU.bat`, `makefileGPU.bat`, `cccc.bat`
（コンパイルフラグ `/MD`・`/MT` の変更は不要。下記参照）

Linux 側（`MakefileCPU` / `MakefileGPU`）は `-ltiff`（システムの共有 libtiff、
多くの環境で既に 4.x）が依存を自動解決するため変更不要。ソース修正の恩恵のみ受ける。

## CRT モデルについて（重要）

- 既存ビルドは `cl`（フラグ無指定）の既定 **/MT（静的CRT）**。nvcc の既定ホスト CRT も /MT。
- そのため vcpkg トリプレットは **`x64-windows-static`（/MT）** を採用。
  - 旧 3.6.0 lib も /MT だったため、コンパイルフラグを一切変えずリンクできる。
  - 生成 exe はスタンドアロン（動的 CRT DLL 不要）。
- `x64-windows-static-md`（/MD）を使うと `__imp_realloc` 等が未解決になり、
  全コンパイルに `/MD` 追加が必要になる。互換性維持のため /MT を選択した。

## libtiff 4.6.0 の入手（vcpkg・再現手順）

**注意**: プロジェクトのパスに日本語・空白が含まれると vcpkg のビルドが失敗する。
必ず **ASCII パス**で vcpkg を運用し、成果物だけを `src/` にコピーする。

```powershell
# 1) ASCII パスに vcpkg を用意（浅いクローン不可: 過去バージョン取得に全履歴が要る）
git clone https://github.com/microsoft/vcpkg C:\Users\<user>\vcpkg
C:\Users\<user>\vcpkg\bootstrap-vcpkg.bat -disableMetrics

# 2) ASCII の作業ディレクトリに vcpkg.json（4.6.0 固定）を置いて install
$env:VCPKG_ROOT = "C:\Users\<user>\vcpkg"
cd C:\Users\<user>\ct-rec-tiff-build       # vcpkg.json を配置
C:\Users\<user>\vcpkg\vcpkg.exe install --triplet x64-windows-static

# 3) 成果物を src/ へコピー
#   include\{tiff.h,tiffio.h,tiffvers.h,tiffconf.h} -> src\
#   lib\tiff.lib -> src\libtiff.lib   （リネーム）
#   lib\{jpeg.lib,lzma.lib,zs.lib}    -> src\
```

`vcpkg.json`（リポジトリ直下）でバージョンを 4.6.0 に固定している:
```json
{ "overrides": [ { "name": "tiff", "version": "4.6.0", "port-version": 5 } ] }
```
別バージョン（例 4.7.x）に上げたい場合は overrides を変更して再 install → 再コピー。

## 検証結果

- 代表 CPU プログラム（`ct_rec_t_r`, `tif_mgf`, `tif_blf`, `ct_prj_f`, `hp2DO`）が
  既定 /MT + 新 libtiff でリンク成功。
- 16bit TIFF の生成→`tif_mgf`（median+gaussian）→読み戻しの round-trip が正常動作。
  出力の幾何・画素値とも妥当（`TIFFGetVersion()` = 4.6.0）。
- GPU（nvcc）版 `ofct_DO_g.exe` も新 libtiff でリンド確認。

## 既知の注意点

- `tif_mgf` は入力 TIFF に ImageDescription タグが無いと、出力時に
  `TIFFSetField(IMAGEDESCRIPTION, NULL)` で落ちる（libtiff 4.x は NULL 文字列で crash）。
  実データは当該タグを持つため通常問題にならないが、潜在バグとして記録。
