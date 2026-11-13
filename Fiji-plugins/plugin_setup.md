# Fiji / ImageJ への SP8CT プラグイン導入手順

2026-11-10(2026-11-13 Fast_Image_Sequence_Hook.js を追記)

## 概要

ct-rec の ImageJ プラグインを導入すると、ImageJ と Fiji のどちらでも HIS / IMG / KIF / FPI ファイルを開け、フィルター類も使えるようになる。ドラッグ＆ドロップで開く仕組みが ImageJ と Fiji で異なるため、置くファイルが少し違う。

| ファイル | 役割 | ImageJ 1.54 | Fiji |
| --- | --- | --- | --- |
| SP8CT_Plugins.jar | 読み込み本体（Open_HIS_IMG）とフィルター類 | plugins/ | plugins/ |
| HandleExtraFileTypes.class | ImageJ でのドラッグ＆ドロップ対応 | plugins/ | 置かない |
| SP8CT_IO-0.1.0.jar | Fiji でのドラッグ＆ドロップ対応（SciJava IOPlugin） | 置かない | jars/ |
| Fast_Image_Sequence_Hook.js | Fiji で連番画像の読み込みを ImageJ1 のリーダーに戻す起動スクリプト | 置かない | scripts/Plugins/AutoRun/ |

Fiji では jar を置くだけでは足りず、**Edit → Options → ImageJ2** で「Use SCIFIO when opening files」をオンにする必要がある。オフのままだと、ドロップした HIS ファイルが SP8CT_IO に届かず、Bio-Formats Import Options が表示される。それでも Bio-Formats が出る場合は、**Plugins → Bio-Formats → Bio-Formats Plugins Configuration** の Formats タブで「Hamamatsu HIS」の Enabled を外す。

## 連番画像の読み込み(Fast_Image_Sequence_Hook.js)

「Use SCIFIO when opening files」をオンにすると、Fiji は連番画像も 1 枚ずつ SCIFIO 経由で読むようになり、読み込みが遅くなる。これを避けるため、`Fiji-plugins/scripts/Plugins/AutoRun/Fast_Image_Sequence_Hook.js` を Fiji の `scripts/Plugins/AutoRun/` に置く(Windows なら `Fiji.app\scripts\Plugins\AutoRun\`。AutoRun フォルダーが無ければ作る)。このフォルダーのスクリプトは Fiji の起動時に自動実行されるので、置いたあと Fiji を再起動する。

スクリプトが効くのは次の 3 つで、いずれも ImageJ1 本来のリーダーで読むようになる。

- File → Import → Image Sequence...
- フォルダーのドラッグ＆ドロップ
- マクロの `File.openSequence` と `run("Image Sequence...")`

単独のファイルを開く場合は何も変えないので、HIS / IMG / KIF / FPI のドラッグ＆ドロップは従来どおり SP8CT_IO に届く。外すときはこのファイルを削除して Fiji を再起動する。組み込みに失敗した場合は Log ウィンドウに `Fast_Image_Sequence_Hook: install failed` と表示される。

ソースは ct-rec リポジトリの `imagej-plugins/src`（SP8CT_Plugins.jar、HandleExtraFileTypes.class）と `Fiji-plugins/src`（SP8CT_IO: HisImgKifReader.java、pom.xml、build.txt）にある。Fast_Image_Sequence_Hook.js は `Fiji-plugins/scripts/Plugins/AutoRun/` にある。