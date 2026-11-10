# Fiji / ImageJ への SP8CT プラグイン導入手順

2026-11-10

## 概要

ct-rec の ImageJ プラグインを導入すると、ImageJ と Fiji のどちらでも HIS / IMG / KIF / FPI ファイルを開け、フィルター類も使えるようになる。ドラッグ＆ドロップで開く仕組みが ImageJ と Fiji で異なるため、置くファイルが少し違う。

| ファイル | 役割 | ImageJ 1.54 | Fiji |
| --- | --- | --- | --- |
| SP8CT_Plugins.jar | 読み込み本体（Open_HIS_IMG）とフィルター類 | plugins/ | plugins/ |
| HandleExtraFileTypes.class | ImageJ でのドラッグ＆ドロップ対応 | plugins/ | 置かない |
| SP8CT_IO-0.1.0.jar | Fiji でのドラッグ＆ドロップ対応（SciJava IOPlugin） | 置かない | jars/ |

Fiji では jar を置くだけでは足りず、**Edit → Options → ImageJ2** で「Use SCIFIO when opening files」をオンにする必要がある。オフのままだと、ドロップした HIS ファイルが SP8CT_IO に届かず、Bio-Formats Import Options が表示される。それでも Bio-Formats が出る場合は、**Plugins → Bio-Formats → Bio-Formats Plugins Configuration** の Formats タブで「Hamamatsu HIS」の Enabled を外す。

ソースは ct-rec リポジトリの `imagej-plugins/src`（SP8CT_Plugins.jar、HandleExtraFileTypes.class）と `Fiji-plugins/src`（SP8CT_IO: HisImgKifReader.java、pom.xml、build.txt）にある。