### ビルド手順（Windows）

cd imagej-plugins\src
build.bat C:\path\to\ij.jar

これで `SP8CT_Plugins.jar` が生成されます。


### インストール

1. `SP8CT_Plugins.jar` を `plugins/` フォルダに配置して再起動

### D&D で開くための手順

**ImageJ の場合:**  
`HandleExtraFileTypes.class` を `build/` から取り出して `plugins/` フォルダのルートに直接配置
（jar の中ではなく単体の `.class` として）。ImageJ は未知の拡張子に対して自動的にこのクラスを呼び出します。

**Fiji の場合:**  
Fiji は既に独自の `HandleExtraFileTypes` を持っているため、上記の方法は使えません。
代わりに `Plugins > SP8CT > IO > Open HIS/IMG...` からファイルを選択してください。

