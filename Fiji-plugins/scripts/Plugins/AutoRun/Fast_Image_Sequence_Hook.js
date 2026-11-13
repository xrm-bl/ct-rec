// Fast_Image_Sequence_Hook.js  (Fiji startup script)
//
// "Use SCIFIO when opening files" をオンにしたまま、
//   - File > Import > Image Sequence...
//   - フォルダのドラッグ&ドロップ
//   - マクロの File.openSequence / run("Image Sequence...")
// だけを ImageJ1 本来のリーダーで読ませる。
//
// 置き場所: Fiji.app/scripts/Plugins/AutoRun/  (Fiji 起動時に自動実行される)
// 外し方  : このファイルを削除して Fiji を再起動
//
// 仕組み: Fiji は連番の1枚ごとに LegacyOpener プラグイン経由で SCIFIO を呼ぶ。
//         その列の先頭に、呼び出し元が ij.plugin.FolderOpener のときだけ
//         ImageJ1 の Opener で直接読んで返すオープナーを差し込む。
//         単独ファイルを開く場合は何もしない (従来どおり SCIFIO)。

importClass(Packages.ij.IJ)
importClass(Packages.ij.io.Opener)
importClass(Packages.java.io.File)
importClass(Packages.java.lang.System)
importClass(Packages.java.lang.Thread)
importClass(Packages.net.imagej.legacy.plugin.LegacyOpener)

function installFastSequenceOpener() {
	var KEY = "fast.sequence.opener";

	var fast = new LegacyOpener({
		open: function(path, planeIndex, displayResult) {
			if (displayResult || planeIndex != -1) return null;
			var st = Thread.currentThread().getStackTrace();
			for (var i = 0; i < st.length; i++) {
				if (st[i].getClassName().equals("ij.plugin.FolderOpener")) {
					var f = new File(path);
					var o = new Opener();
					o.setSilentMode(true);
					return o.openTempImage(f.getParent() + File.separator, f.getName());
				}
			}
			return null;
		}
	});

	var hf = IJ.class.getDeclaredField("_hooks");
	hf.setAccessible(true);
	var hooks = hf.get(null);
	var lf = hooks.getClass().getDeclaredField("legacyOpeners");
	lf.setAccessible(true);
	var openers = lf.get(hooks);

	var old = System.getProperties().get(KEY);  // 二重登録を避ける
	if (old != null) openers.remove(old);
	openers.add(0, fast);
	System.getProperties().put(KEY, fast);
}

try {
	installFastSequenceOpener();
} catch (e) {
	IJ.log("Fast_Image_Sequence_Hook: install failed: " + e);
}

null
