/*
 * HisImgKifReader.java
 *
 * Fiji (ImageJ2 / SciJava) reader for HIS, IMG, KIF and FPI files.
 *
 * Fiji の HandleExtraFileTypes は IO_ プラグインが使っているため、
 * ImageJ 1.x と同じ方法ではドラッグ＆ドロップに対応できない。
 * 代わりに SciJava の IOPlugin として登録し、Fiji の I/O 機構から呼ばれるようにする。
 *
 * 実際の読み込みは SP8CT_Plugins.jar の Open_HIS_IMG.openFile(String) に任せる
 * （デフォルトパッケージのクラスなのでリフレクションで呼ぶ）。
 * そのため、この jar とは別に SP8CT_Plugins.jar も Fiji の plugins/ に必要。
 *
 * 拡張子だけでは判定せず、ファイル先頭のマジックで確認する。
 * （.img は ANALYZE 形式などと重複するため）
 *   .his / .img : "IM"
 *   .kif        : ".KIF"
 *   .fpi        : 4 + width*height*4 == ファイルサイズ
 */
package jp.spring8.sp8ct.io;

import ij.IJ;
import ij.ImagePlus;

import java.io.File;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.net.URI;
import java.util.Locale;

import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.ImgPlus;
import net.imagej.axis.Axes;
import net.imagej.axis.AxisType;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;

import org.scijava.Priority;
import org.scijava.io.AbstractIOPlugin;
import org.scijava.io.IOPlugin;
import org.scijava.io.location.FileLocation;
import org.scijava.io.location.Location;
import org.scijava.log.LogService;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

@Plugin(type = IOPlugin.class, priority = Priority.HIGH)
public class HisImgKifReader extends AbstractIOPlugin<Dataset> {

    private static final String IJ1_CLASS = "Open_HIS_IMG";

    @Parameter
    private DatasetService datasetService;

    @Parameter(required = false)
    private LogService log;

    @Override
    public Class<Dataset> getDataType() {
        return Dataset.class;
    }

    /* ------------------------------------------------------------
     * 判定と読み込み
     *
     * scijava-common のバージョンによって Location 版と String 版の
     * どちらが呼ばれるかが異なるため、両方用意している。
     * どちらかが親クラスに存在しなくてもコンパイルが通るよう、
     * @Override は付けていない。
     * ------------------------------------------------------------ */

    public boolean supportsOpen(final Location source) {
        return isSupported(toFile(source));
    }

    public boolean supportsOpen(final String source) {
        return isSupported(source == null ? null : new File(source));
    }

    public Dataset open(final Location source) throws IOException {
        return openFile(toFile(source));
    }

    public Dataset open(final String source) throws IOException {
        return openFile(new File(source));
    }

    /* ------------------------------------------------------------ */

    private boolean isSupported(final File f) {
        if (f == null || !f.isFile()) return false;
        final String n = f.getName().toLowerCase(Locale.ROOT);
        try {
            if (n.endsWith(".his") || n.endsWith(".img")) return startsWith(f, "IM");
            if (n.endsWith(".kif")) return startsWith(f, ".KIF");
            if (n.endsWith(".fpi")) return isFpi(f);
        } catch (final IOException e) {
            return false;
        }
        return false;
    }

    private Dataset openFile(final File f) throws IOException {
        if (f == null) throw new IOException("Not a local file");

        final ImagePlus imp = openWithIJ1(f.getAbsolutePath());
        if (imp == null) {
            /* ダイアログでキャンセルされた場合もここに来る */
            throw new IOException("Opening was cancelled or failed: " + f.getName());
        }
        if (imp.getType() == ImagePlus.COLOR_RGB) {
            throw new IOException("RGB KIF is not supported by drag & drop. "
                + "Use Plugins > Open_HIS_IMG instead.");
        }
        return toDataset(imp);
    }

    /** SP8CT_Plugins.jar の Open_HIS_IMG.openFile(String) をリフレクションで呼ぶ */
    private ImagePlus openWithIJ1(final String path) throws IOException {
        try {
            ClassLoader cl = IJ.getClassLoader();
            if (cl == null) cl = getClass().getClassLoader();
            final Class<?> c = Class.forName(IJ1_CLASS, true, cl);
            final Object opener = c.getDeclaredConstructor().newInstance();
            final Method m = c.getMethod("openFile", String.class);
            return (ImagePlus) m.invoke(opener, path);
        } catch (final ClassNotFoundException e) {
            throw new IOException(IJ1_CLASS + " not found. "
                + "Put SP8CT_Plugins.jar into Fiji's plugins folder.", e);
        } catch (final InvocationTargetException e) {
            final Throwable cause = e.getCause();
            throw new IOException("Error in " + IJ1_CLASS + ": " + cause, cause);
        } catch (final ReflectiveOperationException e) {
            throw new IOException("Cannot call " + IJ1_CLASS + ".openFile", e);
        }
    }

    /**
     * ImagePlus を ImgLib2 で包んで Dataset にする。
     * 仮想スタックの場合も、ここで全フレームを読み込むことはない。
     */
    @SuppressWarnings({ "rawtypes", "unchecked" })
    private Dataset toDataset(final ImagePlus imp) {
        final Img img = ImageJFunctions.wrap(imp);
        final AxisType[] axes = (img.numDimensions() == 2)
            ? new AxisType[] { Axes.X, Axes.Y }
            : new AxisType[] { Axes.X, Axes.Y, Axes.Z };
        final ImgPlus imgPlus = new ImgPlus(img, imp.getTitle(), axes);
        return datasetService.create(imgPlus);
    }

    /* ------------------------------------------------------------
     * ヘルパー
     * ------------------------------------------------------------ */

    private static File toFile(final Location source) {
        if (source == null) return null;
        if (source instanceof FileLocation) return ((FileLocation) source).getFile();
        final URI uri = source.getURI();
        if (uri != null && "file".equals(uri.getScheme())) return new File(uri);
        return null;
    }

    private static boolean startsWith(final File f, final String magic) throws IOException {
        final byte[] want = magic.getBytes("US-ASCII");
        if (f.length() < want.length) return false;
        final byte[] got = new byte[want.length];
        try (RandomAccessFile raf = new RandomAccessFile(f, "r")) {
            raf.readFully(got);
        }
        for (int i = 0; i < want.length; i++) {
            if (got[i] != want[i]) return false;
        }
        return true;
    }

    /** FPI: short width, short height (little-endian) + float[width*height] */
    private static boolean isFpi(final File f) throws IOException {
        if (f.length() < 4) return false;
        final byte[] h = new byte[4];
        try (RandomAccessFile raf = new RandomAccessFile(f, "r")) {
            raf.readFully(h);
        }
        final int w = (short) ((h[0] & 0xFF) | ((h[1] & 0xFF) << 8));
        final int ht = (short) ((h[2] & 0xFF) | ((h[3] & 0xFF) << 8));
        if (w <= 0 || ht <= 0) return false;
        return f.length() == 4L + (long) w * ht * 4L;
    }
}
