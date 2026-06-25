/*
 * Wavelet_Denoise_3D_Ext.java - 3D Wavelet Denoising ImageJ Plugin
 * (External Process Version - calls tif_wvd_g)
 *
 * Compile: javac -cp ij.jar Wavelet_Denoise_3D_Ext.java
 * Install: Place Wavelet_Denoise_3D_Ext.class and tif_wvd_g in plugins/
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.io.FileSaver;
import ij.plugin.FolderOpener;
import ij.plugin.filter.PlugInFilter;
import ij.process.ImageProcessor;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Path;

public class Wavelet_Denoise_3D_Ext implements PlugInFilter {

    private ImagePlus imp;
    private int levels = 3;
    private double thresholdScale = 1.0;
    private String customExePath = "";
    private boolean keepTempFiles = false;

    private static final String EXE_WIN = "tif_wvd_g.exe";
    private static final String EXE_UNIX = "tif_wvd_g";

    @Override
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        if (imp == null) { IJ.error("Wavelet Denoise 3D Ext", "No image is open"); return DONE; }
        if (imp.getStackSize() < 2) { IJ.error("Wavelet Denoise 3D Ext", "Requires a stack"); return DONE; }
        int bd = imp.getBitDepth();
        if (bd != 8 && bd != 16 && bd != 32) { IJ.error("Wavelet Denoise 3D Ext", "Use 8/16/32-bit images"); return DONE; }
        return DOES_8G | DOES_16 | DOES_32 | STACK_REQUIRED | NO_CHANGES;
    }

    @Override
    public void run(ImageProcessor ip) {
        if (!showDialog()) return;
        long t0 = System.currentTimeMillis();
        IJ.log("=== 3D Wavelet Denoising (External GPU) ===");
        IJ.log("Image: " + imp.getTitle() + " (" + imp.getWidth() + "x" + imp.getHeight() + "x" + imp.getStackSize() + ")");
        IJ.log("Levels: " + levels + ", Threshold scale: " + thresholdScale);
        ImagePlus result = processWithExternalTool();
        IJ.log("Total time: " + ((System.currentTimeMillis() - t0) / 1000.0) + " seconds");
        if (result != null) { result.setTitle(imp.getTitle() + "_WVD3D"); result.show(); }
    }

    private boolean showDialog() {
        GenericDialog gd = new NonBlockingGenericDialog("3D Wavelet Denoising (External GPU)");
        gd.addNumericField("Decomposition levels (1-5):", levels, 0);
        gd.addNumericField("Threshold scale:", thresholdScale, 2);
        gd.addStringField("Custom executable path:", customExePath, 40);
        gd.addCheckbox("Keep temporary files", keepTempFiles);
        gd.addMessage(">1.0 = stronger denoising, <1.0 = weaker\n\n"
                + "Executable: " + (isWindows() ? EXE_WIN : EXE_UNIX));
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        levels = (int) gd.getNextNumber();
        thresholdScale = gd.getNextNumber();
        customExePath = gd.getNextString();
        keepTempFiles = gd.getNextBoolean();
        if (levels < 1 || levels > 5) { IJ.error("Levels must be 1-5"); return false; }
        if (thresholdScale <= 0) { IJ.error("Threshold scale must be positive"); return false; }
        return true;
    }

    private static boolean isWindows() { return System.getProperty("os.name").toLowerCase().contains("win"); }

    private String findExecutable() {
        String exeName = isWindows() ? EXE_WIN : EXE_UNIX;
        if (customExePath != null && !customExePath.trim().isEmpty()) {
            File f = new File(customExePath.trim()); if (f.exists() && f.canExecute()) return f.getAbsolutePath(); }
        String d = IJ.getDirectory("plugins");
        if (d != null) { File f = new File(d, exeName); if (f.exists() && f.canExecute()) return f.getAbsolutePath(); }
        d = IJ.getDirectory("imagej");
        if (d != null) { File f = new File(d, exeName); if (f.exists() && f.canExecute()) return f.getAbsolutePath(); }
        return exeName;
    }

    private File getWorkingDirectory(File tempInputDir) {
        if (imp != null) { ij.io.FileInfo fi = imp.getOriginalFileInfo();
            if (fi != null && fi.directory != null && !fi.directory.isEmpty()) {
                File od = new File(fi.directory); File pd = od.getParentFile();
                if (pd != null && pd.exists()) return pd; if (od.exists()) return od; } }
        File p = tempInputDir.getParentFile(); return (p != null && p.exists()) ? p : tempInputDir;
    }

    private ImagePlus processWithExternalTool() {
        File tid = null, tod = null;
        try {
            Path tb = Files.createTempDirectory("imagej_wvd3d_");
            tid = new File(tb.toFile(), "input"); tod = new File(tb.toFile(), "output");
            if (!tid.mkdirs() || !tod.mkdirs()) { IJ.error("Cannot create temp dirs"); return null; }
            IJ.showStatus("Saving stack..."); if (!saveStack(imp, tid)) { IJ.error("Save failed"); return null; }
            String exe = findExecutable();
            IJ.showStatus("Running Wavelet denoising...");
            int ec = runTool(exe, tid, tod);
            if (ec != 0) { IJ.error("Error code: " + ec); return null; }
            IJ.showStatus("Reading results...");
            return new FolderOpener().openFolder(tod.getAbsolutePath());
        } catch (IOException e) { IJ.error("I/O error: " + e.getMessage()); return null;
        } finally {
            if (!keepTempFiles) { if (tid != null) deleteDir(tid); if (tod != null) deleteDir(tod);
                if (tid != null && tid.getParentFile() != null) tid.getParentFile().delete(); }
            IJ.showStatus(""); IJ.showProgress(1.0);
        }
    }

    private boolean saveStack(ImagePlus imp, File dir) {
        ImageStack s = imp.getStack(); int d = s.getSize();
        int pw = Math.max(4, String.valueOf(d).length()); String fmt = "slice_%0" + pw + "d.tif";
        for (int z = 1; z <= d; z++) {
            if (!new FileSaver(new ImagePlus(String.format(fmt, z), s.getProcessor(z))).saveAsTiff(
                    new File(dir, String.format(fmt, z)).getAbsolutePath())) return false;
            if (z % 50 == 0 || z == d) IJ.showProgress(z, d); } return true;
    }

    private int runTool(String exe, File inDir, File outDir) throws IOException {
        ProcessBuilder pb = new ProcessBuilder(exe, inDir.getAbsolutePath(), outDir.getAbsolutePath(),
                String.valueOf(levels), String.valueOf(thresholdScale));
        pb.redirectErrorStream(true); pb.directory(getWorkingDirectory(inDir));
        IJ.log("Command: " + String.join(" ", pb.command()));
        IJ.log("Working dir: " + pb.directory().getAbsolutePath());
        Process p = pb.start();
        try (BufferedReader r = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
            String line; while ((line = r.readLine()) != null) IJ.log("[ext] " + line); }
        try { return p.waitFor(); } catch (InterruptedException e) { Thread.currentThread().interrupt(); p.destroy(); return -1; }
    }

    private void deleteDir(File d) { if (d == null || !d.exists()) return;
        File[] fs = d.listFiles(); if (fs != null) for (File f : fs) { if (f.isDirectory()) deleteDir(f); else f.delete(); } d.delete(); }
}