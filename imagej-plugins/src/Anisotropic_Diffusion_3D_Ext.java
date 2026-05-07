/*
 * Anisotropic_Diffusion_3D_Ext.java - 3D Anisotropic Diffusion ImageJ Plugin
 * (External Process Version - calls tif_adf_g)
 *
 * Compile: javac -cp ij.jar Anisotropic_Diffusion_3D_Ext.java
 * Install: Place Anisotropic_Diffusion_3D_Ext.class and tif_adf_g in plugins/
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
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

public class Anisotropic_Diffusion_3D_Ext implements PlugInFilter {

    private ImagePlus imp;
    private int iterations = 20;
    private double K = -1.0;
    private double dt = 0.14;
    private int mode = 1;
    private String customExePath = "";
    private boolean keepTempFiles = false;

    private static final String EXE_WIN = "tif_adf_g.exe";
    private static final String EXE_UNIX = "tif_adf_g";

    @Override
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        if (imp == null) { IJ.error("Anisotropic Diffusion 3D Ext", "No image is open"); return DONE; }
        if (imp.getStackSize() < 2) { IJ.error("Anisotropic Diffusion 3D Ext", "Requires a stack"); return DONE; }
        int bd = imp.getBitDepth();
        if (bd != 8 && bd != 16 && bd != 32) { IJ.error("Anisotropic Diffusion 3D Ext", "Use 8/16/32-bit images"); return DONE; }
        return DOES_8G | DOES_16 | DOES_32 | STACK_REQUIRED | NO_CHANGES;
    }

    @Override
    public void run(ImageProcessor ip) {
        if (!showDialog()) return;
        long t0 = System.currentTimeMillis();
        IJ.log("=== 3D Anisotropic Diffusion (External GPU) ===");
        IJ.log("Image: " + imp.getTitle() + " (" + imp.getWidth() + "x" + imp.getHeight() + "x" + imp.getStackSize() + ")");
        IJ.log("Iterations: " + iterations + ", K: " + K + ", dt: " + dt + ", mode: " + mode);
        ImagePlus result = processWithExternalTool();
        IJ.log("Total time: " + ((System.currentTimeMillis() - t0) / 1000.0) + " seconds");
        if (result != null) { result.setTitle(imp.getTitle() + "_ADF3D"); result.show(); }
    }

    private boolean showDialog() {
        GenericDialog gd = new GenericDialog("3D Anisotropic Diffusion (External GPU)");
        gd.addNumericField("Iterations:", iterations, 0);
        gd.addNumericField("K (edge threshold, -1=auto):", K, 2);
        gd.addNumericField("dt (time step, <1/6):", dt, 3);
        gd.addChoice("Mode:", new String[]{"1: 1/(1+(s/K)^2)", "2: exp(-(s/K)^2)"}, "1: 1/(1+(s/K)^2)");
        gd.addStringField("Custom executable path:", customExePath, 40);
        gd.addCheckbox("Keep temporary files", keepTempFiles);
        gd.addMessage("K=-1: auto-estimate from noise\n\n"
                + "Executable: " + (isWindows() ? EXE_WIN : EXE_UNIX));
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        iterations = (int) gd.getNextNumber();
        K = gd.getNextNumber();
        dt = gd.getNextNumber();
        String ms = gd.getNextChoice(); mode = ms.startsWith("1") ? 1 : 2;
        customExePath = gd.getNextString();
        keepTempFiles = gd.getNextBoolean();
        if (iterations < 1 || iterations > 10000) { IJ.error("Iterations must be 1-10000"); return false; }
        if (dt <= 0) { IJ.error("dt must be positive"); return false; }
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
            Path tb = Files.createTempDirectory("imagej_adf3d_");
            tid = new File(tb.toFile(), "input"); tod = new File(tb.toFile(), "output");
            if (!tid.mkdirs() || !tod.mkdirs()) { IJ.error("Cannot create temp dirs"); return null; }
            IJ.showStatus("Saving stack..."); if (!saveStack(imp, tid)) { IJ.error("Save failed"); return null; }
            String exe = findExecutable();
            IJ.showStatus("Running Anisotropic Diffusion...");
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
        java.util.List<String> cmd = new java.util.ArrayList<String>();
        cmd.add(exe); cmd.add(inDir.getAbsolutePath()); cmd.add(outDir.getAbsolutePath());
        cmd.add(String.valueOf(iterations));
        if (K > 0) cmd.add(String.valueOf(K)); else cmd.add("-1");
        cmd.add(String.valueOf(dt)); cmd.add(String.valueOf(mode));

        ProcessBuilder pb = new ProcessBuilder(cmd); pb.redirectErrorStream(true);
        pb.directory(getWorkingDirectory(inDir));
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