/*
 * NLM_Filter_3D_Ext.java - 3D Non-Local Means Filter ImageJ Plugin
 * (External Process Version - calls tif_nlm_g)
 *
 * Compile: javac -cp ij.jar NLM_Filter_3D_Ext.java
 * Install: Place NLM_Filter_3D_Ext.class and tif_nlm_g in plugins/
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

public class NLM_Filter_3D_Ext implements PlugInFilter {

    private ImagePlus imp;
    private int patchRadius = 1;
    private int searchRadius = 3;
    private double h = -1.0;
    private String customExePath = "";
    private boolean keepTempFiles = false;

    private static final String EXE_WIN = "tif_nlm_g.exe";
    private static final String EXE_UNIX = "tif_nlm_g";

    @Override
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        if (imp == null) { IJ.error("NLM Filter 3D Ext", "No image is open"); return DONE; }
        if (imp.getStackSize() < 2) {
            IJ.error("NLM Filter 3D Ext", "Requires an image stack (at least 2 slices)"); return DONE;
        }
        int bd = imp.getBitDepth();
        if (bd != 8 && bd != 16 && bd != 32) {
            IJ.error("NLM Filter 3D Ext", "Use 8-bit, 16-bit, or 32-bit float images."); return DONE;
        }
        return DOES_8G | DOES_16 | DOES_32 | STACK_REQUIRED | NO_CHANGES;
    }

    @Override
    public void run(ImageProcessor ip) {
        if (!showDialog()) return;
        long t0 = System.currentTimeMillis();
        IJ.log("=== 3D NLM Filter (External GPU) ===");
        IJ.log("Image: " + imp.getTitle() + " (" + imp.getWidth() + "x" + imp.getHeight() + "x" + imp.getStackSize() + ")");
        IJ.log("Patch radius: " + patchRadius + ", Search radius: " + searchRadius + ", h: " + h);

        ImagePlus result = processWithExternalTool();
        IJ.log("Total time: " + ((System.currentTimeMillis() - t0) / 1000.0) + " seconds");
        if (result != null) { result.setTitle(imp.getTitle() + "_NLM3D"); result.show(); }
    }

    private boolean showDialog() {
        GenericDialog gd = new GenericDialog("3D NLM Filter (External GPU)");
        gd.addNumericField("Patch radius (1-5):", patchRadius, 0);
        gd.addNumericField("Search radius (1-15):", searchRadius, 0);
        gd.addNumericField("h (filtering strength, -1=auto):", h, 1);
        gd.addStringField("Custom executable path:", customExePath, 40);
        gd.addCheckbox("Keep temporary files", keepTempFiles);
        gd.addMessage("h=-1: auto-estimate from noise level\n\n"
                + "Executable: " + (isWindows() ? EXE_WIN : EXE_UNIX));
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        patchRadius = (int) gd.getNextNumber();
        searchRadius = (int) gd.getNextNumber();
        h = gd.getNextNumber();
        customExePath = gd.getNextString();
        keepTempFiles = gd.getNextBoolean();
        if (patchRadius < 1 || patchRadius > 5) { IJ.error("Patch radius must be 1-5"); return false; }
        if (searchRadius < 1 || searchRadius > 15) { IJ.error("Search radius must be 1-15"); return false; }
        return true;
    }

    private static boolean isWindows() { return System.getProperty("os.name").toLowerCase().contains("win"); }

    private String findExecutable() {
        String exeName = isWindows() ? EXE_WIN : EXE_UNIX;
        if (customExePath != null && !customExePath.trim().isEmpty()) {
            File f = new File(customExePath.trim());
            if (f.exists() && f.canExecute()) return f.getAbsolutePath();
        }
        String d = IJ.getDirectory("plugins");
        if (d != null) { File f = new File(d, exeName); if (f.exists() && f.canExecute()) return f.getAbsolutePath(); }
        d = IJ.getDirectory("imagej");
        if (d != null) { File f = new File(d, exeName); if (f.exists() && f.canExecute()) return f.getAbsolutePath(); }
        return exeName;
    }

    private File getWorkingDirectory(File tempInputDir) {
        if (imp != null) {
            ij.io.FileInfo fi = imp.getOriginalFileInfo();
            if (fi != null && fi.directory != null && !fi.directory.isEmpty()) {
                File od = new File(fi.directory);
                File pd = od.getParentFile();
                if (pd != null && pd.exists()) return pd;
                if (od.exists()) return od;
            }
        }
        File p = tempInputDir.getParentFile();
        return (p != null && p.exists()) ? p : tempInputDir;
    }

    private ImagePlus processWithExternalTool() {
        File tid = null, tod = null;
        try {
            Path tb = Files.createTempDirectory("imagej_nlm3d_");
            tid = new File(tb.toFile(), "input"); tod = new File(tb.toFile(), "output");
            if (!tid.mkdirs() || !tod.mkdirs()) { IJ.error("Cannot create temp dirs"); return null; }
            IJ.showStatus("Saving stack..."); if (!saveStackAsTiffs(imp, tid)) { IJ.error("Failed to save stack"); return null; }
            String exe = findExecutable(); IJ.log("Executable: " + exe);
            IJ.showStatus("Running NLM filter...");
            int exitCode = runExternalTool(exe, tid, tod);
            if (exitCode != 0) { IJ.error("External tool error code: " + exitCode); return null; }
            IJ.showStatus("Reading results...");
            return new FolderOpener().openFolder(tod.getAbsolutePath());
        } catch (IOException e) { IJ.error("I/O error: " + e.getMessage()); return null;
        } finally {
            if (!keepTempFiles) { if (tid != null) deleteDir(tid); if (tod != null) deleteDir(tod);
                if (tid != null && tid.getParentFile() != null) tid.getParentFile().delete(); }
            IJ.showStatus(""); IJ.showProgress(1.0);
        }
    }

    private boolean saveStackAsTiffs(ImagePlus imp, File dir) {
        ImageStack stack = imp.getStack(); int depth = stack.getSize();
        int pw = Math.max(4, String.valueOf(depth).length());
        String fmt = "slice_%0" + pw + "d.tif";
        for (int z = 1; z <= depth; z++) {
            ImagePlus si = new ImagePlus(String.format(fmt, z), stack.getProcessor(z));
            if (!new FileSaver(si).saveAsTiff(new File(dir, String.format(fmt, z)).getAbsolutePath())) return false;
            if (z % 50 == 0 || z == depth) IJ.showProgress(z, depth);
        }
        return true;
    }

    private int runExternalTool(String exe, File inDir, File outDir) throws IOException {
        java.util.List<String> cmd = new java.util.ArrayList<String>();
        cmd.add(exe); cmd.add(inDir.getAbsolutePath()); cmd.add(outDir.getAbsolutePath());
        cmd.add(String.valueOf(patchRadius)); cmd.add(String.valueOf(searchRadius));
        if (h > 0) cmd.add(String.valueOf(h));

        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.redirectErrorStream(true);
        pb.directory(getWorkingDirectory(inDir));
        IJ.log("Command: " + String.join(" ", pb.command()));
        IJ.log("Working dir: " + pb.directory().getAbsolutePath());
        Process p = pb.start();
        try (BufferedReader r = new BufferedReader(new InputStreamReader(p.getInputStream()))) {
            String line; while ((line = r.readLine()) != null) IJ.log("[ext] " + line);
        }
        try { return p.waitFor(); } catch (InterruptedException e) { Thread.currentThread().interrupt(); p.destroy(); return -1; }
    }

    private void deleteDir(File d) {
        if (d == null || !d.exists()) return;
        File[] fs = d.listFiles();
        if (fs != null) for (File f : fs) { if (f.isDirectory()) deleteDir(f); else f.delete(); }
        d.delete();
    }
}