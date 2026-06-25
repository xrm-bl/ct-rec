/*
 * Gaussian_Filter_3D_Ext.java - 3D Gaussian Filter ImageJ Plugin
 * (External Process Version)
 *
 * Compile:
 *   javac -cp ij.jar Gaussian_Filter_3D_Ext.java
 *
 * Install:
 *   Place Gaussian_Filter_3D_Ext.class and tif_gsf/tif_gsf_g in plugins/
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

public class Gaussian_Filter_3D_Ext implements PlugInFilter {

    private ImagePlus imp;
    private double sigma = 2.0;
    private boolean useGpu = false;
    private String customExePath = "";
    private boolean keepTempFiles = false;

    private static final String CPU_EXE_WIN = "tif_gsf.exe";
    private static final String CPU_EXE_UNIX = "tif_gsf";
    private static final String GPU_EXE_WIN = "tif_gsf_g.exe";
    private static final String GPU_EXE_UNIX = "tif_gsf_g";

    @Override
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        if (imp == null) {
            IJ.error("Gaussian Filter 3D Ext", "No image is open");
            return DONE;
        }
        if (imp.getStackSize() < 2) {
            IJ.error("Gaussian Filter 3D Ext",
                "This plugin requires an image stack (at least 2 slices)");
            return DONE;
        }
        int bitDepth = imp.getBitDepth();
        if (bitDepth != 8 && bitDepth != 16 && bitDepth != 32) {
            IJ.error("Gaussian Filter 3D Ext",
                "Unsupported bit depth. Use 8-bit, 16-bit, or 32-bit float images.");
            return DONE;
        }
        return DOES_8G | DOES_16 | DOES_32 | STACK_REQUIRED | NO_CHANGES;
    }

    @Override
    public void run(ImageProcessor ip) {
        if (!showDialog()) return;

        long startTime = System.currentTimeMillis();
        int kernelSize = 2 * (int) Math.ceil(3.0 * sigma) + 1;

        IJ.log("=== 3D Gaussian Filter (External) ===");
        IJ.log("Image: " + imp.getTitle());
        IJ.log("Dimensions: " + imp.getWidth() + " x " + imp.getHeight()
                + " x " + imp.getStackSize());
        IJ.log("Sigma: " + sigma);
        IJ.log("Auto kernel size: " + kernelSize);
        IJ.log("Mode: " + (useGpu ? "GPU (CUDA)" : "CPU (OpenMP)"));

        ImagePlus result = processWithExternalTool();

        long elapsed = System.currentTimeMillis() - startTime;
        IJ.log("Total time: " + (elapsed / 1000.0) + " seconds");

        if (result != null) {
            result.setTitle(imp.getTitle() + "_GSF3D");
            result.show();
        }
    }

    private boolean showDialog() {
        GenericDialog gd = new GenericDialog("3D Gaussian Filter (External)");
        gd.addCheckbox("Use GPU (CUDA) version", useGpu);
        gd.addNumericField("Sigma:", sigma, 2);
        gd.addStringField("Custom executable path (optional):", customExePath, 40);
        gd.addCheckbox("Keep temporary files (for debugging)", keepTempFiles);
        gd.addMessage("Kernel size is auto-calculated: 2*ceil(3*sigma)+1\n\n"
                + "Required executables (in plugins/ or PATH):\n"
                + "  CPU: " + (isWindows() ? CPU_EXE_WIN : CPU_EXE_UNIX) + "\n"
                + "  GPU: " + (isWindows() ? GPU_EXE_WIN : GPU_EXE_UNIX));

        gd.showDialog();
        if (gd.wasCanceled()) return false;

        useGpu = gd.getNextBoolean();
        sigma = gd.getNextNumber();
        customExePath = gd.getNextString();
        keepTempFiles = gd.getNextBoolean();

        if (sigma <= 0) {
            IJ.error("Gaussian Filter 3D Ext", "Sigma must be positive");
            return false;
        }
        return true;
    }

    private static boolean isWindows() {
        return System.getProperty("os.name").toLowerCase().contains("win");
    }

    private String findExecutable() {
        String exeName = useGpu
            ? (isWindows() ? GPU_EXE_WIN : GPU_EXE_UNIX)
            : (isWindows() ? CPU_EXE_WIN : CPU_EXE_UNIX);

        if (customExePath != null && !customExePath.trim().isEmpty()) {
            File f = new File(customExePath.trim());
            if (f.exists() && f.canExecute()) return f.getAbsolutePath();
        }
        String pluginsDir = IJ.getDirectory("plugins");
        if (pluginsDir != null) {
            File f = new File(pluginsDir, exeName);
            if (f.exists() && f.canExecute()) return f.getAbsolutePath();
        }
        String imagejDir = IJ.getDirectory("imagej");
        if (imagejDir != null) {
            File f = new File(imagejDir, exeName);
            if (f.exists() && f.canExecute()) return f.getAbsolutePath();
        }
        return exeName;
    }

    private File getWorkingDirectory(File tempInputDir) {
        if (imp != null) {
            ij.io.FileInfo fi = imp.getOriginalFileInfo();
            if (fi != null && fi.directory != null && !fi.directory.isEmpty()) {
                File originalDir = new File(fi.directory);
                File parentDir = originalDir.getParentFile();
                if (parentDir != null && parentDir.exists()) return parentDir;
                if (originalDir.exists()) return originalDir;
            }
        }
        File parent = tempInputDir.getParentFile();
        return (parent != null && parent.exists()) ? parent : tempInputDir;
    }

    private ImagePlus processWithExternalTool() {
        File tempInputDir = null;
        File tempOutputDir = null;
        try {
            Path tempBase = Files.createTempDirectory("imagej_gsf3d_");
            tempInputDir = new File(tempBase.toFile(), "input");
            tempOutputDir = new File(tempBase.toFile(), "output");
            if (!tempInputDir.mkdirs() || !tempOutputDir.mkdirs()) {
                IJ.error("Gaussian Filter 3D Ext", "Cannot create temporary directories");
                return null;
            }
            IJ.showStatus("Saving stack as TIFF files...");
            if (!saveStackAsTiffs(imp, tempInputDir)) {
                IJ.error("Gaussian Filter 3D Ext", "Failed to save stack");
                return null;
            }
            String exePath = findExecutable();
            IJ.log("Executable: " + exePath);
            IJ.showStatus("Running external Gaussian filter...");
            int exitCode = runExternalTool(exePath, tempInputDir, tempOutputDir);
            if (exitCode != 0) {
                IJ.error("Gaussian Filter 3D Ext",
                    "External tool exited with error code: " + exitCode);
                return null;
            }
            IJ.showStatus("Reading filtered TIFF files...");
            FolderOpener folderOpener = new FolderOpener();
            return folderOpener.openFolder(tempOutputDir.getAbsolutePath());
        } catch (IOException e) {
            IJ.error("Gaussian Filter 3D Ext", "I/O error: " + e.getMessage());
            return null;
        } finally {
            if (!keepTempFiles) {
                if (tempInputDir != null) deleteDirectory(tempInputDir);
                if (tempOutputDir != null) deleteDirectory(tempOutputDir);
                if (tempInputDir != null && tempInputDir.getParentFile() != null)
                    tempInputDir.getParentFile().delete();
            }
            IJ.showStatus("");
            IJ.showProgress(1.0);
        }
    }

    private boolean saveStackAsTiffs(ImagePlus imp, File dir) {
        ImageStack stack = imp.getStack();
        int depth = stack.getSize();
        int padWidth = String.valueOf(depth).length();
        if (padWidth < 4) padWidth = 4;
        String formatStr = "slice_%0" + padWidth + "d.tif";
        for (int z = 1; z <= depth; z++) {
            ImageProcessor ip = stack.getProcessor(z);
            ImagePlus sliceImp = new ImagePlus(String.format(formatStr, z), ip);
            File outFile = new File(dir, String.format(formatStr, z));
            FileSaver saver = new FileSaver(sliceImp);
            if (!saver.saveAsTiff(outFile.getAbsolutePath())) return false;
            if (z % 50 == 0 || z == depth) IJ.showProgress(z, depth);
        }
        return true;
    }

    private int runExternalTool(String exePath, File inputDir, File outputDir)
            throws IOException {
        ProcessBuilder pb = new ProcessBuilder(
            exePath,
            inputDir.getAbsolutePath(),
            outputDir.getAbsolutePath(),
            String.valueOf(sigma)
        );
        pb.redirectErrorStream(true);

        File workingDir = getWorkingDirectory(inputDir);
        pb.directory(workingDir);

        IJ.log("Command: " + String.join(" ", pb.command()));
        IJ.log("Working directory: " + workingDir.getAbsolutePath());

        Process process = pb.start();
        try (BufferedReader reader = new BufferedReader(
                new InputStreamReader(process.getInputStream()))) {
            String line;
            while ((line = reader.readLine()) != null) IJ.log("[ext] " + line);
        }
        try {
            return process.waitFor();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            process.destroy();
            return -1;
        }
    }

    private void deleteDirectory(File dir) {
        if (dir == null || !dir.exists()) return;
        File[] files = dir.listFiles();
        if (files != null) {
            for (File f : files) {
                if (f.isDirectory()) deleteDirectory(f);
                else f.delete();
            }
        }
        dir.delete();
    }
}