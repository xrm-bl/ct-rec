/*
 * Bilateral_Filter_3D_Ext.java - 3D Bilateral Filter ImageJ Plugin
 * (External Process Version)
 *
 * Compile:
 *   javac -cp ij.jar Bilateral_Filter_3D_Ext.java
 *
 * Install:
 *   Place Bilateral_Filter_3D_Ext.class and tif_blf/tif_blf_g in plugins/
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

public class Bilateral_Filter_3D_Ext implements PlugInFilter {

    private ImagePlus imp;
    private int kernelSize = 5;
    private double spatialSigma = 2.0;
    private double intensitySigma = -1.0;
    private boolean useGpu = false;
    private String customExePath = "";
    private boolean keepTempFiles = false;

    private static final String CPU_EXE_WIN = "tif_blf.exe";
    private static final String CPU_EXE_UNIX = "tif_blf";
    private static final String GPU_EXE_WIN = "tif_blf_g.exe";
    private static final String GPU_EXE_UNIX = "tif_blf_g";

    @Override
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        if (imp == null) {
            IJ.error("Bilateral Filter 3D Ext", "No image is open");
            return DONE;
        }
        if (imp.getStackSize() < 2) {
            IJ.error("Bilateral Filter 3D Ext",
                "This plugin requires an image stack (at least 2 slices)");
            return DONE;
        }
        int bitDepth = imp.getBitDepth();
        if (bitDepth != 8 && bitDepth != 16 && bitDepth != 32) {
            IJ.error("Bilateral Filter 3D Ext",
                "Unsupported bit depth. Use 8-bit, 16-bit, or 32-bit float images.");
            return DONE;
        }
        return DOES_8G | DOES_16 | DOES_32 | STACK_REQUIRED | NO_CHANGES;
    }

    @Override
    public void run(ImageProcessor ip) {
        if (!showDialog()) return;

        long startTime = System.currentTimeMillis();
        IJ.log("=== 3D Bilateral Filter (External) ===");
        IJ.log("Image: " + imp.getTitle());
        IJ.log("Dimensions: " + imp.getWidth() + " x " + imp.getHeight()
                + " x " + imp.getStackSize());
        IJ.log("Bit depth: " + imp.getBitDepth());
        IJ.log("Mode: " + (useGpu ? "GPU (CUDA)" : "CPU (OpenMP)"));
        IJ.log("Kernel size: " + kernelSize);
        IJ.log("Spatial sigma: " + spatialSigma);
        IJ.log("Intensity sigma: " + intensitySigma);

        ImagePlus result = processWithExternalTool();

        long elapsed = System.currentTimeMillis() - startTime;
        IJ.log("Total time: " + (elapsed / 1000.0) + " seconds");

        if (result != null) {
            result.setTitle(imp.getTitle() + "_BLF3D");
            result.show();
        }
    }

    private boolean showDialog() {
        double defaultIntensitySigma = getDefaultIntensitySigma(imp);
        if (intensitySigma <= 0) intensitySigma = defaultIntensitySigma;

        GenericDialog gd = new GenericDialog("3D Bilateral Filter (External)");
        gd.addCheckbox("Use GPU (CUDA) version", useGpu);
        gd.addNumericField("Kernel size (odd, 3-21):", kernelSize, 0);
        gd.addNumericField("Spatial sigma:", spatialSigma, 2);
        gd.addNumericField("Intensity sigma:", intensitySigma, 3);
        gd.addStringField("Custom executable path (optional):", customExePath, 40);
        gd.addCheckbox("Keep temporary files (for debugging)", keepTempFiles);
        gd.addMessage("Required executables (in plugins/ or PATH):\n"
                + "  CPU: " + (isWindows() ? CPU_EXE_WIN : CPU_EXE_UNIX) + "\n"
                + "  GPU: " + (isWindows() ? GPU_EXE_WIN : GPU_EXE_UNIX));

        gd.showDialog();
        if (gd.wasCanceled()) return false;

        useGpu = gd.getNextBoolean();
        kernelSize = (int) gd.getNextNumber();
        spatialSigma = gd.getNextNumber();
        intensitySigma = gd.getNextNumber();
        customExePath = gd.getNextString();
        keepTempFiles = gd.getNextBoolean();

        if (kernelSize < 3 || kernelSize > 21 || kernelSize % 2 == 0) {
            IJ.error("Bilateral Filter 3D Ext", "Kernel size must be odd and between 3 and 21");
            return false;
        }
        if (spatialSigma <= 0 || intensitySigma <= 0) {
            IJ.error("Bilateral Filter 3D Ext", "Sigma values must be positive");
            return false;
        }
        return true;
    }

    private double getDefaultIntensitySigma(ImagePlus imp) {
        int bitDepth = imp.getBitDepth();
        if (bitDepth == 8) return 256.0 / 3.0;
        else if (bitDepth == 16) return 65536.0 / 3.0;
        else if (bitDepth == 32) return getFloatDynamicRange(imp) / 3.0;
        return 50.0;
    }

    private double getFloatDynamicRange(ImagePlus imp) {
        ImageStack stack = imp.getStack();
        int depth = stack.getSize();
        int sampleCount = Math.max(1, Math.min(10, depth / 20));
        int sampleInterval = Math.max(1, depth / sampleCount);
        float minVal = Float.MAX_VALUE, maxVal = -Float.MAX_VALUE;
        int width = imp.getWidth(), height = imp.getHeight();
        int startY = height / 4, endY = 3 * height / 4;
        int startX = width / 4, endX = 3 * width / 4;
        int yStep = Math.max(1, (endY - startY) / 10);
        for (int i = 0; i < sampleCount; i++) {
            int sliceIdx = i * sampleInterval + 1;
            if (sliceIdx > depth) break;
            ImageProcessor ip = stack.getProcessor(sliceIdx);
            Object pixels = ip.getPixels();
            if (!(pixels instanceof float[])) continue;
            float[] floatPixels = (float[]) pixels;
            for (int y = startY; y < endY; y += yStep) {
                for (int x = startX; x < endX; x++) {
                    float val = floatPixels[y * width + x];
                    if (!Float.isNaN(val) && !Float.isInfinite(val)) {
                        if (val < minVal) minVal = val;
                        if (val > maxVal) maxVal = val;
                    }
                }
            }
        }
        double range = maxVal - minVal;
        return range > 0.0 ? range : 1.0;
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
            Path tempBase = Files.createTempDirectory("imagej_blf3d_");
            tempInputDir = new File(tempBase.toFile(), "input");
            tempOutputDir = new File(tempBase.toFile(), "output");
            if (!tempInputDir.mkdirs() || !tempOutputDir.mkdirs()) {
                IJ.error("Bilateral Filter 3D Ext", "Cannot create temporary directories");
                return null;
            }
            IJ.showStatus("Saving stack as TIFF files...");
            if (!saveStackAsTiffs(imp, tempInputDir)) {
                IJ.error("Bilateral Filter 3D Ext", "Failed to save stack");
                return null;
            }
            String exePath = findExecutable();
            IJ.log("Executable: " + exePath);
            IJ.showStatus("Running external bilateral filter...");
            int exitCode = runExternalTool(exePath, tempInputDir, tempOutputDir);
            if (exitCode != 0) {
                IJ.error("Bilateral Filter 3D Ext",
                    "External tool exited with error code: " + exitCode);
                return null;
            }
            IJ.showStatus("Reading filtered TIFF files...");
            FolderOpener folderOpener = new FolderOpener();
            return folderOpener.openFolder(tempOutputDir.getAbsolutePath());
        } catch (IOException e) {
            IJ.error("Bilateral Filter 3D Ext", "I/O error: " + e.getMessage());
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
            String.valueOf(kernelSize),
            String.valueOf(spatialSigma),
            String.valueOf(intensitySigma)
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