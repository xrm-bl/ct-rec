/*
 * Bilateral_Filter_2D.java - 2D Bilateral Filter ImageJ Plugin with Preview
 *
 * Compile:
 *   javac -cp ij.jar Bilateral_Filter_2D.java
 *
 * Install:
 *   Place Bilateral_Filter_2D.class in ImageJ/plugins/
 *
 * Usage:
 *   1. Open an image in ImageJ
 *   2. Plugins > Bilateral Filter 2D
 *   3. Adjust parameters with live preview
 *   4. Use Image > Adjust > Brightness/Contrast during preview if needed
 */

import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;
import ij.process.FloatProcessor;

import java.awt.AWTEvent;

public class Bilateral_Filter_2D implements ExtendedPlugInFilter, DialogListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;

    private ImagePlus imp;
    private int nPasses = 1;
    private int kernelSize = 5;
    private double spatialSigma = 2.0;
    private double intensitySigma = -1.0;

    @Override
    public int setup(String arg, ImagePlus imp) {
        if ("final".equals(arg)) return DONE;
        this.imp = imp;
        if (imp == null) {
            IJ.error("Bilateral Filter 2D", "No image is open");
            return DONE;
        }
        /* Initialize default intensity sigma based on bit depth */
        if (intensitySigma <= 0) {
            intensitySigma = getDefaultIntensitySigma(imp);
        }
        return FLAGS;
    }

    @Override
    public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
        double defaultIntensitySigma = getDefaultIntensitySigma(imp);
        if (intensitySigma <= 0) {
            intensitySigma = defaultIntensitySigma;
        }

        GenericDialog gd = new NonBlockingGenericDialog("2D Bilateral Filter");
        gd.addNumericField("Kernel size (odd, 3-21):", kernelSize, 0);
        gd.addNumericField("Spatial sigma:", spatialSigma, 2);
        gd.addNumericField("Intensity sigma:", intensitySigma, 3);
        gd.addMessage("Auto-suggested intensity sigma:\n"
                + "  8-bit:  " + String.format("%.1f", 256.0 / 3.0) + "\n"
                + "  16-bit: " + String.format("%.1f", 65536.0 / 3.0) + "\n"
                + "  32-bit: based on image dynamic range / 3\n\n"
                + "Larger values = more smoothing\n"
                + "Smaller values = more edge preservation\n\n"
                + "Tip: Use Image > Adjust > Brightness/Contrast (Ctrl+Shift+C)\n"
                + "during preview to check filter effect.\n"
                + "Contrast changes affect display only, not the filter result.");

        gd.addPreviewCheckbox(pfr);
        gd.addDialogListener(this);
        gd.showDialog();

        if (gd.wasCanceled()) return DONE;

        IJ.register(this.getClass());
        return IJ.setupDialog(imp, FLAGS);
    }

    @Override
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        kernelSize = (int) gd.getNextNumber();
        spatialSigma = gd.getNextNumber();
        intensitySigma = gd.getNextNumber();

        if (gd.invalidNumber()) return false;
        if (kernelSize < 3 || kernelSize > 21 || kernelSize % 2 == 0) return false;
        if (spatialSigma <= 0 || intensitySigma <= 0) return false;

        return true;
    }

    @Override
    public void setNPasses(int nPasses) {
        this.nPasses = nPasses;
    }

    @Override
    public void run(ImageProcessor ip) {
        applyBilateralFilter2D(ip);
    }

    private double getDefaultIntensitySigma(ImagePlus imp) {
        int bitDepth = imp.getBitDepth();
        if (bitDepth == 8) {
            return 256.0 / 3.0;
        } else if (bitDepth == 16) {
            return 65536.0 / 3.0;
        } else if (bitDepth == 32) {
            return getFloatDynamicRange(imp) / 3.0;
        }
        return 50.0;
    }

    private double getFloatDynamicRange(ImagePlus imp) {
        ImageProcessor ip = imp.getProcessor();
        Object pixels = ip.getPixels();
        if (!(pixels instanceof float[])) return 1.0;

        float[] floatPixels = (float[]) pixels;
        int width = ip.getWidth();
        int height = ip.getHeight();
        int startY = height / 4, endY = 3 * height / 4;
        int startX = width / 4, endX = 3 * width / 4;
        int yStep = Math.max(1, (endY - startY) / 20);
        int xStep = Math.max(1, (endX - startX) / 20);

        float minVal = Float.MAX_VALUE, maxVal = -Float.MAX_VALUE;
        for (int y = startY; y < endY; y += yStep) {
            for (int x = startX; x < endX; x += xStep) {
                float val = floatPixels[y * width + x];
                if (!Float.isNaN(val) && !Float.isInfinite(val)) {
                    if (val < minVal) minVal = val;
                    if (val > maxVal) maxVal = val;
                }
            }
        }
        double range = maxVal - minVal;
        return range > 0.0 ? range : 1.0;
    }

    private void applyBilateralFilter2D(ImageProcessor ip) {
        int width = ip.getWidth();
        int height = ip.getHeight();
        Object pixels = ip.getPixels();

        java.awt.Rectangle roi = ip.getRoi();
        final int roiX = (roi != null) ? roi.x : 0;
        final int roiY = (roi != null) ? roi.y : 0;
        final int roiW = (roi != null) ? roi.width : width;
        final int roiH = (roi != null) ? roi.height : height;

        if (pixels instanceof byte[]) {
            byte[] input = (byte[]) pixels;
            byte[] output = new byte[input.length];
            System.arraycopy(input, 0, output, 0, input.length);
            applyBilateral8bit(input, output, width, height, roiX, roiY, roiW, roiH);
            System.arraycopy(output, 0, input, 0, input.length);
        } else if (pixels instanceof short[]) {
            short[] input = (short[]) pixels;
            short[] output = new short[input.length];
            System.arraycopy(input, 0, output, 0, input.length);
            applyBilateral16bit(input, output, width, height, roiX, roiY, roiW, roiH);
            System.arraycopy(output, 0, input, 0, input.length);
        } else if (pixels instanceof float[]) {
            float[] input = (float[]) pixels;
            float[] output = new float[input.length];
            System.arraycopy(input, 0, output, 0, input.length);
            applyBilateral32bit(input, output, width, height, roiX, roiY, roiW, roiH);
            System.arraycopy(output, 0, input, 0, input.length);
        }
    }

    private void applyBilateral8bit(byte[] input, byte[] output,
            int width, int height, int roiX, int roiY, int roiW, int roiH) {

        final int half = kernelSize / 2;
        final double spatialSigmaSq2 = 2.0 * spatialSigma * spatialSigma;
        final double intensitySigmaSq2 = 2.0 * intensitySigma * intensitySigma;
        final int fw = width, fh = height;
        final byte[] fi = input;
        final byte[] fo = output;
        final int rx = roiX, ry = roiY, rw = roiW, rh = roiH;

        int numThreads = Runtime.getRuntime().availableProcessors();
        Thread[] threads = new Thread[numThreads];

        for (int t = 0; t < numThreads; t++) {
            final int tid = t, tt = numThreads;
            threads[t] = new Thread(new Runnable() {
                @Override
                public void run() {
                    for (int y = ry + tid; y < ry + rh; y += tt) {
                        if (y < 0 || y >= fh) continue;
                        for (int x = rx; x < rx + rw; x++) {
                            if (x < 0 || x >= fw) continue;
                            double centerValue = (double)(fi[y * fw + x] & 0xFF);
                            double weightedSum = 0.0, weightSum = 0.0;
                            for (int ky = -half; ky <= half; ky++) {
                                int ny = y + ky;
                                if (ny < 0 || ny >= fh) continue;
                                for (int kx = -half; kx <= half; kx++) {
                                    int nx = x + kx;
                                    if (nx < 0 || nx >= fw) continue;
                                    double neighborValue = (double)(fi[ny * fw + nx] & 0xFF);
                                    double spatialDistSq = kx*kx + ky*ky;
                                    double spatialWeight = Math.exp(-spatialDistSq / spatialSigmaSq2);
                                    double intensityDiff = neighborValue - centerValue;
                                    double intensityWeight = Math.exp(-(intensityDiff * intensityDiff) / intensitySigmaSq2);
                                    double weight = spatialWeight * intensityWeight;
                                    weightedSum += neighborValue * weight;
                                    weightSum += weight;
                                }
                            }
                            int result;
                            if (weightSum > 0.0) {
                                result = (int) Math.round(weightedSum / weightSum);
                                if (result < 0) result = 0;
                                if (result > 255) result = 255;
                            } else {
                                result = (int) centerValue;
                            }
                            fo[y * fw + x] = (byte) result;
                        }
                    }
                }
            });
            threads[t].start();
        }
        for (Thread thread : threads) {
            try { thread.join(); } catch (InterruptedException e) { Thread.currentThread().interrupt(); }
        }
    }

    private void applyBilateral16bit(short[] input, short[] output,
            int width, int height, int roiX, int roiY, int roiW, int roiH) {

        final int half = kernelSize / 2;
        final double spatialSigmaSq2 = 2.0 * spatialSigma * spatialSigma;
        final double intensitySigmaSq2 = 2.0 * intensitySigma * intensitySigma;
        final int fw = width, fh = height;
        final short[] fi = input;
        final short[] fo = output;
        final int rx = roiX, ry = roiY, rw = roiW, rh = roiH;

        int numThreads = Runtime.getRuntime().availableProcessors();
        Thread[] threads = new Thread[numThreads];

        for (int t = 0; t < numThreads; t++) {
            final int tid = t, tt = numThreads;
            threads[t] = new Thread(new Runnable() {
                @Override
                public void run() {
                    for (int y = ry + tid; y < ry + rh; y += tt) {
                        if (y < 0 || y >= fh) continue;
                        for (int x = rx; x < rx + rw; x++) {
                            if (x < 0 || x >= fw) continue;
                            double centerValue = (double)(fi[y * fw + x] & 0xFFFF);
                            double weightedSum = 0.0, weightSum = 0.0;
                            for (int ky = -half; ky <= half; ky++) {
                                int ny = y + ky;
                                if (ny < 0 || ny >= fh) continue;
                                for (int kx = -half; kx <= half; kx++) {
                                    int nx = x + kx;
                                    if (nx < 0 || nx >= fw) continue;
                                    double neighborValue = (double)(fi[ny * fw + nx] & 0xFFFF);
                                    double spatialDistSq = kx*kx + ky*ky;
                                    double spatialWeight = Math.exp(-spatialDistSq / spatialSigmaSq2);
                                    double intensityDiff = neighborValue - centerValue;
                                    double intensityWeight = Math.exp(-(intensityDiff * intensityDiff) / intensitySigmaSq2);
                                    double weight = spatialWeight * intensityWeight;
                                    weightedSum += neighborValue * weight;
                                    weightSum += weight;
                                }
                            }
                            int result;
                            if (weightSum > 0.0) {
                                result = (int) Math.round(weightedSum / weightSum);
                                if (result < 0) result = 0;
                                if (result > 65535) result = 65535;
                            } else {
                                result = (int) centerValue;
                            }
                            fo[y * fw + x] = (short) result;
                        }
                    }
                }
            });
            threads[t].start();
        }
        for (Thread thread : threads) {
            try { thread.join(); } catch (InterruptedException e) { Thread.currentThread().interrupt(); }
        }
    }

    private void applyBilateral32bit(float[] input, float[] output,
            int width, int height, int roiX, int roiY, int roiW, int roiH) {

        final int half = kernelSize / 2;
        final double spatialSigmaSq2 = 2.0 * spatialSigma * spatialSigma;
        final double intensitySigmaSq2 = 2.0 * intensitySigma * intensitySigma;
        final int fw = width, fh = height;
        final float[] fi = input;
        final float[] fo = output;
        final int rx = roiX, ry = roiY, rw = roiW, rh = roiH;

        int numThreads = Runtime.getRuntime().availableProcessors();
        Thread[] threads = new Thread[numThreads];

        for (int t = 0; t < numThreads; t++) {
            final int tid = t, tt = numThreads;
            threads[t] = new Thread(new Runnable() {
                @Override
                public void run() {
                    for (int y = ry + tid; y < ry + rh; y += tt) {
                        if (y < 0 || y >= fh) continue;
                        for (int x = rx; x < rx + rw; x++) {
                            if (x < 0 || x >= fw) continue;
                            double centerValue = (double) fi[y * fw + x];
                            double weightedSum = 0.0, weightSum = 0.0;
                            for (int ky = -half; ky <= half; ky++) {
                                int ny = y + ky;
                                if (ny < 0 || ny >= fh) continue;
                                for (int kx = -half; kx <= half; kx++) {
                                    int nx = x + kx;
                                    if (nx < 0 || nx >= fw) continue;
                                    double neighborValue = (double) fi[ny * fw + nx];
                                    double spatialDistSq = kx*kx + ky*ky;
                                    double spatialWeight = Math.exp(-spatialDistSq / spatialSigmaSq2);
                                    double intensityDiff = neighborValue - centerValue;
                                    double intensityWeight = Math.exp(-(intensityDiff * intensityDiff) / intensitySigmaSq2);
                                    double weight = spatialWeight * intensityWeight;
                                    weightedSum += neighborValue * weight;
                                    weightSum += weight;
                                }
                            }
                            fo[y * fw + x] = (weightSum > 0.0)
                                ? (float)(weightedSum / weightSum)
                                : (float) centerValue;
                        }
                    }
                }
            });
            threads[t].start();
        }
        for (Thread thread : threads) {
            try { thread.join(); } catch (InterruptedException e) { Thread.currentThread().interrupt(); }
        }
    }
}