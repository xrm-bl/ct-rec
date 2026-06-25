/*
 * Gaussian_Filter_2D.java - 2D Gaussian Filter ImageJ Plugin with Preview
 *
 * Compile:
 *   javac -cp ij.jar Gaussian_Filter_2D.java
 *
 * Install:
 *   Place Gaussian_Filter_2D.class in ImageJ/plugins/
 */

import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;

import java.awt.AWTEvent;

public class Gaussian_Filter_2D implements ExtendedPlugInFilter, DialogListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;

    private ImagePlus imp;
    private int nPasses = 1;
    private double sigma = 2.0;
    private int kernelSize;
    private double[] kernelWeights;

    @Override
    public int setup(String arg, ImagePlus imp) {
        if ("final".equals(arg)) return DONE;
        this.imp = imp;
        if (imp == null) {
            IJ.error("Gaussian Filter 2D", "No image is open");
            return DONE;
        }
        return FLAGS;
    }

    @Override
    public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
        GenericDialog gd = new GenericDialog("2D Gaussian Filter");
        gd.addNumericField("Sigma:", sigma, 2);
        gd.addMessage("Kernel size: auto-calculated as 2*ceil(3*sigma)+1\n\n"
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
        sigma = gd.getNextNumber();
        if (gd.invalidNumber()) return false;
        if (sigma <= 0) return false;

        /* Recalculate kernel */
        kernelSize = 2 * (int) Math.ceil(3.0 * sigma) + 1;
        kernelWeights = precomputeGaussianKernel2D(kernelSize, sigma);
        return true;
    }

    @Override
    public void setNPasses(int nPasses) {
        this.nPasses = nPasses;
    }

    @Override
    public void run(ImageProcessor ip) {
        if (kernelWeights == null) {
            kernelSize = 2 * (int) Math.ceil(3.0 * sigma) + 1;
            kernelWeights = precomputeGaussianKernel2D(kernelSize, sigma);
        }
        applyGaussianFilter2D(ip);
    }

    private double[] precomputeGaussianKernel2D(int size, double sigma) {
        int half = size / 2;
        double[] weights = new double[size * size];
        double sigmaSq2 = 2.0 * sigma * sigma;
        double sum = 0.0;
        int idx = 0;
        for (int ky = -half; ky <= half; ky++) {
            for (int kx = -half; kx <= half; kx++) {
                double distSq = kx * kx + ky * ky;
                weights[idx] = Math.exp(-distSq / sigmaSq2);
                sum += weights[idx];
                idx++;
            }
        }
        for (idx = 0; idx < weights.length; idx++) weights[idx] /= sum;
        return weights;
    }

    private void applyGaussianFilter2D(ImageProcessor ip) {
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
            applyKernel8bit(input, output, width, height, roiX, roiY, roiW, roiH);
            System.arraycopy(output, 0, input, 0, input.length);
        } else if (pixels instanceof short[]) {
            short[] input = (short[]) pixels;
            short[] output = new short[input.length];
            System.arraycopy(input, 0, output, 0, input.length);
            applyKernel16bit(input, output, width, height, roiX, roiY, roiW, roiH);
            System.arraycopy(output, 0, input, 0, input.length);
        } else if (pixels instanceof float[]) {
            float[] input = (float[]) pixels;
            float[] output = new float[input.length];
            System.arraycopy(input, 0, output, 0, input.length);
            applyKernel32bit(input, output, width, height, roiX, roiY, roiW, roiH);
            System.arraycopy(output, 0, input, 0, input.length);
        }
    }

    private void applyKernel8bit(byte[] input, byte[] output,
            int width, int height, int roiX, int roiY, int roiW, int roiH) {

        final int half = kernelSize / 2;
        final double[] w = kernelWeights;
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
                            double sum = 0.0;
                            int widx = 0;
                            for (int ky = -half; ky <= half; ky++) {
                                int ny = y + ky;
                                for (int kx = -half; kx <= half; kx++) {
                                    int nx = x + kx;
                                    if (ny >= 0 && ny < fh && nx >= 0 && nx < fw) {
                                        sum += (double)(fi[ny * fw + nx] & 0xFF) * w[widx];
                                    }
                                    widx++;
                                }
                            }
                            int result = (int) Math.round(sum);
                            if (result < 0) result = 0;
                            if (result > 255) result = 255;
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

    private void applyKernel16bit(short[] input, short[] output,
            int width, int height, int roiX, int roiY, int roiW, int roiH) {

        final int half = kernelSize / 2;
        final double[] w = kernelWeights;
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
                            double sum = 0.0;
                            int widx = 0;
                            for (int ky = -half; ky <= half; ky++) {
                                int ny = y + ky;
                                for (int kx = -half; kx <= half; kx++) {
                                    int nx = x + kx;
                                    if (ny >= 0 && ny < fh && nx >= 0 && nx < fw) {
                                        sum += (double)(fi[ny * fw + nx] & 0xFFFF) * w[widx];
                                    }
                                    widx++;
                                }
                            }
                            int result = (int) Math.round(sum);
                            if (result < 0) result = 0;
                            if (result > 65535) result = 65535;
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

    private void applyKernel32bit(float[] input, float[] output,
            int width, int height, int roiX, int roiY, int roiW, int roiH) {

        final int half = kernelSize / 2;
        final double[] w = kernelWeights;
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
                            double sum = 0.0;
                            int widx = 0;
                            for (int ky = -half; ky <= half; ky++) {
                                int ny = y + ky;
                                for (int kx = -half; kx <= half; kx++) {
                                    int nx = x + kx;
                                    if (ny >= 0 && ny < fh && nx >= 0 && nx < fw) {
                                        sum += (double) fi[ny * fw + nx] * w[widx];
                                    }
                                    widx++;
                                }
                            }
                            fo[y * fw + x] = (float) sum;
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