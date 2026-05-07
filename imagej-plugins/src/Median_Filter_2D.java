/*
 * Median_Filter_2D.java - 2D Median Filter ImageJ Plugin with Preview
 *
 * Compile:
 *   javac -cp ij.jar Median_Filter_2D.java
 *
 * Install:
 *   Place Median_Filter_2D.class in ImageJ/plugins/
 */

import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;

import java.awt.AWTEvent;
import java.util.Arrays;

public class Median_Filter_2D implements ExtendedPlugInFilter, DialogListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;

    private ImagePlus imp;
    private int nPasses = 1;
    private int kernelSize = 3;

    @Override
    public int setup(String arg, ImagePlus imp) {
        if ("final".equals(arg)) return DONE;
        this.imp = imp;
        if (imp == null) {
            IJ.error("Median Filter 2D", "No image is open");
            return DONE;
        }
        return FLAGS;
    }

    @Override
    public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
        GenericDialog gd = new GenericDialog("2D Median Filter");
        gd.addNumericField("Kernel size (odd, 3-21):", kernelSize, 0);
        gd.addMessage("Tip: Use Image > Adjust > Brightness/Contrast (Ctrl+Shift+C)\n"
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
        if (gd.invalidNumber()) return false;
        if (kernelSize < 3 || kernelSize > 21 || kernelSize % 2 == 0) return false;
        return true;
    }

    @Override
    public void setNPasses(int nPasses) {
        this.nPasses = nPasses;
    }

    @Override
    public void run(ImageProcessor ip) {
        applyMedianFilter2D(ip);
    }

    private void applyMedianFilter2D(ImageProcessor ip) {
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
            applyMedian8bit(input, output, width, height, roiX, roiY, roiW, roiH);
            System.arraycopy(output, 0, input, 0, input.length);
        } else if (pixels instanceof short[]) {
            short[] input = (short[]) pixels;
            short[] output = new short[input.length];
            System.arraycopy(input, 0, output, 0, input.length);
            applyMedian16bit(input, output, width, height, roiX, roiY, roiW, roiH);
            System.arraycopy(output, 0, input, 0, input.length);
        } else if (pixels instanceof float[]) {
            float[] input = (float[]) pixels;
            float[] output = new float[input.length];
            System.arraycopy(input, 0, output, 0, input.length);
            applyMedian32bit(input, output, width, height, roiX, roiY, roiW, roiH);
            System.arraycopy(output, 0, input, 0, input.length);
        }
    }

    private void applyMedian8bit(byte[] input, byte[] output,
            int width, int height, int roiX, int roiY, int roiW, int roiH) {

        final int half = kernelSize / 2;
        final int maxNeighbors = kernelSize * kernelSize;
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
                    int[] neighbors = new int[maxNeighbors];
                    for (int y = ry + tid; y < ry + rh; y += tt) {
                        if (y < 0 || y >= fh) continue;
                        for (int x = rx; x < rx + rw; x++) {
                            if (x < 0 || x >= fw) continue;
                            int count = 0;
                            for (int ky = -half; ky <= half; ky++) {
                                int ny = y + ky;
                                if (ny < 0 || ny >= fh) continue;
                                for (int kx = -half; kx <= half; kx++) {
                                    int nx = x + kx;
                                    if (nx < 0 || nx >= fw) continue;
                                    neighbors[count++] = fi[ny * fw + nx] & 0xFF;
                                }
                            }
                            Arrays.sort(neighbors, 0, count);
                            fo[y * fw + x] = (byte) neighbors[count / 2];
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

    private void applyMedian16bit(short[] input, short[] output,
            int width, int height, int roiX, int roiY, int roiW, int roiH) {

        final int half = kernelSize / 2;
        final int maxNeighbors = kernelSize * kernelSize;
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
                    int[] neighbors = new int[maxNeighbors];
                    for (int y = ry + tid; y < ry + rh; y += tt) {
                        if (y < 0 || y >= fh) continue;
                        for (int x = rx; x < rx + rw; x++) {
                            if (x < 0 || x >= fw) continue;
                            int count = 0;
                            for (int ky = -half; ky <= half; ky++) {
                                int ny = y + ky;
                                if (ny < 0 || ny >= fh) continue;
                                for (int kx = -half; kx <= half; kx++) {
                                    int nx = x + kx;
                                    if (nx < 0 || nx >= fw) continue;
                                    neighbors[count++] = fi[ny * fw + nx] & 0xFFFF;
                                }
                            }
                            Arrays.sort(neighbors, 0, count);
                            fo[y * fw + x] = (short) neighbors[count / 2];
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

    private void applyMedian32bit(float[] input, float[] output,
            int width, int height, int roiX, int roiY, int roiW, int roiH) {

        final int half = kernelSize / 2;
        final int maxNeighbors = kernelSize * kernelSize;
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
                    float[] neighbors = new float[maxNeighbors];
                    for (int y = ry + tid; y < ry + rh; y += tt) {
                        if (y < 0 || y >= fh) continue;
                        for (int x = rx; x < rx + rw; x++) {
                            if (x < 0 || x >= fw) continue;
                            int count = 0;
                            for (int ky = -half; ky <= half; ky++) {
                                int ny = y + ky;
                                if (ny < 0 || ny >= fh) continue;
                                for (int kx = -half; kx <= half; kx++) {
                                    int nx = x + kx;
                                    if (nx < 0 || nx >= fw) continue;
                                    neighbors[count++] = fi[ny * fw + nx];
                                }
                            }
                            Arrays.sort(neighbors, 0, count);
                            fo[y * fw + x] = neighbors[count / 2];
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