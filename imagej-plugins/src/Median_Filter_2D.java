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
import ij.ImageListener;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;

import java.awt.AWTEvent;
import java.awt.EventQueue;

import javax.swing.Timer;
import java.util.Arrays;

public class Median_Filter_2D implements ExtendedPlugInFilter, DialogListener, ImageListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;

    private ImagePlus imp;
    private int nPasses = 1;
    private int kernelSize = 3;

    /* Preserve the user-set display range (contrast) across preview on/off and
     * Cancel. ImageProcessor.snapshot()/reset() also save & restore min/max, so
     * the runner's reset() on preview toggle (or the Cancel cleanup) would
     * otherwise revert the contrast.
     *   keepMin/keepMax - the contrast to preserve (tracks the user's changes)
     *   origMin/origMax - the snapshot range that reset() reverts to */
    private double keepMin, keepMax;
    private double origMin, origMax;
    private boolean rangeValid = false;
    private boolean enforcing = false;
    private boolean pendingRender = false;
    private boolean closing = false;

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
        /* Remember the current display range and watch for redraws so the
         * user-set contrast survives preview toggling and Cancel. */
        keepMin = origMin = imp.getDisplayRangeMin();
        keepMax = origMax = imp.getDisplayRangeMax();
        rangeValid = true;
        pendingRender = false;
        closing = false;
        ImagePlus.addImageListener(this);

        GenericDialog gd = new NonBlockingGenericDialog("2D Median Filter");
        gd.addNumericField("Kernel size (odd, 3-21):", kernelSize, 0);
        gd.addMessage("Tip: Use Image > Adjust > Brightness/Contrast (Ctrl+Shift+C)\n"
                + "during preview to check filter effect.\n"
                + "Contrast changes affect display only, not the filter result.");

        gd.addPreviewCheckbox(pfr);
        gd.addDialogListener(this);
        gd.showDialog();

        /* Keep the listener attached briefly: the runner's post-dialog reset
         * (Cancel with preview on) happens on another thread after this method
         * returns. The listener re-applies the contrast when it sees that
         * reset; this timer then detaches as a fallback. */
        closing = true;
        Timer detachTimer = new Timer(1000, e -> detach());
        detachTimer.setRepeats(false);
        detachTimer.start();

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

    /* ImageListener: keep the user's display range (contrast) across preview
     * renders, preview off, and Cancel (the runner's reset() reverts min/max). */
    @Override
    public void imageOpened(ImagePlus imp) {}

    @Override
    public void imageClosed(ImagePlus imp) {}

    @Override
    public void imageUpdated(ImagePlus updated) {
        if (updated != this.imp || !rangeValid || enforcing) return;

        double curMin = this.imp.getDisplayRangeMin();
        double curMax = this.imp.getDisplayRangeMax();

        if (pendingRender) {
            pendingRender = false;
            applyKeep(curMin, curMax);
        } else if (curMin == origMin && curMax == origMax
                && (keepMin != origMin || keepMax != origMax)) {
            applyKeep(curMin, curMax);
            if (closing) EventQueue.invokeLater(this::detach);
        } else {
            keepMin = curMin;
            keepMax = curMax;
        }
    }

    private void applyKeep(double curMin, double curMax) {
        if (curMin != keepMin || curMax != keepMax) {
            enforcing = true;
            this.imp.setDisplayRange(keepMin, keepMax);
            this.imp.updateAndDraw();
            enforcing = false;
        }
    }

    private void detach() {
        if (!rangeValid) return;
        rangeValid = false;
        ImagePlus.removeImageListener(this);
    }

    @Override
    public void run(ImageProcessor ip) {
        applyMedianFilter2D(ip);
        pendingRender = true;
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