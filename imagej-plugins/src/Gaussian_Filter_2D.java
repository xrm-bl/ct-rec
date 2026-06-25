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

public class Gaussian_Filter_2D implements ExtendedPlugInFilter, DialogListener, ImageListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;

    private ImagePlus imp;
    private int nPasses = 1;
    private double sigma = 2.0;
    private int kernelSize;
    private double[] kernelWeights;

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
            IJ.error("Gaussian Filter 2D", "No image is open");
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

        GenericDialog gd = new NonBlockingGenericDialog("2D Gaussian Filter");
        gd.addNumericField("Sigma:", sigma, 2);
        gd.addMessage("Kernel size: auto-calculated as 2*ceil(3*sigma)+1\n\n"
                + "Tip: Use Image > Adjust > Brightness/Contrast (Ctrl+Shift+C)\n"
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

    /* ImageListener: re-apply the user's display range whenever a preview
     * render or the preview-off reset would otherwise revert it. */
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
            /* Redraw right after a filter pass: the runner's reset() reverted
             * min/max to the snapshot range, so re-apply the preserved one. */
            pendingRender = false;
            applyKeep(curMin, curMax);
        } else if (curMin == origMin && curMax == origMax
                && (keepMin != origMin || keepMax != origMax)) {
            /* Redraw that reverted to the snapshot range without a filter pass:
             * preview turned off, or the Cancel cleanup restoring the image.
             * Re-apply the preserved contrast instead. */
            applyKeep(curMin, curMax);
            if (closing) EventQueue.invokeLater(this::detach);
        } else {
            /* Any other redraw is a genuine user contrast change (e.g.
             * Brightness/Contrast): adopt it as the contrast to preserve. */
            keepMin = curMin;
            keepMax = curMax;
        }
    }

    /* Re-apply the preserved display range if the current one drifted. */
    private void applyKeep(double curMin, double curMax) {
        if (curMin != keepMin || curMax != keepMax) {
            enforcing = true;
            this.imp.setDisplayRange(keepMin, keepMax);
            this.imp.updateAndDraw();
            enforcing = false;
        }
    }

    /* Stop watching the image. */
    private void detach() {
        if (!rangeValid) return;
        rangeValid = false;
        ImagePlus.removeImageListener(this);
    }

    @Override
    public void run(ImageProcessor ip) {
        if (kernelWeights == null) {
            kernelSize = 2 * (int) Math.ceil(3.0 * sigma) + 1;
            kernelWeights = precomputeGaussianKernel2D(kernelSize, sigma);
        }
        applyGaussianFilter2D(ip);
        /* Mark that the redraw following this filter pass should keep the
         * preserved contrast (the runner's reset() reverts min/max). */
        pendingRender = true;
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