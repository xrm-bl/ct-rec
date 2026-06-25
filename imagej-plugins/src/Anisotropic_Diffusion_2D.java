/*
 * Anisotropic_Diffusion_2D.java - 2D Perona-Malik Anisotropic Diffusion
 * ImageJ Plugin with Preview
 *
 * Compile: javac -cp ij.jar Anisotropic_Diffusion_2D.java
 * Install: Place Anisotropic_Diffusion_2D.class in ImageJ/plugins/
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

public class Anisotropic_Diffusion_2D implements ExtendedPlugInFilter, DialogListener, ImageListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;

    private ImagePlus imp;
    private int nPasses = 1;
    private int iterations = 20;
    private double K = -1.0;
    private double dt = 0.2;
    private int mode = 1;

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
        if (imp == null) { IJ.error("Anisotropic Diffusion 2D", "No image is open"); return DONE; }
        if (K <= 0) K = estimateK(imp);
        return FLAGS;
    }

    @Override
    public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
        if (K <= 0) K = estimateK(imp);

        /* Remember the current display range and watch for redraws so the
         * user-set contrast survives preview toggling and Cancel. */
        keepMin = origMin = imp.getDisplayRangeMin();
        keepMax = origMax = imp.getDisplayRangeMax();
        rangeValid = true;
        pendingRender = false;
        closing = false;
        ImagePlus.addImageListener(this);

        GenericDialog gd = new NonBlockingGenericDialog("2D Anisotropic Diffusion (Perona-Malik)");
        gd.addNumericField("Iterations:", iterations, 0);
        gd.addNumericField("K (edge threshold):", K, 2);
        gd.addNumericField("dt (time step, <0.25):", dt, 3);
        gd.addChoice("Mode:", new String[]{"1: 1/(1+(s/K)^2)", "2: exp(-(s/K)^2)"}, "1: 1/(1+(s/K)^2)");
        gd.addMessage("More iterations = more smoothing\n"
                + "Larger K = less edge preservation\n\n"
                + "Tip: Ctrl+Shift+C for Brightness/Contrast during preview.");

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
        iterations = (int) gd.getNextNumber();
        K = gd.getNextNumber();
        dt = gd.getNextNumber();
        String modeStr = gd.getNextChoice();
        mode = modeStr.startsWith("1") ? 1 : 2;
        if (gd.invalidNumber()) return false;
        if (iterations < 1 || iterations > 5000) return false;
        if (K <= 0 || dt <= 0 || dt >= 0.25) return false;
        return true;
    }

    @Override
    public void setNPasses(int nPasses) { this.nPasses = nPasses; }

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
        applyAD2D(ip);
        pendingRender = true;
    }

    private double estimateK(ImagePlus imp) {
        ImageProcessor ip = imp.getProcessor();
        int w = ip.getWidth(), h = ip.getHeight();
        int sx = w / 4, ex = 3 * w / 4, sy = h / 4, ey = 3 * h / 4;
        int step = Math.max(1, (ey - sy) / 20);
        double sumDiffSq = 0; long count = 0;
        for (int y = sy; y < ey; y += step) {
            for (int x = sx; x < ex - 1; x++) {
                double diff = ip.getf(x + 1, y) - ip.getf(x, y);
                sumDiffSq += diff * diff; count++;
            }
        }
        if (count < 2) return 50.0;
        double sigma = Math.sqrt(sumDiffSq / count / 2.0);
        return 2.0 * sigma;
    }

    private void applyAD2D(ImageProcessor ip) {
        int w = ip.getWidth(), h = ip.getHeight();
        float[] u = new float[w * h];
        for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
                u[y * w + x] = ip.getf(x, y);

        float[] uNew = new float[w * h];
        float kSqInv = (float) (1.0 / (K * K));
        float fdt = (float) dt;

        for (int iter = 0; iter < iterations; iter++) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    int idx = y * w + x;
                    float uc = u[idx];
                    float dE = (x + 1 < w) ? u[idx + 1] - uc : 0;
                    float dW = (x - 1 >= 0) ? u[idx - 1] - uc : 0;
                    float dS = (y + 1 < h) ? u[idx + w] - uc : 0;
                    float dN = (y - 1 >= 0) ? u[idx - w] - uc : 0;
                    float gE, gW, gS, gN;
                    if (mode == 1) {
                        gE = 1.0f / (1.0f + dE * dE * kSqInv);
                        gW = 1.0f / (1.0f + dW * dW * kSqInv);
                        gS = 1.0f / (1.0f + dS * dS * kSqInv);
                        gN = 1.0f / (1.0f + dN * dN * kSqInv);
                    } else {
                        gE = (float) Math.exp(-dE * dE * kSqInv);
                        gW = (float) Math.exp(-dW * dW * kSqInv);
                        gS = (float) Math.exp(-dS * dS * kSqInv);
                        gN = (float) Math.exp(-dN * dN * kSqInv);
                    }
                    uNew[idx] = uc + fdt * (gE * dE + gW * dW + gS * dS + gN * dN);
                }
            }
            System.arraycopy(uNew, 0, u, 0, u.length);
        }

        Object pixels = ip.getPixels();
        if (pixels instanceof byte[]) {
            byte[] out = (byte[]) pixels;
            for (int i = 0; i < out.length; i++) {
                int v = Math.round(u[i]); if (v < 0) v = 0; if (v > 255) v = 255;
                out[i] = (byte) v;
            }
        } else if (pixels instanceof short[]) {
            short[] out = (short[]) pixels;
            for (int i = 0; i < out.length; i++) {
                int v = Math.round(u[i]); if (v < 0) v = 0; if (v > 65535) v = 65535;
                out[i] = (short) v;
            }
        } else if (pixels instanceof float[]) {
            System.arraycopy(u, 0, (float[]) pixels, 0, u.length);
        }
    }
}