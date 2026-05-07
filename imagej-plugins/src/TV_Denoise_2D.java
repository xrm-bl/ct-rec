/*
 * TV_Denoise_2D.java - 2D Total Variation Denoising ImageJ Plugin with Preview
 *
 * Algorithm: Chambolle-Pock (Primal-Dual) for 2D isotropic TV
 *
 * Compile: javac -cp ij.jar TV_Denoise_2D.java
 * Install: Place TV_Denoise_2D.class in ImageJ/plugins/
 */

import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;

import java.awt.AWTEvent;

public class TV_Denoise_2D implements ExtendedPlugInFilter, DialogListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;

    private ImagePlus imp;
    private int nPasses = 1;
    private double lambda = 10.0;
    private int iterations = 50;

    @Override
    public int setup(String arg, ImagePlus imp) {
        if ("final".equals(arg)) return DONE;
        this.imp = imp;
        if (imp == null) { IJ.error("TV Denoise 2D", "No image is open"); return DONE; }
        return FLAGS;
    }

    @Override
    public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
        GenericDialog gd = new GenericDialog("2D Total Variation Denoising");
        gd.addNumericField("Lambda (larger=less denoising):", lambda, 2);
        gd.addNumericField("Iterations:", iterations, 0);
        gd.addMessage("Lambda: 1-5=strong, 5-20=standard, 20-100=weak\n\n"
                + "Tip: Ctrl+Shift+C for Brightness/Contrast during preview.");

        gd.addPreviewCheckbox(pfr);
        gd.addDialogListener(this);
        gd.showDialog();
        if (gd.wasCanceled()) return DONE;
        IJ.register(this.getClass());
        return IJ.setupDialog(imp, FLAGS);
    }

    @Override
    public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        lambda = gd.getNextNumber();
        iterations = (int) gd.getNextNumber();
        if (gd.invalidNumber()) return false;
        if (lambda <= 0 || iterations < 1 || iterations > 5000) return false;
        return true;
    }

    @Override
    public void setNPasses(int nPasses) { this.nPasses = nPasses; }

    @Override
    public void run(ImageProcessor ip) {
        applyTV2D(ip);
    }

    private void applyTV2D(ImageProcessor ip) {
        int w = ip.getWidth(), h = ip.getHeight();

        /* Convert to float for processing */
        float[] f = new float[w * h];
        for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
                f[y * w + x] = ip.getf(x, y);

        float[] u = f.clone();
        float[] uBar = f.clone();
        float[] px = new float[w * h];
        float[] py = new float[w * h];

        /* Chambolle-Pock step sizes for 2D: L=sqrt(8) */
        float L = (float) Math.sqrt(8.0);
        float tau = 1.0f / L;
        float sigma = 1.0f / L;
        float theta = 1.0f;
        float lam = (float) lambda;

        for (int iter = 0; iter < iterations; iter++) {
            /* Dual update */
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    int idx = y * w + x;
                    float ub = uBar[idx];
                    float gx = (x + 1 < w) ? uBar[idx + 1] - ub : 0;
                    float gy = (y + 1 < h) ? uBar[idx + w] - ub : 0;
                    float npx = px[idx] + sigma * gx;
                    float npy = py[idx] + sigma * gy;
                    float norm = (float) Math.sqrt(npx * npx + npy * npy);
                    if (norm > 1.0f) { npx /= norm; npy /= norm; }
                    px[idx] = npx;
                    py[idx] = npy;
                }
            }
            /* Primal update */
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    int idx = y * w + x;
                    float divPx = (x == 0) ? px[idx] : (x == w - 1) ? -px[idx - 1] : px[idx] - px[idx - 1];
                    float divPy = (y == 0) ? py[idx] : (y == h - 1) ? -py[idx - w] : py[idx] - py[idx - w];
                    float divP = divPx + divPy;
                    float uOld = u[idx];
                    float uNew = (uOld + tau * divP + tau * lam * f[idx]) / (1.0f + tau * lam);
                    u[idx] = uNew;
                    uBar[idx] = uNew + theta * (uNew - uOld);
                }
            }
        }

        /* Write back */
        Object pixels = ip.getPixels();
        if (pixels instanceof byte[]) {
            byte[] out = (byte[]) pixels;
            for (int i = 0; i < out.length; i++) {
                int v = Math.round(u[i]);
                if (v < 0) v = 0; if (v > 255) v = 255;
                out[i] = (byte) v;
            }
        } else if (pixels instanceof short[]) {
            short[] out = (short[]) pixels;
            for (int i = 0; i < out.length; i++) {
                int v = Math.round(u[i]);
                if (v < 0) v = 0; if (v > 65535) v = 65535;
                out[i] = (short) v;
            }
        } else if (pixels instanceof float[]) {
            float[] out = (float[]) pixels;
            System.arraycopy(u, 0, out, 0, out.length);
        }
    }
}