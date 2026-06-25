/*
 * Anisotropic_Diffusion_2D.java - 2D Perona-Malik Anisotropic Diffusion
 * ImageJ Plugin with Preview
 *
 * Compile: javac -cp ij.jar Anisotropic_Diffusion_2D.java
 * Install: Place Anisotropic_Diffusion_2D.class in ImageJ/plugins/
 */

import ij.IJ;
import ij.ImagePlus;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.filter.ExtendedPlugInFilter;
import ij.plugin.filter.PlugInFilterRunner;
import ij.process.ImageProcessor;

import java.awt.AWTEvent;

public class Anisotropic_Diffusion_2D implements ExtendedPlugInFilter, DialogListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;

    private ImagePlus imp;
    private int nPasses = 1;
    private int iterations = 20;
    private double K = -1.0;
    private double dt = 0.2;
    private int mode = 1;

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

    @Override
    public void run(ImageProcessor ip) {
        applyAD2D(ip);
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