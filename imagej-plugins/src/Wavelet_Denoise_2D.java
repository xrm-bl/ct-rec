/*
 * Wavelet_Denoise_2D.java - 2D Wavelet Denoising ImageJ Plugin with Preview
 *
 * Algorithm: 2D Haar DWT + BayesShrink soft thresholding
 *
 * Compile: javac -cp ij.jar Wavelet_Denoise_2D.java
 * Install: Place Wavelet_Denoise_2D.class in ImageJ/plugins/
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

public class Wavelet_Denoise_2D implements ExtendedPlugInFilter, DialogListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;

    private ImagePlus imp;
    private int nPasses = 1;
    private int levels = 3;
    private double thresholdScale = 1.0;

    @Override
    public int setup(String arg, ImagePlus imp) {
        if ("final".equals(arg)) return DONE;
        this.imp = imp;
        if (imp == null) { IJ.error("Wavelet Denoise 2D", "No image is open"); return DONE; }
        return FLAGS;
    }

    @Override
    public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
        GenericDialog gd = new GenericDialog("2D Wavelet Denoising (Haar + BayesShrink)");
        gd.addNumericField("Decomposition levels (1-5):", levels, 0);
        gd.addNumericField("Threshold scale:", thresholdScale, 2);
        gd.addMessage(">1.0 = stronger denoising, <1.0 = less denoising\n\n"
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
        levels = (int) gd.getNextNumber();
        thresholdScale = gd.getNextNumber();
        if (gd.invalidNumber()) return false;
        if (levels < 1 || levels > 5) return false;
        if (thresholdScale <= 0) return false;
        return true;
    }

    @Override
    public void setNPasses(int nPasses) { this.nPasses = nPasses; }

    @Override
    public void run(ImageProcessor ip) {
        applyWavelet2D(ip);
    }

    private static final double SQRT2_INV = 1.0 / Math.sqrt(2.0);

    private void haarForwardX(float[] data, float[] temp, int w, int h, int len) {
        int half = len / 2;
        for (int y = 0; y < h; y++) {
            int base = y * w;
            for (int i = 0; i < half; i++) {
                float a = data[base + 2 * i];
                float b = (2 * i + 1 < len) ? data[base + 2 * i + 1] : a;
                temp[base + i] = (float) ((a + b) * SQRT2_INV);
                temp[base + half + i] = (float) ((a - b) * SQRT2_INV);
            }
            for (int i = 0; i < len; i++) data[base + i] = temp[base + i];
        }
    }

    private void haarForwardY(float[] data, float[] temp, int w, int h, int len) {
        int half = len / 2;
        for (int x = 0; x < w; x++) {
            for (int i = 0; i < half; i++) {
                float a = data[(2 * i) * w + x];
                float b = (2 * i + 1 < len) ? data[(2 * i + 1) * w + x] : a;
                temp[i * w + x] = (float) ((a + b) * SQRT2_INV);
                temp[(half + i) * w + x] = (float) ((a - b) * SQRT2_INV);
            }
            for (int i = 0; i < len; i++) data[i * w + x] = temp[i * w + x];
        }
    }

    private void haarInverseX(float[] data, float[] temp, int w, int h, int len) {
        int half = len / 2;
        for (int y = 0; y < h; y++) {
            int base = y * w;
            for (int i = 0; i < half; i++) {
                float lo = data[base + i], hi = data[base + half + i];
                temp[base + 2 * i] = (float) ((lo + hi) * SQRT2_INV);
                if (2 * i + 1 < len)
                    temp[base + 2 * i + 1] = (float) ((lo - hi) * SQRT2_INV);
            }
            for (int i = 0; i < len; i++) data[base + i] = temp[base + i];
        }
    }

    private void haarInverseY(float[] data, float[] temp, int w, int h, int len) {
        int half = len / 2;
        for (int x = 0; x < w; x++) {
            for (int i = 0; i < half; i++) {
                float lo = data[i * w + x], hi = data[(half + i) * w + x];
                temp[(2 * i) * w + x] = (float) ((lo + hi) * SQRT2_INV);
                if (2 * i + 1 < len)
                    temp[(2 * i + 1) * w + x] = (float) ((lo - hi) * SQRT2_INV);
            }
            for (int i = 0; i < len; i++) data[i * w + x] = temp[i * w + x];
        }
    }

    private float estimateNoiseSigma(float[] data, int x0, int y0, int x1, int y1, int w) {
        double sum = 0; int count = 0;
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                sum += Math.abs(data[y * w + x]);
                count++;
            }
        }
        if (count == 0) return 1.0f;
        return (float) (sum / count * 1.2533141373); /* mean|X| * sqrt(pi/2) */
    }

    private float bayesThreshold(float[] data, int x0, int y0, int x1, int y1, int w, float sigmaN) {
        double sqSum = 0; int count = 0;
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                float v = data[y * w + x];
                sqSum += v * v; count++;
            }
        }
        if (count == 0) return 0.0f;
        double sigmaYsq = sqSum / count;
        double sigmaNsq = sigmaN * sigmaN;
        double sigmaXsq = sigmaYsq - sigmaNsq;
        if (sigmaXsq <= 0) return (float) (Math.sqrt(sigmaYsq) * 10.0);
        return (float) (sigmaNsq / Math.sqrt(sigmaXsq));
    }

    private void softThreshold(float[] data, int x0, int y0, int x1, int y1, int w, float t) {
        for (int y = y0; y < y1; y++) {
            for (int x = x0; x < x1; x++) {
                int idx = y * w + x;
                float v = data[idx];
                if (v > t) data[idx] = v - t;
                else if (v < -t) data[idx] = v + t;
                else data[idx] = 0;
            }
        }
    }

    private void applyWavelet2D(ImageProcessor ip) {
        int w = ip.getWidth(), h = ip.getHeight();
        float[] data = new float[w * h];
        for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
                data[y * w + x] = ip.getf(x, y);

        float[] temp = new float[w * h];

        /* Forward 2D Haar DWT */
        int lw = w, lh = h;
        int actualLevels = 0;
        for (int lv = 0; lv < levels; lv++) {
            if (lw < 2 || lh < 2) break;
            haarForwardX(data, temp, w, h, lw);
            haarForwardY(data, temp, w, h, lh);
            lw /= 2; lh /= 2;
            actualLevels++;
        }

        /* Estimate noise from HH subband at finest level */
        int hw0 = w / 2, hh0 = h / 2;
        float sigmaN = estimateNoiseSigma(data, hw0, hh0, w, h, w);

        /* BayesShrink on each detail subband */
        int cw = w, ch = h;
        for (int lv = 0; lv < actualLevels; lv++) {
            int halfW = cw / 2, halfH = ch / 2;
            /* LH */ float tLH = bayesThreshold(data, 0, halfH, halfW, ch, w, sigmaN) * (float) thresholdScale;
            softThreshold(data, 0, halfH, halfW, ch, w, tLH);
            /* HL */ float tHL = bayesThreshold(data, halfW, 0, cw, halfH, w, sigmaN) * (float) thresholdScale;
            softThreshold(data, halfW, 0, cw, halfH, w, tHL);
            /* HH */ float tHH = bayesThreshold(data, halfW, halfH, cw, ch, w, sigmaN) * (float) thresholdScale;
            softThreshold(data, halfW, halfH, cw, ch, w, tHH);
            cw = halfW; ch = halfH;
        }

        /* Inverse 2D Haar DWT */
        int[] sizesW = new int[actualLevels + 1], sizesH = new int[actualLevels + 1];
        sizesW[0] = w; sizesH[0] = h;
        for (int lv = 1; lv <= actualLevels; lv++) { sizesW[lv] = sizesW[lv-1]/2; sizesH[lv] = sizesH[lv-1]/2; }
        for (int lv = actualLevels - 1; lv >= 0; lv--) {
            haarInverseY(data, temp, w, h, sizesH[lv]);
            haarInverseX(data, temp, w, h, sizesW[lv]);
        }

        /* Write back */
        Object pixels = ip.getPixels();
        if (pixels instanceof byte[]) {
            byte[] out = (byte[]) pixels;
            for (int i = 0; i < out.length; i++) {
                int v = Math.round(data[i]); if (v < 0) v = 0; if (v > 255) v = 255;
                out[i] = (byte) v;
            }
        } else if (pixels instanceof short[]) {
            short[] out = (short[]) pixels;
            for (int i = 0; i < out.length; i++) {
                int v = Math.round(data[i]); if (v < 0) v = 0; if (v > 65535) v = 65535;
                out[i] = (short) v;
            }
        } else if (pixels instanceof float[]) {
            System.arraycopy(data, 0, (float[]) pixels, 0, data.length);
        }
    }
}