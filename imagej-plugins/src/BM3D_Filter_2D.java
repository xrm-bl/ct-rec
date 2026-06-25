/*
 * BM3D_Filter_2D.java - Simplified 2D BM3D Filter ImageJ Plugin with Preview
 *
 * Algorithm: Block-matching with collaborative weighted averaging (simplified BM3D)
 *
 * Compile: javac -cp ij.jar BM3D_Filter_2D.java
 * Install: Place BM3D_Filter_2D.class in ImageJ/plugins/
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

public class BM3D_Filter_2D implements ExtendedPlugInFilter, DialogListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;
    private static final int MAX_MATCHES = 16;

    private ImagePlus imp;
    private int nPasses = 1;
    private int blockRadius = 2;
    private int searchRadius = 3;
    private double sigma = -1.0;

    @Override
    public int setup(String arg, ImagePlus imp) {
        if ("final".equals(arg)) return DONE;
        this.imp = imp;
        if (imp == null) { IJ.error("BM3D Filter 2D", "No image is open"); return DONE; }
        if (sigma <= 0) sigma = estimateSigma(imp);
        return FLAGS;
    }

    @Override
    public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
        if (sigma <= 0) sigma = estimateSigma(imp);

        GenericDialog gd = new NonBlockingGenericDialog("2D BM3D Filter (Simplified)");
        gd.addNumericField("Block radius (1-4):", blockRadius, 0);
        gd.addNumericField("Search radius (1-10):", searchRadius, 0);
        gd.addNumericField("Sigma (noise std):", sigma, 2);
        gd.addMessage("Sigma is auto-estimated from image noise.\n\n"
                + "Note: 2D version processes each slice independently.\n"
                + "For 3D collaborative filtering use tif_bm4d_g.\n\n"
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
        blockRadius = (int) gd.getNextNumber();
        searchRadius = (int) gd.getNextNumber();
        sigma = gd.getNextNumber();
        if (gd.invalidNumber()) return false;
        if (blockRadius < 1 || blockRadius > 4) return false;
        if (searchRadius < 1 || searchRadius > 10) return false;
        if (sigma <= 0) return false;
        return true;
    }

    @Override
    public void setNPasses(int nPasses) { this.nPasses = nPasses; }

    @Override
    public void run(ImageProcessor ip) {
        applyBM3D2D(ip);
    }

    private double estimateSigma(ImagePlus imp) {
        ImageProcessor ip = imp.getProcessor();
        int w = ip.getWidth(), h = ip.getHeight();
        int sx = w / 4, ex = 3 * w / 4, sy = h / 4, ey = 3 * h / 4;
        int step = Math.max(1, (ey - sy) / 20);
        double sumDiffSq = 0; long count = 0;
        for (int y = sy; y < ey; y += step)
            for (int x = sx; x < ex - 1; x++) {
                double d = ip.getf(x + 1, y) - ip.getf(x, y);
                sumDiffSq += d * d; count++;
            }
        if (count < 2) return 50.0;
        return Math.sqrt(sumDiffSq / count / 2.0);
    }

    private void applyBM3D2D(ImageProcessor ip) {
        int w = ip.getWidth(), h = ip.getHeight();
        float[] input = new float[w * h];
        for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
                input[y * w + x] = ip.getf(x, y);

        float[] numer = new float[w * h];
        float[] denom = new float[w * h];

        java.awt.Rectangle roi = ip.getRoi();
        final int roiX = (roi != null) ? roi.x : 0;
        final int roiY = (roi != null) ? roi.y : 0;
        final int roiW = (roi != null) ? roi.width : w;
        final int roiH = (roi != null) ? roi.height : h;

        final float[] fin = input;
        final float[] fnum = numer;
        final float[] fden = denom;
        final int fw = w, fh = h;
        final int br = blockRadius, sr = searchRadius;
        final float sigmaSq = (float) (sigma * sigma);
        final float matchThresh = 2.7f * sigmaSq;
        final float weightNorm = 2.0f * sigmaSq;

        int numThreads = Runtime.getRuntime().availableProcessors();
        Thread[] threads = new Thread[numThreads];

        for (int t = 0; t < numThreads; t++) {
            final int tid = t, tt = numThreads;
            threads[t] = new Thread(new Runnable() {
                @Override
                public void run() {
                    float[] mDist = new float[MAX_MATCHES];
                    int[] mX = new int[MAX_MATCHES], mY = new int[MAX_MATCHES];
                    /* thread-local accumulators to reduce atomic contention */
                    float[] localNum = new float[fw * fh];
                    float[] localDen = new float[fw * fh];

                    for (int y = roiY + tid; y < roiY + roiH; y += tt) {
                        if (y < 0 || y >= fh) continue;
                        for (int x = roiX; x < roiX + roiW; x++) {
                            if (x < 0 || x >= fw) continue;

                            int nm = 0;
                            mDist[0] = 0; mX[0] = x; mY[0] = y; nm = 1;

                            for (int dy = -sr; dy <= sr; dy++) {
                                int ny = y + dy; if (ny < 0 || ny >= fh) continue;
                                for (int dx = -sr; dx <= sr; dx++) {
                                    int nx = x + dx; if (nx < 0 || nx >= fw) continue;
                                    if (dx == 0 && dy == 0) continue;

                                    float ssd = 0; int cnt = 0;
                                    for (int by = -br; by <= br; by++) {
                                        int ry = y + by, cy = ny + by;
                                        if (ry < 0 || ry >= fh || cy < 0 || cy >= fh) continue;
                                        for (int bx = -br; bx <= br; bx++) {
                                            int rx = x + bx, cx = nx + bx;
                                            if (rx < 0 || rx >= fw || cx < 0 || cx >= fw) continue;
                                            float d = fin[ry * fw + rx] - fin[cy * fw + cx];
                                            ssd += d * d; cnt++;
                                        }
                                    }
                                    float nd = (cnt > 0) ? ssd / cnt : Float.MAX_VALUE;
                                    if (nd > matchThresh) continue;

                                    if (nm < MAX_MATCHES) {
                                        int pos = nm;
                                        while (pos > 0 && mDist[pos - 1] > nd) {
                                            if (pos < MAX_MATCHES) { mDist[pos] = mDist[pos-1]; mX[pos] = mX[pos-1]; mY[pos] = mY[pos-1]; }
                                            pos--;
                                        }
                                        mDist[pos] = nd; mX[pos] = nx; mY[pos] = ny;
                                        nm++;
                                    } else if (nd < mDist[MAX_MATCHES - 1]) {
                                        int pos = MAX_MATCHES - 1;
                                        while (pos > 0 && mDist[pos - 1] > nd) {
                                            mDist[pos] = mDist[pos-1]; mX[pos] = mX[pos-1]; mY[pos] = mY[pos-1];
                                            pos--;
                                        }
                                        mDist[pos] = nd; mX[pos] = nx; mY[pos] = ny;
                                    }
                                }
                            }

                            for (int by = -br; by <= br; by++) {
                                int ry = y + by; if (ry < 0 || ry >= fh) continue;
                                for (int bx = -br; bx <= br; bx++) {
                                    int rx = x + bx; if (rx < 0 || rx >= fw) continue;
                                    int oidx = ry * fw + rx;
                                    float vs = 0, ws = 0;
                                    for (int m = 0; m < nm; m++) {
                                        int cx = mX[m] + bx, cy = mY[m] + by;
                                        if (cx < 0 || cx >= fw || cy < 0 || cy >= fh) continue;
                                        float wt = (float) Math.exp(-mDist[m] / weightNorm);
                                        vs += wt * fin[cy * fw + cx];
                                        ws += wt;
                                    }
                                    if (ws > 0) { localNum[oidx] += vs; localDen[oidx] += ws; }
                                }
                            }
                        }
                    }

                    /* Merge thread-local into shared arrays */
                    synchronized (fnum) {
                        for (int i = 0; i < fw * fh; i++) {
                            fnum[i] += localNum[i];
                            fden[i] += localDen[i];
                        }
                    }
                }
            });
            threads[t].start();
        }
        for (Thread th : threads) { try { th.join(); } catch (InterruptedException e) { Thread.currentThread().interrupt(); } }

        /* Final division */
        float[] result = new float[w * h];
        for (int i = 0; i < w * h; i++) {
            result[i] = (denom[i] > 0) ? numer[i] / denom[i] : input[i];
        }

        /* Write back */
        Object pixels = ip.getPixels();
        if (pixels instanceof byte[]) {
            byte[] out = (byte[]) pixels;
            for (int i = 0; i < out.length; i++) {
                int v = Math.round(result[i]); if (v < 0) v = 0; if (v > 255) v = 255;
                out[i] = (byte) v;
            }
        } else if (pixels instanceof short[]) {
            short[] out = (short[]) pixels;
            for (int i = 0; i < out.length; i++) {
                int v = Math.round(result[i]); if (v < 0) v = 0; if (v > 65535) v = 65535;
                out[i] = (short) v;
            }
        } else if (pixels instanceof float[]) {
            System.arraycopy(result, 0, (float[]) pixels, 0, result.length);
        }
    }
}