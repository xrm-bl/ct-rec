/*
 * NLM_Filter_2D.java - 2D Non-Local Means Filter ImageJ Plugin with Preview
 *
 * Compile: javac -cp ij.jar NLM_Filter_2D.java
 * Install: Place NLM_Filter_2D.class in ImageJ/plugins/
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

public class NLM_Filter_2D implements ExtendedPlugInFilter, DialogListener {

    private static final int FLAGS = DOES_8G | DOES_16 | DOES_32 | KEEP_PREVIEW
            | PARALLELIZE_STACKS | FINAL_PROCESSING;

    private ImagePlus imp;
    private int nPasses = 1;
    private int patchRadius = 1;
    private int searchRadius = 3;
    private double h = 200.0;

    @Override
    public int setup(String arg, ImagePlus imp) {
        if ("final".equals(arg)) return DONE;
        this.imp = imp;
        if (imp == null) {
            IJ.error("NLM Filter 2D", "No image is open");
            return DONE;
        }
        if (h <= 0) h = estimateDefaultH(imp);
        return FLAGS;
    }

    @Override
    public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
        if (h <= 0) h = estimateDefaultH(imp);

        GenericDialog gd = new NonBlockingGenericDialog("2D Non-Local Means Filter");
        gd.addNumericField("Patch radius (1-5):", patchRadius, 0);
        gd.addNumericField("Search radius (1-15):", searchRadius, 0);
        gd.addNumericField("h (filtering strength):", h, 1);
        gd.addMessage("Smaller h = sharper, Larger h = smoother\n\n"
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
        patchRadius = (int) gd.getNextNumber();
        searchRadius = (int) gd.getNextNumber();
        h = gd.getNextNumber();
        if (gd.invalidNumber()) return false;
        if (patchRadius < 1 || patchRadius > 5) return false;
        if (searchRadius < 1 || searchRadius > 15) return false;
        if (h <= 0) return false;
        return true;
    }

    @Override
    public void setNPasses(int nPasses) { this.nPasses = nPasses; }

    @Override
    public void run(ImageProcessor ip) {
        applyNLM2D(ip);
    }

    private double estimateDefaultH(ImagePlus imp) {
        ImageProcessor ip = imp.getProcessor();
        int width = ip.getWidth(), height = ip.getHeight();
        int sx = width / 4, ex = 3 * width / 4;
        int sy = height / 4, ey = 3 * height / 4;
        int step = Math.max(1, (ey - sy) / 20);
        double sumDiffSq = 0; long count = 0;

        for (int y = sy; y < ey; y += step) {
            for (int x = sx; x < ex - 1; x++) {
                double diff = ip.getf(x + 1, y) - ip.getf(x, y);
                sumDiffSq += diff * diff;
                count++;
            }
        }
        if (count < 2) return 50.0;
        double sigma = Math.sqrt(sumDiffSq / count / 2.0);
        return 1.2 * sigma;
    }

    private void applyNLM2D(ImageProcessor ip) {
        int width = ip.getWidth(), height = ip.getHeight();
        Object pixels = ip.getPixels();

        java.awt.Rectangle roi = ip.getRoi();
        final int rx = (roi != null) ? roi.x : 0;
        final int ry = (roi != null) ? roi.y : 0;
        final int rw = (roi != null) ? roi.width : width;
        final int rh = (roi != null) ? roi.height : height;

        if (pixels instanceof byte[]) {
            byte[] in = (byte[]) pixels;
            byte[] out = new byte[in.length];
            System.arraycopy(in, 0, out, 0, in.length);
            nlm8bit(in, out, width, height, rx, ry, rw, rh);
            System.arraycopy(out, 0, in, 0, in.length);
        } else if (pixels instanceof short[]) {
            short[] in = (short[]) pixels;
            short[] out = new short[in.length];
            System.arraycopy(in, 0, out, 0, in.length);
            nlm16bit(in, out, width, height, rx, ry, rw, rh);
            System.arraycopy(out, 0, in, 0, in.length);
        } else if (pixels instanceof float[]) {
            float[] in = (float[]) pixels;
            float[] out = new float[in.length];
            System.arraycopy(in, 0, out, 0, in.length);
            nlm32bit(in, out, width, height, rx, ry, rw, rh);
            System.arraycopy(out, 0, in, 0, in.length);
        }
    }

    private void nlm8bit(byte[] in, byte[] out, int w, int h,
            int roiX, int roiY, int roiW, int roiH) {
        final int pr = patchRadius, sr = searchRadius;
        final double hSqInv = 1.0 / (this.h * this.h);
        final int fw = w, fh = h;
        final byte[] fi = in; final byte[] fo = out;

        int numThreads = Runtime.getRuntime().availableProcessors();
        Thread[] threads = new Thread[numThreads];
        for (int t = 0; t < numThreads; t++) {
            final int tid = t, tt = numThreads;
            threads[t] = new Thread(new Runnable() {
                @Override
                public void run() {
                    for (int y = roiY + tid; y < roiY + roiH; y += tt) {
                        if (y < 0 || y >= fh) continue;
                        for (int x = roiX; x < roiX + roiW; x++) {
                            if (x < 0 || x >= fw) continue;
                            double wSum = 0, vSum = 0, maxW = 0;
                            for (int sy = -sr; sy <= sr; sy++) {
                                int ny = y + sy; if (ny < 0 || ny >= fh) continue;
                                for (int sx = -sr; sx <= sr; sx++) {
                                    int nx = x + sx; if (nx < 0 || nx >= fw) continue;
                                    if (sx == 0 && sy == 0) continue;
                                    double ssd = 0; int cnt = 0;
                                    for (int py = -pr; py <= pr; py++) {
                                        int cy = y + py, qy = ny + py;
                                        if (cy < 0 || cy >= fh || qy < 0 || qy >= fh) continue;
                                        for (int px = -pr; px <= pr; px++) {
                                            int cx = x + px, qx = nx + px;
                                            if (cx < 0 || cx >= fw || qx < 0 || qx >= fw) continue;
                                            double d = (double)(fi[cy*fw+cx]&0xFF) - (double)(fi[qy*fw+qx]&0xFF);
                                            ssd += d * d; cnt++;
                                        }
                                    }
                                    double norm = (cnt > 0) ? ssd / cnt : 0;
                                    double wt = Math.exp(-norm * hSqInv);
                                    if (wt > maxW) maxW = wt;
                                    vSum += wt * (double)(fi[ny*fw+nx]&0xFF);
                                    wSum += wt;
                                }
                            }
                            vSum += maxW * (double)(fi[y*fw+x]&0xFF);
                            wSum += maxW;
                            int r = (wSum > 0) ? (int) Math.round(vSum / wSum) : (fi[y*fw+x]&0xFF);
                            if (r < 0) r = 0; if (r > 255) r = 255;
                            fo[y*fw+x] = (byte) r;
                        }
                    }
                }
            });
            threads[t].start();
        }
        for (Thread th : threads) { try { th.join(); } catch (InterruptedException e) { Thread.currentThread().interrupt(); } }
    }

    private void nlm16bit(short[] in, short[] out, int w, int h,
            int roiX, int roiY, int roiW, int roiH) {
        final int pr = patchRadius, sr = searchRadius;
        final double hSqInv = 1.0 / (this.h * this.h);
        final int fw = w, fh = h;
        final short[] fi = in; final short[] fo = out;

        int numThreads = Runtime.getRuntime().availableProcessors();
        Thread[] threads = new Thread[numThreads];
        for (int t = 0; t < numThreads; t++) {
            final int tid = t, tt = numThreads;
            threads[t] = new Thread(new Runnable() {
                @Override
                public void run() {
                    for (int y = roiY + tid; y < roiY + roiH; y += tt) {
                        if (y < 0 || y >= fh) continue;
                        for (int x = roiX; x < roiX + roiW; x++) {
                            if (x < 0 || x >= fw) continue;
                            double wSum = 0, vSum = 0, maxW = 0;
                            for (int sy = -sr; sy <= sr; sy++) {
                                int ny = y + sy; if (ny < 0 || ny >= fh) continue;
                                for (int sx = -sr; sx <= sr; sx++) {
                                    int nx = x + sx; if (nx < 0 || nx >= fw) continue;
                                    if (sx == 0 && sy == 0) continue;
                                    double ssd = 0; int cnt = 0;
                                    for (int py = -pr; py <= pr; py++) {
                                        int cy = y + py, qy = ny + py;
                                        if (cy < 0 || cy >= fh || qy < 0 || qy >= fh) continue;
                                        for (int px = -pr; px <= pr; px++) {
                                            int cx = x + px, qx = nx + px;
                                            if (cx < 0 || cx >= fw || qx < 0 || qx >= fw) continue;
                                            double d = (double)(fi[cy*fw+cx]&0xFFFF) - (double)(fi[qy*fw+qx]&0xFFFF);
                                            ssd += d * d; cnt++;
                                        }
                                    }
                                    double norm = (cnt > 0) ? ssd / cnt : 0;
                                    double wt = Math.exp(-norm * hSqInv);
                                    if (wt > maxW) maxW = wt;
                                    vSum += wt * (double)(fi[ny*fw+nx]&0xFFFF);
                                    wSum += wt;
                                }
                            }
                            vSum += maxW * (double)(fi[y*fw+x]&0xFFFF);
                            wSum += maxW;
                            int r = (wSum > 0) ? (int) Math.round(vSum / wSum) : (fi[y*fw+x]&0xFFFF);
                            if (r < 0) r = 0; if (r > 65535) r = 65535;
                            fo[y*fw+x] = (short) r;
                        }
                    }
                }
            });
            threads[t].start();
        }
        for (Thread th : threads) { try { th.join(); } catch (InterruptedException e) { Thread.currentThread().interrupt(); } }
    }

    private void nlm32bit(float[] in, float[] out, int w, int h,
            int roiX, int roiY, int roiW, int roiH) {
        final int pr = patchRadius, sr = searchRadius;
        final double hSqInv = 1.0 / (this.h * this.h);
        final int fw = w, fh = h;
        final float[] fi = in; final float[] fo = out;

        int numThreads = Runtime.getRuntime().availableProcessors();
        Thread[] threads = new Thread[numThreads];
        for (int t = 0; t < numThreads; t++) {
            final int tid = t, tt = numThreads;
            threads[t] = new Thread(new Runnable() {
                @Override
                public void run() {
                    for (int y = roiY + tid; y < roiY + roiH; y += tt) {
                        if (y < 0 || y >= fh) continue;
                        for (int x = roiX; x < roiX + roiW; x++) {
                            if (x < 0 || x >= fw) continue;
                            double wSum = 0, vSum = 0, maxW = 0;
                            for (int sy = -sr; sy <= sr; sy++) {
                                int ny = y + sy; if (ny < 0 || ny >= fh) continue;
                                for (int sx = -sr; sx <= sr; sx++) {
                                    int nx = x + sx; if (nx < 0 || nx >= fw) continue;
                                    if (sx == 0 && sy == 0) continue;
                                    double ssd = 0; int cnt = 0;
                                    for (int py = -pr; py <= pr; py++) {
                                        int cy = y + py, qy = ny + py;
                                        if (cy < 0 || cy >= fh || qy < 0 || qy >= fh) continue;
                                        for (int px = -pr; px <= pr; px++) {
                                            int cx = x + px, qx = nx + px;
                                            if (cx < 0 || cx >= fw || qx < 0 || qx >= fw) continue;
                                            double d = (double)fi[cy*fw+cx] - (double)fi[qy*fw+qx];
                                            ssd += d * d; cnt++;
                                        }
                                    }
                                    double norm = (cnt > 0) ? ssd / cnt : 0;
                                    double wt = Math.exp(-norm * hSqInv);
                                    if (wt > maxW) maxW = wt;
                                    vSum += wt * (double)fi[ny*fw+nx];
                                    wSum += wt;
                                }
                            }
                            vSum += maxW * (double)fi[y*fw+x];
                            wSum += maxW;
                            fo[y*fw+x] = (wSum > 0) ? (float)(vSum / wSum) : fi[y*fw+x];
                        }
                    }
                }
            });
            threads[t].start();
        }
        for (Thread th : threads) { try { th.join(); } catch (InterruptedException e) { Thread.currentThread().interrupt(); } }
    }
}