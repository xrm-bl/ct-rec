/*
 * Cylinder_Unwrap_3D.java - cylindrical (theta-Z) unwrapping of a CT volume
 *
 * Unrolls a cylindrical sample whose axis is parallel to the CT rotation axis
 * (i.e. the stack's Z axis) into a stack of theta-Z maps, one per radius:
 *
 *     output X = theta   (Ntheta samples over the angular span)
 *     output Y = Z       (input slice index, first..last)
 *     output slice = R   (Rmin, Rmin+dR, ... Rmax)
 *
 * The sampling point in the input volume is
 *
 *     x = Xc + R*cos(theta),  y = Yc + R*sin(theta),  z = z
 *
 * so theta = 0 points along +X. Theta increases toward +Y in array coordinates,
 * which looks clockwise on screen because ImageJ's Y axis points down. Use the
 * "Reverse theta direction" option to flip it.
 *
 * 8-bit, 16-bit unsigned and 32-bit float grayscale stacks are supported; the
 * output keeps the input's bit depth.
 *
 * Provenance: every run records its parameters, including a ready-to-run macro
 * call, into the output's "Info" property (Image > Show Info) and the ImageJ Log.
 * ImageJ keeps that property and the per-slice R labels inside TIFF stacks it
 * writes, so the record travels with the file. Saving as an image sequence drops
 * it - tick "Write parameter file" to get the same record as a .txt as well.
 *
 * From a macro (the labels' first words are the option keys):
 *   run("Cylinder Unwrap 3D", "x=1023.5 y=1024.25 size=Radius outer=800 inner=200 "
 *       + "radial=1 theta=5027 start=0 angular=360 first=1 last=1500 "
 *       + "interpolation=Bilinear background=0");
 * Optional flags: "reverse" (reverse theta), "write path=[C:/dir/params.txt]".
 *
 * Compile:
 *   javac -cp ij.jar Cylinder_Unwrap_3D.java
 *
 * Install:
 *   Place Cylinder_Unwrap_3D.class in ImageJ/plugins/
 */

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.Macro;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.io.SaveDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import ij.process.ShortProcessor;

import java.awt.Rectangle;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

public class Cylinder_Unwrap_3D implements PlugIn {

    private static final String TITLE = "Cylinder Unwrap 3D";

    /* Remembered between invocations. NaN / <=0 means "derive from the image". */
    private static double xc = Double.NaN, yc = Double.NaN;
    private static boolean useDiameter = false;
    private static double outerSize = Double.NaN;   // radius or diameter, see useDiameter
    private static double rMin = 0.0;
    private static double dR = 1.0;
    private static int nThetaPref = 0;              // 0 = auto: round(span/360 * 2*PI*Rmax)
    private static double thetaStart = 0.0;         // deg
    private static double thetaSpan = 360.0;        // deg
    private static boolean reverseTheta = false;
    private static boolean bilinear = true;
    private static double bgValue = 0.0;
    private static boolean writeParamFile = false;

    /* Not remembered - always defaults to the whole stack. */
    private int zFirst, zLast;

    private volatile boolean aborted = false;

    public void run(String arg) {
        ImagePlus imp = IJ.getImage();          // reports "no image" itself
        if (imp == null) return;

        int bd = imp.getBitDepth();
        if (bd != 8 && bd != 16 && bd != 32) {
            IJ.error(TITLE, "8-bit, 16-bit or 32-bit grayscale required.\n"
                    + "RGB images are not supported.");
            return;
        }
        if (imp.getNChannels() > 1 || imp.getNFrames() > 1) {
            IJ.error(TITLE, "A single-channel, single-frame stack is required.\n"
                    + "Split the hyperstack first (Image > Hyperstacks > Stack to Hyperstack / Split).");
            return;
        }

        final int w = imp.getWidth();
        final int h = imp.getHeight();
        final int nz = imp.getStackSize();
        if (nz < 2) {
            IJ.error(TITLE, "A stack of at least 2 slices is required.");
            return;
        }

        applyDefaults(imp, w, h);
        zFirst = 1;
        zLast = nz;

        if (!showDialog(imp, w, h, nz, bd)) return;

        process(imp, bd, w, h);
    }

    /* ---------------------------------------------------------------- dialog */

    /** Seeds the centre and the radius from the ROI (if any) or the image size. */
    private void applyDefaults(ImagePlus imp, int w, int h) {
        double defXc = (w - 1) / 2.0;
        double defYc = (h - 1) / 2.0;
        double defR = Math.min(w, h) / 2.0 - 1.0;

        Roi roi = imp.getRoi();
        if (roi != null) {
            Rectangle b = roi.getBounds();
            defXc = b.x + (b.width - 1) / 2.0;
            defYc = b.y + (b.height - 1) / 2.0;
            defR = Math.min(b.width, b.height) / 2.0;
        }

        /* A ROI is an explicit gesture - let it override the remembered values. */
        if (roi != null || Double.isNaN(xc) || Double.isNaN(yc)) {
            xc = defXc;
            yc = defYc;
        }
        if (roi != null || Double.isNaN(outerSize) || outerSize <= 0) {
            outerSize = useDiameter ? 2.0 * defR : defR;
        }
    }

    private boolean showDialog(ImagePlus imp, int w, int h, int nz, int bd) {
        GenericDialog gd = new GenericDialog(TITLE);
        gd.addMessage("Input: " + w + " x " + h + " x " + nz + ", " + bd + "-bit\n"
                + "Coordinates are 0-based pixel indices (as shown in the ImageJ status bar).");
        /* Each label's FIRST WORD becomes its macro key, and a duplicated first
         * word makes every field after the first unreachable from a macro - so
         * keep the leading words distinct:
         *   x y size outer inner radial theta start angular first last
         *   interpolation background reverse */
        gd.addNumericField("X center:", xc, 2, 10, "px");
        gd.addNumericField("Y center:", yc, 2, 10, "px");
        gd.addChoice("Size given as:", new String[] { "Radius", "Diameter" },
                useDiameter ? "Diameter" : "Radius");
        gd.addNumericField("Outer size:", outerSize, 2, 10, "px");
        gd.addNumericField("Inner radius:", rMin, 2, 10, "px");
        gd.addNumericField("Radial step:", dR, 3, 10, "px");
        gd.addMessage("");
        gd.addNumericField("Theta samples:", nThetaPref, 0, 10, "0 = auto (1 px at outer radius)");
        gd.addNumericField("Start angle:", thetaStart, 2, 10, "deg");
        gd.addNumericField("Angular span:", thetaSpan, 2, 10, "deg");
        gd.addCheckbox("Reverse theta direction", reverseTheta);
        gd.addMessage("");
        gd.addNumericField("First slice:", zFirst, 0, 10, "");
        gd.addNumericField("Last slice:", zLast, 0, 10, "");
        gd.addChoice("Interpolation:", new String[] { "Bilinear", "Nearest neighbor" },
                bilinear ? "Bilinear" : "Nearest neighbor");
        gd.addNumericField("Background value:", bgValue, 3, 10, "outside the image");
        gd.addCheckbox("Write parameter file", writeParamFile);
        gd.showDialog();
        if (gd.wasCanceled()) return false;

        xc = gd.getNextNumber();
        yc = gd.getNextNumber();
        useDiameter = gd.getNextChoiceIndex() == 1;
        outerSize = gd.getNextNumber();
        rMin = gd.getNextNumber();
        dR = gd.getNextNumber();
        nThetaPref = (int) gd.getNextNumber();
        thetaStart = gd.getNextNumber();
        thetaSpan = gd.getNextNumber();
        reverseTheta = gd.getNextBoolean();
        zFirst = (int) gd.getNextNumber();
        zLast = (int) gd.getNextNumber();
        bilinear = gd.getNextChoiceIndex() == 0;
        bgValue = gd.getNextNumber();
        writeParamFile = gd.getNextBoolean();

        if (gd.invalidNumber()) {
            IJ.error(TITLE, "Invalid number in the dialog.");
            return false;
        }
        return validate(nz, bd);
    }

    private boolean validate(int nz, int bd) {
        double rMax = rMax();
        if (dR <= 0) {
            IJ.error(TITLE, "Radial step must be positive.");
            return false;
        }
        if (rMin < 0) {
            IJ.error(TITLE, "Inner radius must not be negative.");
            return false;
        }
        if (rMax <= rMin) {
            IJ.error(TITLE, "Outer radius (" + IJ.d2s(rMax, 3) + " px) must be larger than the "
                    + "inner radius (" + IJ.d2s(rMin, 3) + " px).");
            return false;
        }
        if (thetaSpan == 0) {
            IJ.error(TITLE, "Theta span must not be zero.");
            return false;
        }
        if (zFirst < 1 || zLast > nz || zFirst > zLast) {
            IJ.error(TITLE, "Slice range must satisfy 1 <= first <= last <= " + nz + ".");
            return false;
        }
        if (nThetaPref < 0) {
            IJ.error(TITLE, "Theta samples must not be negative.");
            return false;
        }

        long total = (long) nTheta() * (zLast - zFirst + 1) * nR() * bytesPerPixel(bd);
        Runtime rt = Runtime.getRuntime();
        long avail = rt.maxMemory() - (rt.totalMemory() - rt.freeMemory());
        if (total > avail * 8 / 10) {
            IJ.error(TITLE, "The output needs about " + (total >> 20) + " MB but only about "
                    + (avail >> 20) + " MB is available.\n\n"
                    + "Raise the ImageJ heap (Edit > Options > Memory & Threads), or reduce\n"
                    + "the radial/angular/slice range.");
            return false;
        }
        return true;
    }

    /* ------------------------------------------------------------- geometry */

    private double rMax() {
        return useDiameter ? outerSize / 2.0 : outerSize;
    }

    /** Number of radial samples, inclusive of both rMin and the last step <= rMax. */
    private int nR() {
        return (int) Math.floor((rMax() - rMin) / dR + 1e-9) + 1;
    }

    /** Resolved number of angular samples: one pixel of arc at rMax when auto. */
    private int nTheta() {
        if (nThetaPref > 0) return nThetaPref;
        double arc = Math.abs(thetaSpan) / 360.0 * 2.0 * Math.PI * rMax();
        return Math.max(1, (int) Math.round(arc));
    }

    private static int bytesPerPixel(int bd) {
        return bd == 8 ? 1 : bd == 16 ? 2 : 4;
    }

    /* ------------------------------------------------------------ processing */

    private void process(ImagePlus imp, final int bd, final int w, final int h) {
        final int nT = nTheta();
        final int nRad = nR();
        final int nzOut = zLast - zFirst + 1;
        final double rMax = rMax();

        /* Angular sampling grid: dTheta = span / nTheta, so a full 360 deg span
         * wraps without repeating the first column. */
        final double[] cosT = new double[nT];
        final double[] sinT = new double[nT];
        final double dTheta = Math.toRadians(thetaSpan) / nT;
        final double a0 = Math.toRadians(thetaStart);
        final double sgn = reverseTheta ? -1.0 : 1.0;
        for (int i = 0; i < nT; i++) {
            double a = a0 + sgn * i * dTheta;
            cosT[i] = Math.cos(a);
            sinT[i] = Math.sin(a);
        }

        /* One output plane per radius, already in the final pixel type. */
        final Object[] planes = new Object[nRad];
        for (int i = 0; i < nRad; i++) {
            int n = nT * nzOut;
            planes[i] = bd == 8 ? (Object) new byte[n] : bd == 16 ? (Object) new short[n] : (Object) new float[n];
        }

        final ImageStack in = imp.getStack();
        final AtomicInteger done = new AtomicInteger(0);
        final AtomicLong outside = new AtomicLong(0);
        final long t0 = System.currentTimeMillis();

        IJ.showStatus(TITLE + ": " + nT + " x " + nzOut + " x " + nRad + " ...");
        aborted = false;

        /* Each z fills row (z - zFirst) of every output plane, so slices are
         * independent. VirtualStack readers are not safe to share, so those run
         * single-threaded. */
        int nThreads = in.isVirtual() ? 1
                : Math.max(1, Math.min(Runtime.getRuntime().availableProcessors(), nzOut));
        ExecutorService pool = Executors.newFixedThreadPool(nThreads);
        try {
            for (int z = zFirst; z <= zLast; z++) {
                final int zi = z;
                pool.submit(new Runnable() {
                    public void run() {
                        if (aborted) return;
                        try {
                            float[] src = (float[]) in.getProcessor(zi).toFloat(0, null).getPixels();
                            long oob = unwrapSlice(src, w, h, zi - zFirst, nT, nzOut, nRad,
                                    cosT, sinT, planes, bd);
                            if (oob > 0) outside.addAndGet(oob);
                        } catch (RuntimeException e) {
                            aborted = true;
                            IJ.handleException(e);
                            return;
                        }
                        int d = done.incrementAndGet();
                        IJ.showProgress(d, nzOut);
                        if (IJ.escapePressed()) aborted = true;
                    }
                });
            }
        } finally {
            pool.shutdown();
        }
        try {
            pool.awaitTermination(365, TimeUnit.DAYS);
        } catch (InterruptedException e) {
            aborted = true;
        }
        IJ.showProgress(1.0);

        if (aborted) {
            IJ.showStatus(TITLE + ": aborted");
            IJ.error(TITLE, "Aborted before completion.");
            return;
        }

        /* ---- assemble the output stack ---- */
        ImageStack os = new ImageStack(nT, nzOut);
        for (int i = 0; i < nRad; i++) {
            double r = rMin + i * dR;
            ImageProcessor ip;
            if (bd == 8) {
                ip = new ByteProcessor(nT, nzOut, (byte[]) planes[i], null);
            } else if (bd == 16) {
                ip = new ShortProcessor(nT, nzOut, (short[]) planes[i], null);
            } else {
                ip = new FloatProcessor(nT, nzOut, (float[]) planes[i]);
            }
            os.addSlice("R=" + IJ.d2s(r, 3) + " px", ip);
            planes[i] = null;
        }

        ImagePlus out = new ImagePlus(imp.getShortTitle() + "-unwrap", os);
        out.setDisplayRange(imp.getDisplayRangeMin(), imp.getDisplayRangeMax());

        /* The theta axis has a different arc length on every slice, so no single
         * pixel size is correct - leave the output uncalibrated and record the
         * geometry as text instead. ImageJ stores the "Info" property inside the
         * TIFF (private tags 50838/50839) together with the slice labels, so the
         * parameters travel with the file as long as it is saved as a TIFF stack
         * from ImageJ. Saving as an image sequence drops them - hence the
         * optional parameter file. */
        String params = buildParamText(imp, bd, nT, nzOut, nRad, outside.get());
        out.setProperty("Info", params);
        out.show();

        long ms = System.currentTimeMillis() - t0;
        IJ.showStatus(TITLE + ": " + nT + " x " + nzOut + " x " + nRad + " in " + ms + " ms");
        IJ.log(TITLE + ": " + imp.getTitle() + " -> " + out.getTitle() + "  "
                + nT + " x " + nzOut + " x " + nRad + " (" + bd + "-bit), " + ms + " ms");
        IJ.log("  " + macroLine(nT));
        if (outside.get() > 0) {
            IJ.log(TITLE + ": WARNING - " + outside.get() + " sample(s) fell outside the image "
                    + "and were set to " + IJ.d2s(bgValue, 4) + ".");
        }
        if (writeParamFile) saveParamFile(out.getShortTitle(), params);
    }

    /* ------------------------------------------------------------- provenance */

    /**
     * The full parameter record: stored as the output's "Info" property (visible
     * via Image &gt; Show Info, preserved inside ImageJ-written TIFF stacks) and
     * optionally written next to the saved image as a text file.
     */
    private String buildParamText(ImagePlus imp, int bd, int nT, int nzOut, int nRad, long outside) {
        Calibration cal = imp.getCalibration();
        String unit = cal.getUnit();
        boolean calibrated = unit != null && !"pixel".equals(unit) && cal.pixelWidth > 0;

        StringBuilder s = new StringBuilder();
        s.append("Cylinder Unwrap 3D (SP8CT ImageJ plugin)\n");
        s.append("date=").append(new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date())).append("\n");
        s.append("imagej=").append(IJ.getVersion()).append("\n");
        s.append("source=").append(imp.getTitle()).append("\n");
        FileInfo ofi = imp.getOriginalFileInfo();
        if (ofi != null && ofi.directory != null) {
            s.append("source_path=").append(ofi.directory)
                    .append(ofi.fileName == null ? "" : ofi.fileName).append("\n");
        }
        s.append("source_size=").append(imp.getWidth()).append("x").append(imp.getHeight())
                .append("x").append(imp.getStackSize()).append(" ").append(bd).append("-bit\n");
        s.append("\n");
        s.append("center_x=").append(num(xc)).append("\n");
        s.append("center_y=").append(num(yc)).append("\n");
        s.append("r_min=").append(num(rMin)).append("\n");
        s.append("r_max=").append(num(rMax())).append("\n");
        s.append("r_step=").append(num(dR)).append("\n");
        s.append("r_slices=").append(nRad).append("\n");
        s.append("theta_samples=").append(nT).append("\n");
        s.append("theta_auto=").append(nThetaPref <= 0).append("\n");
        s.append("theta_start_deg=").append(num(thetaStart)).append("\n");
        s.append("theta_span_deg=").append(num(thetaSpan)).append("\n");
        s.append("theta_step_deg=").append(num(thetaSpan / nT)).append("\n");
        s.append("theta_reversed=").append(reverseTheta).append("\n");
        s.append("z_first=").append(zFirst).append("\n");
        s.append("z_last=").append(zLast).append("\n");
        s.append("interpolation=").append(bilinear ? "bilinear" : "nearest").append("\n");
        s.append("background=").append(num(bgValue)).append("\n");
        s.append("out_of_range_samples=").append(outside).append("\n");
        s.append("\n");
        s.append("output=").append(nT).append("x").append(nzOut).append("x").append(nRad)
                .append(" ").append(bd).append("-bit\n");
        s.append("axes=X:theta Y:z slice:R\n");
        s.append("theta_0_direction=+X, theta increases toward +Y in array coordinates")
                .append(reverseTheta ? " (reversed here)" : "").append("\n");
        if (calibrated) {
            s.append("source_pixel_xy=").append(num(cal.pixelWidth)).append(" ").append(unit).append("\n");
            s.append("source_pixel_z=").append(num(cal.pixelDepth)).append(" ").append(unit).append("\n");
            s.append("output_y_spacing=").append(num(cal.pixelDepth)).append(" ").append(unit).append("\n");
            s.append("output_x_arc_at_R=R * ").append(num(cal.pixelWidth)).append(" ").append(unit)
                    .append(" * ").append(num(Math.toRadians(thetaSpan / nT))).append(" rad\n");
        }
        s.append("\n");
        s.append("# re-run with:\n");
        s.append(macroLine(nT)).append("\n");
        return s.toString();
    }

    /**
     * An ImageJ macro call reproducing this run. The resolved theta count is
     * recorded rather than the "0 = auto" placeholder so the result is
     * reproducible independently of how auto is derived.
     */
    private String macroLine(int nT) {
        StringBuilder s = new StringBuilder();
        s.append("run(\"Cylinder Unwrap 3D\", \"");
        s.append("x=").append(num(xc));
        s.append(" y=").append(num(yc));
        s.append(" size=").append(useDiameter ? "Diameter" : "Radius");
        s.append(" outer=").append(num(outerSize));
        s.append(" inner=").append(num(rMin));
        s.append(" radial=").append(num(dR));
        s.append(" theta=").append(nT);
        s.append(" start=").append(num(thetaStart));
        s.append(" angular=").append(num(thetaSpan));
        s.append(" first=").append(zFirst);
        s.append(" last=").append(zLast);
        s.append(" interpolation=").append(bilinear ? "Bilinear" : "[Nearest neighbor]");
        s.append(" background=").append(num(bgValue));
        if (reverseTheta) s.append(" reverse");
        s.append("\");");
        return s.toString();
    }

    /**
     * Writes the parameter record as plain UTF-8 text. Interactively the path is
     * asked for; from a macro it must be given as {@code path=[...]} so the run
     * never blocks on a file dialog.
     */
    private void saveParamFile(String defaultName, String params) {
        File f;
        String macroOptions = Macro.getOptions();
        if (macroOptions != null) {
            String path = Macro.getValue(macroOptions, "path", "");
            if (path.length() == 0) {
                IJ.log(TITLE + ": parameter file skipped - add path=[...] to the macro options.");
                return;
            }
            f = new File(path);
        } else {
            SaveDialog sd = new SaveDialog("Save unwrap parameters", defaultName + "-params", ".txt");
            if (sd.getDirectory() == null || sd.getFileName() == null) {
                IJ.log(TITLE + ": parameter file not written (cancelled).");
                return;
            }
            f = new File(sd.getDirectory(), sd.getFileName());
        }
        Writer wr = null;
        try {
            wr = new OutputStreamWriter(new FileOutputStream(f), "UTF-8");
            wr.write(params);
            IJ.log(TITLE + ": parameters written to " + f.getAbsolutePath());
        } catch (IOException e) {
            IJ.error(TITLE, "Could not write the parameter file:\n" + e.getMessage());
        } finally {
            if (wr != null) try { wr.close(); } catch (IOException ignored) { }
        }
    }

    /** Locale-independent compact number formatting, safe to paste into a macro. */
    private static String num(double v) {
        if (v == Math.rint(v) && Math.abs(v) < 1e15) return String.valueOf((long) v);
        String s = IJ.d2s(v, 6);
        while (s.endsWith("0")) s = s.substring(0, s.length() - 1);
        if (s.endsWith(".")) s = s.substring(0, s.length() - 1);
        return s;
    }

    /**
     * Fills row {@code row} of every output plane from one input slice.
     *
     * @return the number of samples that fell outside the input image.
     */
    private long unwrapSlice(float[] src, int w, int h, int row, int nT, int nzOut, int nRad,
                             double[] cosT, double[] sinT, Object[] planes, int bd) {
        final float bg = (float) bgValue;
        final int base = row * nT;
        long oob = 0;

        for (int ri = 0; ri < nRad; ri++) {
            double r = rMin + ri * dR;
            Object plane = planes[ri];
            for (int ti = 0; ti < nT; ti++) {
                double x = xc + r * cosT[ti];
                double y = yc + r * sinT[ti];
                float v;
                if (x < -0.5 || y < -0.5 || x > w - 0.5 || y > h - 0.5) {
                    v = bg;
                    oob++;
                } else if (bilinear) {
                    v = bilinear(src, w, h, x, y);
                } else {
                    int xi = (int) Math.floor(x + 0.5);
                    int yi = (int) Math.floor(y + 0.5);
                    if (xi < 0) xi = 0; else if (xi > w - 1) xi = w - 1;
                    if (yi < 0) yi = 0; else if (yi > h - 1) yi = h - 1;
                    v = src[yi * w + xi];
                }
                int i = base + ti;
                if (bd == 8) {
                    ((byte[]) plane)[i] = (byte) clamp(v, 255);
                } else if (bd == 16) {
                    ((short[]) plane)[i] = (short) clamp(v, 65535);
                } else {
                    ((float[]) plane)[i] = v;
                }
            }
        }
        return oob;
    }

    /** Bilinear sample with edge clamping; x,y are pixel-index coordinates. */
    private static float bilinear(float[] px, int w, int h, double x, double y) {
        int x0 = (int) Math.floor(x);
        int y0 = (int) Math.floor(y);
        double fx = x - x0;
        double fy = y - y0;
        int x1 = x0 + 1;
        int y1 = y0 + 1;
        if (x0 < 0) { x0 = 0; if (x1 < 0) x1 = 0; }
        if (y0 < 0) { y0 = 0; if (y1 < 0) y1 = 0; }
        if (x1 > w - 1) { x1 = w - 1; if (x0 > w - 1) x0 = w - 1; }
        if (y1 > h - 1) { y1 = h - 1; if (y0 > h - 1) y0 = h - 1; }

        int r0 = y0 * w;
        int r1 = y1 * w;
        double v00 = px[r0 + x0];
        double v10 = px[r0 + x1];
        double v01 = px[r1 + x0];
        double v11 = px[r1 + x1];
        double top = v00 + fx * (v10 - v00);
        double bot = v01 + fx * (v11 - v01);
        return (float) (top + fy * (bot - top));
    }

    /** Rounds to the nearest integer and clamps to [0, max]. */
    private static int clamp(float v, int max) {
        if (v != v) return 0;                       // NaN
        int i = (int) Math.floor(v + 0.5);
        return i < 0 ? 0 : i > max ? max : i;
    }
}
