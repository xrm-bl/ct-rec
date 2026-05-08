/*
 * Open_HIS_IMG.java
 *
 * ImageJ / Fiji plugin for opening HIS format and IMG (HiPic) format files
 * used at SPring-8 synchrotron radiation facility.
 *
 * HIS format: Hamamatsu Image Sequence format
 *   - 64-byte header per frame (magic "IM")
 *   - Variable-length comment
 *   - type=2: 16-bit unsigned short data
 *   - type=6: 12-bit packed data (3 bytes -> 2 pixels)
 *   - Multiple frames stored sequentially
 *
 * IMG format: HiPic single-image format
 *   - 64-byte header (magic "IM")
 *   - Variable-length comment
 *   - 16-bit unsigned short data
 *   - Single frame
 *
 * Based on his2tif6.c and img_ave.c from xrm-bl/ct-rec
 *
 * Author: Generated for SPring-8 BL47XU
 * No external dependencies other than ij.jar
 */

import ij.*;
import ij.io.*;
import ij.gui.*;
import ij.plugin.*;
import ij.process.*;
import java.io.*;
import java.nio.*;
import java.nio.ByteOrder;

public class Open_HIS_IMG implements PlugIn {

    /* Header size in bytes */
    private static final int HEADER_SIZE = 64;

    /* Data type constants */
    private static final int TYPE_16BIT = 2;
    private static final int TYPE_12BIT = 6;

    /**
     * PlugIn entry point. Called from menu (Plugins > SP8CT > IO > Open HIS/IMG...).
     * Opens the image and displays it.
     */
    @Override
    public void run(String arg) {
        String path = arg;

        if (path == null || path.isEmpty()) {
            /* Interactive mode: show file open dialog */
            OpenDialog od = new OpenDialog("Open HIS/IMG File", null);
            String dir = od.getDirectory();
            String name = od.getFileName();
            if (name == null) return;
            path = dir + name;
        }

        ImagePlus imp = openFile(path);
        if (imp != null) {
            imp.show();
        }
    }

    /**
     * Open a HIS or IMG file and return the ImagePlus without showing it.
     * Called by HandleExtraFileTypes for D&D / File > Open support.
     *
     * @param path  absolute path to the file
     * @return ImagePlus on success, null on failure or cancel
     */
    public ImagePlus openFile(String path) {
        File file = new File(path);
        if (!file.exists()) {
            IJ.error("Open HIS/IMG", "File not found: " + path);
            return null;
        }

        String lowerName = file.getName().toLowerCase();

        try {
            if (lowerName.endsWith(".his")) {
                return readHIS(file);
            } else if (lowerName.endsWith(".img")) {
                return readIMG(file);
            } else {
                /* Try to detect from magic bytes */
                return tryOpenByMagic(file);
            }
        } catch (IOException e) {
            IJ.error("Open HIS/IMG", "Error reading file:\n" + e.getMessage());
            return null;
        }
    }

    /**
     * Try to open file by checking the "IM" magic bytes.
     * If detected, ask user which format to use.
     */
    private ImagePlus tryOpenByMagic(File file) throws IOException {
        RandomAccessFile raf = new RandomAccessFile(file, "r");
        try {
            if (raf.length() < HEADER_SIZE) return null;
            byte[] magic = new byte[2];
            raf.readFully(magic);
            if (magic[0] != 'I' || magic[1] != 'M') return null;
        } finally {
            raf.close();
        }

        /* Magic matches "IM" - ask user */
        GenericDialog gd = new GenericDialog("Select Format");
        gd.addMessage("File has 'IM' header. Select format:");
        String[] formats = {"HIS (multi-frame sequence)", "IMG (single frame HiPic)"};
        gd.addChoice("Format:", formats, formats[0]);
        gd.showDialog();
        if (gd.wasCanceled()) return null;

        if (gd.getNextChoiceIndex() == 0) {
            return readHIS(file);
        } else {
            return readIMG(file);
        }
    }

    /* ================================================================
     * HIS format reader (multi-frame)
     * ================================================================ */

    private ImagePlus readHIS(File file) throws IOException {
        IJ.showStatus("Reading HIS file: " + file.getName());

        RandomAccessFile raf = new RandomAccessFile(file, "r");
        try {
            /* Read the first header to get image dimensions */
            byte[] headerBuf = new byte[HEADER_SIZE];
            raf.readFully(headerBuf);
            ByteBuffer hdr = ByteBuffer.wrap(headerBuf).order(ByteOrder.LITTLE_ENDIAN);

            /* Check magic */
            byte m0 = hdr.get(0);
            byte m1 = hdr.get(1);
            if (m0 != 'I' || m1 != 'M') {
                IJ.error("Open HIS", "Not a valid HIS file (missing 'IM' header).");
                return null;
            }

            int commentLength = hdr.getShort(2) & 0xFFFF;
            int width         = hdr.getShort(4) & 0xFFFF;
            int height        = hdr.getShort(6) & 0xFFFF;
            int type          = hdr.getShort(12) & 0xFFFF;
            int nImage1       = hdr.getShort(14) & 0xFFFF;
            int nImage2       = hdr.getShort(16) & 0xFFFF;
            long nImages      = nImage1 + 65536L * nImage2;

            if (width <= 0 || height <= 0) {
                IJ.error("Open HIS", "Invalid image dimensions: " + width + " x " + height);
                return null;
            }

            String typeStr;
            if (type == TYPE_16BIT) {
                typeStr = "16-bit";
            } else if (type == TYPE_12BIT) {
                typeStr = "12-bit packed";
            } else {
                IJ.error("Open HIS", "Unsupported data type: " + type);
                return null;
            }

            IJ.log("HIS file: " + file.getName());
            IJ.log("  Dimensions: " + width + " x " + height);
            IJ.log("  Data type: " + typeStr + " (type=" + type + ")");
            IJ.log("  Declared frames: " + nImages);

            /* Skip comment of first header */
            if (commentLength > 0) {
                raf.skipBytes(commentLength);
            }

            /* Now read all frames */
            ImageStack stack = new ImageStack(width, height);
            int pixelCount = width * height;
            int frameIndex = 0;

            /* Read the first frame data */
            short[] pixels = readFrameData(raf, type, pixelCount);
            if (pixels != null) {
                frameIndex++;
                stack.addSlice("Frame " + frameIndex, new ShortProcessor(width, height, pixels, null));
                IJ.showProgress(frameIndex, (int) Math.max(nImages, 1));
            }

            /* Read subsequent frames */
            while (raf.getFilePointer() < raf.length()) {
                /* Read next header */
                headerBuf = new byte[HEADER_SIZE];
                try {
                    raf.readFully(headerBuf);
                } catch (EOFException e) {
                    break;
                }

                hdr = ByteBuffer.wrap(headerBuf).order(ByteOrder.LITTLE_ENDIAN);

                /* Check magic for each frame */
                if (hdr.get(0) != 'I' || hdr.get(1) != 'M') {
                    IJ.log("  Warning: Frame header mismatch at frame " + (frameIndex + 1) + ", stopping.");
                    break;
                }

                int cl = hdr.getShort(2) & 0xFFFF;
                int fw = hdr.getShort(4) & 0xFFFF;
                int fh = hdr.getShort(6) & 0xFFFF;
                int ft = hdr.getShort(12) & 0xFFFF;

                /* Skip comment */
                if (cl > 0) {
                    raf.skipBytes(cl);
                }

                /* Read frame data */
                pixels = readFrameData(raf, ft, fw * fh);
                if (pixels == null) break;

                frameIndex++;
                stack.addSlice("Frame " + frameIndex, new ShortProcessor(fw, fh, pixels, null));
                IJ.showProgress(frameIndex, (int) Math.max(nImages, 1));
            }

            IJ.log("  Frames read: " + frameIndex);

            if (stack.getSize() == 0) {
                IJ.error("Open HIS", "No valid frames found in file.");
                return null;
            }

            ImagePlus imp = new ImagePlus(file.getName(), stack);
            imp.setProperty("Info",
                "Format: HIS (Hamamatsu Image Sequence)\n" +
                "Data type: " + typeStr + "\n" +
                "Frames: " + frameIndex + "\n" +
                "Source: " + file.getAbsolutePath());
            IJ.showStatus("HIS file loaded: " + frameIndex + " frames");
            return imp;

        } finally {
            raf.close();
        }
    }

    /**
     * Read one frame of pixel data from current file position.
     * Returns 16-bit pixel array (short[]) or null on error/EOF.
     */
    private short[] readFrameData(RandomAccessFile raf, int type, int pixelCount)
            throws IOException {

        if (type == TYPE_16BIT) {
            /* 16-bit unsigned short: 2 bytes per pixel, little-endian */
            int byteCount = pixelCount * 2;
            byte[] raw = new byte[byteCount];
            try {
                raf.readFully(raw);
            } catch (EOFException e) {
                return null;
            }

            short[] pixels = new short[pixelCount];
            ByteBuffer bb = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
            for (int i = 0; i < pixelCount; i++) {
                pixels[i] = bb.getShort();
            }
            return pixels;

        } else if (type == TYPE_12BIT) {
            /* 12-bit packed: 3 bytes -> 2 pixels
             * Packing (from his2tif6.c):
             *   byte[0]*16 + byte[1]/16           -> pixel A
             *   (byte[1]%16)*256 + byte[2]        -> pixel B
             */
            int packedBytes = (pixelCount * 3) / 2;
            byte[] raw = new byte[packedBytes];
            try {
                raf.readFully(raw);
            } catch (EOFException e) {
                return null;
            }

            short[] pixels = new short[pixelCount];
            int pixIdx = 0;

            for (int i = 0; i < packedBytes; i += 3) {
                if (i + 2 >= packedBytes) break;

                int b0 = raw[i]     & 0xFF;
                int b1 = raw[i + 1] & 0xFF;
                int b2 = raw[i + 2] & 0xFF;

                /* Pixel A: b0*16 + b1/16 */
                int valA = b0 * 16 + b1 / 16;
                /* Pixel B: (b1%16)*256 + b2 */
                int valB = (b1 % 16) * 256 + b2;

                if (pixIdx < pixelCount) pixels[pixIdx++] = (short) valA;
                if (pixIdx < pixelCount) pixels[pixIdx++] = (short) valB;
            }
            return pixels;

        } else {
            IJ.log("  Unsupported frame type: " + type);
            return null;
        }
    }

    /* ================================================================
     * IMG format reader (single frame, HiPic format)
     * ================================================================ */

    private ImagePlus readIMG(File file) throws IOException {
        IJ.showStatus("Reading IMG file: " + file.getName());

        RandomAccessFile raf = new RandomAccessFile(file, "r");
        try {
            byte[] headerBuf = new byte[HEADER_SIZE];
            raf.readFully(headerBuf);
            ByteBuffer hdr = ByteBuffer.wrap(headerBuf).order(ByteOrder.LITTLE_ENDIAN);

            /* Check magic */
            byte m0 = hdr.get(0);
            byte m1 = hdr.get(1);
            if (m0 != 'I' || m1 != 'M') {
                IJ.error("Open IMG", "Not a valid IMG file (missing 'IM' header).");
                return null;
            }

            int commentLength = hdr.getShort(2) & 0xFFFF;
            int width         = hdr.getShort(4) & 0xFFFF;
            int height        = hdr.getShort(6) & 0xFFFF;
            int xOffset       = hdr.getShort(8) & 0xFFFF;
            int yOffset       = hdr.getShort(10) & 0xFFFF;
            int type          = hdr.getShort(12) & 0xFFFF;

            if (width <= 0 || height <= 0) {
                IJ.error("Open IMG", "Invalid image dimensions: " + width + " x " + height);
                return null;
            }

            /* Read comment */
            String comment = "";
            if (commentLength > 0) {
                byte[] commentBuf = new byte[commentLength];
                raf.readFully(commentBuf);
                comment = new String(commentBuf, "ISO-8859-1").trim();
            }

            IJ.log("IMG file: " + file.getName());
            IJ.log("  Dimensions: " + width + " x " + height);
            IJ.log("  Offset: (" + xOffset + ", " + yOffset + ")");
            IJ.log("  Type: " + type);
            if (!comment.isEmpty()) {
                IJ.log("  Comment: " + comment);
            }

            /* Read pixel data - IMG is always 16-bit unsigned short */
            int pixelCount = width * height;
            short[] pixels = readFrameData(raf, TYPE_16BIT, pixelCount);
            if (pixels == null) {
                IJ.error("Open IMG", "Error reading pixel data.");
                return null;
            }

            ShortProcessor sp = new ShortProcessor(width, height, pixels, null);
            ImagePlus imp = new ImagePlus(file.getName(), sp);
            imp.setProperty("Info",
                "Format: IMG (HiPic)\n" +
                "Type: " + type + "\n" +
                "Offset: (" + xOffset + ", " + yOffset + ")\n" +
                "Comment: " + comment + "\n" +
                "Source: " + file.getAbsolutePath());
            IJ.showStatus("IMG file loaded: " + width + " x " + height);
            return imp;

        } finally {
            raf.close();
        }
    }
}
