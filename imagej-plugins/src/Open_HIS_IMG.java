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
 *   - Supports Virtual Stack for large files
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
import java.util.ArrayList;

public class Open_HIS_IMG implements PlugIn {

    /* Header size in bytes */
    private static final int HEADER_SIZE = 64;

    /* Data type constants */
    private static final int TYPE_16BIT = 2;
    private static final int TYPE_12BIT = 6;

    /** When true, the interactive dialog was shown and user chose virtual. */
    private boolean useVirtualStack = false;

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

        File file = new File(path);
        if (!file.exists()) {
            IJ.error("Open HIS/IMG", "File not found: " + path);
            return;
        }

        /* Show options dialog for HIS files */
        String lowerName = file.getName().toLowerCase();
        if (lowerName.endsWith(".his")) {
            if (!showHISOptionsDialog(file)) return;  /* cancelled */
        }

        ImagePlus imp = openFileInternal(file);
        if (imp != null) {
            imp.show();
        }
    }

    /**
     * Open a HIS or IMG file and return the ImagePlus without showing it.
     * Called by HandleExtraFileTypes for D&D / File > Open / double-click.
     * For HIS files, shows the options dialog to allow Virtual Stack choice.
     *
     * @param path  absolute path to the file
     * @return ImagePlus on success, null on failure or cancel
     */
    public ImagePlus openFile(String path) {
        useVirtualStack = false;
        File file = new File(path);
        if (!file.exists()) {
            IJ.error("Open HIS/IMG", "File not found: " + path);
            return null;
        }

        /* Show options dialog for HIS files */
        String lowerName = file.getName().toLowerCase();
        if (lowerName.endsWith(".his")) {
            if (!showHISOptionsDialog(file)) return null;  /* cancelled */
        }

        return openFileInternal(file);
    }

    /**
     * Internal file open dispatcher.
     */
    private ImagePlus openFileInternal(File file) {
        String lowerName = file.getName().toLowerCase();

        try {
            if (lowerName.endsWith(".his")) {
                if (useVirtualStack) {
                    return readHISVirtual(file);
                } else {
                    return readHIS(file);
                }
            } else if (lowerName.endsWith(".img")) {
                return readIMG(file);
            } else {
                return tryOpenByMagic(file);
            }
        } catch (IOException e) {
            IJ.error("Open HIS/IMG", "Error reading file:\n" + e.getMessage());
            return null;
        }
    }

    /**
     * Show options dialog for HIS files.
     * Returns false if user cancelled.
     */
    private boolean showHISOptionsDialog(File file) {
        long fileSize = file.length();
        String sizeStr;
        if (fileSize > 1024L * 1024 * 1024) {
            sizeStr = String.format("%.1f GB", fileSize / (1024.0 * 1024 * 1024));
        } else {
            sizeStr = String.format("%.1f MB", fileSize / (1024.0 * 1024));
        }

        /* Scan frame count quickly */
        int nFrames = 0;
        int width = 0, height = 0, type = 0;
        try {
            RandomAccessFile raf = new RandomAccessFile(file, "r");
            try {
                HISFrameInfo info = readFrameHeader(raf);
                if (info != null) {
                    width = info.width;
                    height = info.height;
                    type = info.type;
                    nFrames = countFrames(raf, file.length());
                }
            } finally {
                raf.close();
            }
        } catch (IOException e) {
            /* If scan fails, just show dialog without frame info */
        }

        String typeStr = (type == TYPE_12BIT) ? "12-bit packed" : "16-bit";
        long memNeeded = (long) width * height * 2 * nFrames;
        String memStr;
        if (memNeeded > 1024L * 1024 * 1024) {
            memStr = String.format("%.1f GB", memNeeded / (1024.0 * 1024 * 1024));
        } else {
            memStr = String.format("%.1f MB", memNeeded / (1024.0 * 1024));
        }

        GenericDialog gd = new GenericDialog("Open HIS File");
        gd.addMessage("File: " + file.getName());
        if (nFrames > 0) {
            gd.addMessage(width + " x " + height + " x " + nFrames
                + " frames (" + typeStr + ")");
            gd.addMessage("File size: " + sizeStr
                + "  /  Memory if loaded: " + memStr);
        } else {
            gd.addMessage("File size: " + sizeStr);
        }
        gd.addCheckbox("Use Virtual Stack (memory-efficient)", false);
        gd.showDialog();

        if (gd.wasCanceled()) return false;

        useVirtualStack = gd.getNextBoolean();
        return true;
    }

    /**
     * Quickly count frames in a HIS file by scanning headers.
     * File position is reset to the beginning on return.
     */
    private int countFrames(RandomAccessFile raf, long fileLength) throws IOException {
        raf.seek(0);
        int count = 0;

        while (raf.getFilePointer() < fileLength) {
            HISFrameInfo info = readFrameHeader(raf);
            if (info == null) break;

            /* Skip past comment + pixel data */
            long dataSize;
            if (info.type == TYPE_16BIT) {
                dataSize = (long) info.width * info.height * 2;
            } else if (info.type == TYPE_12BIT) {
                dataSize = ((long) info.width * info.height * 3) / 2;
            } else {
                break;
            }

            long skipAmount = info.commentLength + dataSize;
            long newPos = raf.getFilePointer() + skipAmount;
            if (newPos > fileLength) break;
            raf.seek(newPos);
            count++;
        }

        raf.seek(0);
        return count;
    }

    /**
     * Try to open file by checking the "IM" magic bytes.
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
     * Frame header info
     * ================================================================ */

    /** Parsed frame header data. */
    private static class HISFrameInfo {
        int commentLength;
        int width;
        int height;
        int type;
        int nImage1;
        int nImage2;
    }

    /**
     * Read a 64-byte frame header at the current file position.
     * Returns null if not a valid "IM" header or EOF.
     * File position advances past the 64-byte header (but NOT past comment).
     */
    private static HISFrameInfo readFrameHeader(RandomAccessFile raf) throws IOException {
        byte[] buf = new byte[HEADER_SIZE];
        try {
            raf.readFully(buf);
        } catch (EOFException e) {
            return null;
        }

        ByteBuffer hdr = ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN);

        if (hdr.get(0) != 'I' || hdr.get(1) != 'M') return null;

        HISFrameInfo info = new HISFrameInfo();
        info.commentLength = hdr.getShort(2) & 0xFFFF;
        info.width         = hdr.getShort(4) & 0xFFFF;
        info.height        = hdr.getShort(6) & 0xFFFF;
        info.type          = hdr.getShort(12) & 0xFFFF;
        info.nImage1       = hdr.getShort(14) & 0xFFFF;
        info.nImage2       = hdr.getShort(16) & 0xFFFF;
        return info;
    }

    /* ================================================================
     * HIS format reader - normal (all frames in memory)
     * ================================================================ */

    private ImagePlus readHIS(File file) throws IOException {
        IJ.showStatus("Reading HIS file: " + file.getName());

        RandomAccessFile raf = new RandomAccessFile(file, "r");
        try {
            HISFrameInfo firstInfo = readFrameHeader(raf);
            if (firstInfo == null) {
                IJ.error("Open HIS", "Not a valid HIS file (missing 'IM' header).");
                return null;
            }

            int width  = firstInfo.width;
            int height = firstInfo.height;
            int type   = firstInfo.type;
            long nImages = firstInfo.nImage1 + 65536L * firstInfo.nImage2;

            if (width <= 0 || height <= 0) {
                IJ.error("Open HIS", "Invalid image dimensions: " + width + " x " + height);
                return null;
            }

            String typeStr = (type == TYPE_16BIT) ? "16-bit"
                           : (type == TYPE_12BIT) ? "12-bit packed" : "type=" + type;
            if (type != TYPE_16BIT && type != TYPE_12BIT) {
                IJ.error("Open HIS", "Unsupported data type: " + type);
                return null;
            }

            IJ.log("HIS file: " + file.getName());
            IJ.log("  Dimensions: " + width + " x " + height);
            IJ.log("  Data type: " + typeStr);
            IJ.log("  Declared frames: " + nImages);

            /* Skip comment of first header */
            if (firstInfo.commentLength > 0) {
                raf.skipBytes(firstInfo.commentLength);
            }

            ImageStack stack = new ImageStack(width, height);
            int pixelCount = width * height;
            int frameIndex = 0;

            /* Read first frame */
            short[] pixels = readFrameData(raf, type, pixelCount);
            if (pixels != null) {
                frameIndex++;
                stack.addSlice("Frame " + frameIndex,
                    new ShortProcessor(width, height, pixels, null));
                IJ.showProgress(frameIndex, (int) Math.max(nImages, 1));
            }

            /* Read subsequent frames */
            while (raf.getFilePointer() < raf.length()) {
                HISFrameInfo info = readFrameHeader(raf);
                if (info == null) {
                    IJ.log("  Warning: Frame header mismatch at frame "
                        + (frameIndex + 1) + ", stopping.");
                    break;
                }

                if (info.commentLength > 0) {
                    raf.skipBytes(info.commentLength);
                }

                pixels = readFrameData(raf, info.type, info.width * info.height);
                if (pixels == null) break;

                frameIndex++;
                stack.addSlice("Frame " + frameIndex,
                    new ShortProcessor(info.width, info.height, pixels, null));
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

    /* ================================================================
     * HIS format reader - Virtual Stack
     * ================================================================ */

    /**
     * Open HIS file as a Virtual Stack.
     * First pass: scan all frame headers to build an offset table.
     * Frames are then read on-demand from disk.
     */
    private ImagePlus readHISVirtual(File file) throws IOException {
        IJ.showStatus("Scanning HIS file: " + file.getName());

        /* First pass: build frame offset table */
        RandomAccessFile raf = new RandomAccessFile(file, "r");
        long fileLength = raf.length();

        ArrayList/*<Long>*/ offsets = new ArrayList/*<Long>*/();
        ArrayList/*<Integer>*/ types = new ArrayList/*<Integer>*/();
        int width = 0, height = 0;
        int firstType = 0;

        try {
            int idx = 0;
            while (raf.getFilePointer() < fileLength) {
                long frameStart = raf.getFilePointer();

                HISFrameInfo info = readFrameHeader(raf);
                if (info == null) break;

                if (idx == 0) {
                    width = info.width;
                    height = info.height;
                    firstType = info.type;
                }

                if (info.type != TYPE_16BIT && info.type != TYPE_12BIT) break;
                if (info.width <= 0 || info.height <= 0) break;

                /* Offset of pixel data = current pos + commentLength */
                long dataOffset = raf.getFilePointer() + info.commentLength;
                offsets.add(new Long(dataOffset));
                types.add(new Integer(info.type));

                /* Skip to next frame */
                long dataSize;
                if (info.type == TYPE_16BIT) {
                    dataSize = (long) info.width * info.height * 2;
                } else {
                    dataSize = ((long) info.width * info.height * 3) / 2;
                }
                long nextFrame = dataOffset + dataSize;
                if (nextFrame > fileLength) break;
                raf.seek(nextFrame);
                idx++;

                if (idx % 100 == 0) {
                    IJ.showStatus("Scanning HIS: " + idx + " frames...");
                }
            }
        } finally {
            raf.close();
        }

        int nFrames = offsets.size();
        if (nFrames == 0) {
            IJ.error("Open HIS", "No valid frames found in file.");
            return null;
        }

        String typeStr = (firstType == TYPE_16BIT) ? "16-bit" : "12-bit packed";

        IJ.log("HIS file (Virtual Stack): " + file.getName());
        IJ.log("  Dimensions: " + width + " x " + height);
        IJ.log("  Data type: " + typeStr);
        IJ.log("  Frames found: " + nFrames);

        /* Convert ArrayLists to arrays for the VirtualStack */
        long[] offsetArray = new long[nFrames];
        int[] typeArray = new int[nFrames];
        for (int i = 0; i < nFrames; i++) {
            offsetArray[i] = ((Long) offsets.get(i)).longValue();
            typeArray[i] = ((Integer) types.get(i)).intValue();
        }

        /* Create the virtual stack */
        HISVirtualStack vstack = new HISVirtualStack(
            file, width, height, offsetArray, typeArray);

        ImagePlus imp = new ImagePlus(file.getName() + " (Virtual)", vstack);
        imp.setProperty("Info",
            "Format: HIS (Hamamatsu Image Sequence) - Virtual Stack\n" +
            "Data type: " + typeStr + "\n" +
            "Frames: " + nFrames + "\n" +
            "Source: " + file.getAbsolutePath());
        IJ.showStatus("HIS Virtual Stack: " + nFrames + " frames");
        return imp;
    }

    /**
     * VirtualStack implementation for HIS files.
     * Each frame is read from disk on demand via RandomAccessFile.
     */
    static class HISVirtualStack extends VirtualStack {

        private final File file;
        private final int width;
        private final int height;
        private final long[] dataOffsets;   /* byte offset of pixel data for each frame */
        private final int[] dataTypes;      /* TYPE_16BIT or TYPE_12BIT per frame */
        private final int nFrames;

        HISVirtualStack(File file, int width, int height,
                         long[] dataOffsets, int[] dataTypes) {
            super(width, height, dataOffsets.length);
            this.file = file;
            this.width = width;
            this.height = height;
            this.dataOffsets = dataOffsets;
            this.dataTypes = dataTypes;
            this.nFrames = dataOffsets.length;
        }

        @Override
        public ImageProcessor getProcessor(int n) {
            /* n is 1-based */
            if (n < 1 || n > nFrames) {
                IJ.error("HIS Virtual Stack",
                    "Frame index out of range: " + n);
                return new ShortProcessor(width, height);
            }

            int idx = n - 1;
            int pixelCount = width * height;
            short[] pixels = null;

            try {
                RandomAccessFile raf = new RandomAccessFile(file, "r");
                try {
                    raf.seek(dataOffsets[idx]);
                    pixels = readFrameDataStatic(raf, dataTypes[idx], pixelCount);
                } finally {
                    raf.close();
                }
            } catch (IOException e) {
                IJ.log("HIS Virtual Stack: Error reading frame " + n
                    + ": " + e.getMessage());
            }

            if (pixels == null) {
                /* Return blank frame on error */
                return new ShortProcessor(width, height);
            }

            return new ShortProcessor(width, height, pixels, null);
        }

        @Override
        public int getSize() {
            return nFrames;
        }

        @Override
        public String getSliceLabel(int n) {
            return "Frame " + n;
        }
    }

    /* ================================================================
     * Frame data readers
     * ================================================================ */

    /**
     * Read one frame of pixel data (instance method for normal reading).
     */
    private short[] readFrameData(RandomAccessFile raf, int type, int pixelCount)
            throws IOException {
        return readFrameDataStatic(raf, type, pixelCount);
    }

    /**
     * Read one frame of pixel data from current file position.
     * Static method so it can be called from the VirtualStack inner class.
     * Returns 16-bit pixel array (short[]) or null on error/EOF.
     */
    static short[] readFrameDataStatic(RandomAccessFile raf, int type, int pixelCount)
            throws IOException {

        if (type == TYPE_16BIT) {
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

                int valA = b0 * 16 + b1 / 16;
                int valB = (b1 % 16) * 256 + b2;

                if (pixIdx < pixelCount) pixels[pixIdx++] = (short) valA;
                if (pixIdx < pixelCount) pixels[pixIdx++] = (short) valB;
            }
            return pixels;

        } else {
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