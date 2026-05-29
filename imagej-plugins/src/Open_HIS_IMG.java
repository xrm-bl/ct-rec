/*
 * Open_HIS_IMG.java
 *
 * ImageJ / Fiji plugin for opening HIS, IMG, and KIF format files
 * used at SPring-8 synchrotron radiation facility.
 *
 * HIS format: Hamamatsu Image Sequence format
 *   - 64-byte header per frame (magic "IM")
 *   - type=2: 16-bit unsigned, type=6: 12-bit packed
 *   - Multiple frames, Virtual Stack support
 *
 * IMG format: HiPic single-image format
 *   - 64-byte header (magic "IM"), 16-bit unsigned, single frame
 *
 * KIF format: Hamamatsu KAISHIN Image File format (v1.0/v2.0)
 *   - 96-byte file header (magic ".KIF") + 32-byte frame headers ("FH")
 *   - 8/16/32-bit unsigned/signed/float, grayscale/RGB
 *   - DateTimeOffset timestamps, Virtual Stack support
 *
 * Based on his2tif6.c, img_ave.c, his_time.c (xrm-bl/ct-rec)
 * and Hamamatsu KIF File Format Specification v2.0 (2021-08-27)
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
import java.util.Calendar;
import java.util.TimeZone;

public class Open_HIS_IMG implements PlugIn {

    /* ================================================================
     * Constants - HIS/IMG
     * ================================================================ */
    private static final int HIS_HEADER_SIZE = 64;
    private static final int TYPE_16BIT = 2;
    private static final int TYPE_12BIT = 6;

    /* ================================================================
     * Constants - KIF
     * ================================================================ */
    private static final int KIF_FILE_HEADER_SIZE = 96;
    private static final int KIF_FRAME_HEADER_SIZE = 32;
    /* .NET epoch offset: ticks from 0001-01-01 to 1970-01-01 */
    private static final long DOTNET_EPOCH_OFFSET = 621355968000000000L;

    /* KIF sample format IDs */
    private static final int KIF_SAMPLE_UNSIGNED = 1;
    private static final int KIF_SAMPLE_SIGNED   = 2;
    private static final int KIF_SAMPLE_FLOAT    = 3;

    /* KIF color space IDs */
    private static final int KIF_COLOR_GRAY = 1;
    private static final int KIF_COLOR_RGB  = 2;

    /** Virtual stack flag set by options dialog */
    private boolean useVirtualStack = false;

    /* ================================================================
     * PlugIn entry point
     * ================================================================ */

    @Override
    public void run(String arg) {
        String path = arg;

        if (path == null || path.isEmpty()) {
            OpenDialog od = new OpenDialog("Open HIS/IMG/KIF File", null);
            String dir = od.getDirectory();
            String name = od.getFileName();
            if (name == null) return;
            path = dir + name;
        }

        File file = new File(path);
        if (!file.exists()) {
            IJ.error("Open HIS/IMG/KIF", "File not found: " + path);
            return;
        }

        String lowerName = file.getName().toLowerCase();
        if (lowerName.endsWith(".his")) {
            if (!showHISOptionsDialog(file)) return;
        } else if (lowerName.endsWith(".kif")) {
            if (!showKIFOptionsDialog(file)) return;
        }

        ImagePlus imp = openFileInternal(file);
        if (imp != null) {
            imp.show();
        }
    }

    /**
     * Open a file and return ImagePlus without showing it.
     * Called by HandleExtraFileTypes for D&D / File > Open.
     */
    public ImagePlus openFile(String path) {
        useVirtualStack = false;
        File file = new File(path);
        if (!file.exists()) {
            IJ.error("Open HIS/IMG/KIF", "File not found: " + path);
            return null;
        }

        String lowerName = file.getName().toLowerCase();
        if (lowerName.endsWith(".his")) {
            if (!showHISOptionsDialog(file)) return null;
        } else if (lowerName.endsWith(".kif")) {
            if (!showKIFOptionsDialog(file)) return null;
        }

        return openFileInternal(file);
    }

    private ImagePlus openFileInternal(File file) {
        String lowerName = file.getName().toLowerCase();

        try {
            if (lowerName.endsWith(".his")) {
                return useVirtualStack ? readHISVirtual(file) : readHIS(file);
            } else if (lowerName.endsWith(".img")) {
                return readIMG(file);
            } else if (lowerName.endsWith(".kif")) {
                return useVirtualStack ? readKIFVirtual(file) : readKIF(file);
            } else {
                return tryOpenByMagic(file);
            }
        } catch (IOException e) {
            IJ.error("Open HIS/IMG/KIF", "Error reading file:\n" + e.getMessage());
            return null;
        }
    }

    private ImagePlus tryOpenByMagic(File file) throws IOException {
        RandomAccessFile raf = new RandomAccessFile(file, "r");
        try {
            if (raf.length() < 4) return null;
            byte[] magic = new byte[4];
            raf.readFully(magic);

            /* Check for KIF magic ".KIF" */
            if (magic[0] == '.' && magic[1] == 'K' && magic[2] == 'I' && magic[3] == 'F') {
                raf.close();
                if (!showKIFOptionsDialog(file)) return null;
                return useVirtualStack ? readKIFVirtual(file) : readKIF(file);
            }

            /* Check for HIS/IMG magic "IM" */
            if (magic[0] == 'I' && magic[1] == 'M') {
                raf.close();
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

            return null;
        } finally {
            try { raf.close(); } catch (IOException e) { /* already closed */ }
        }
    }


    /* ================================================================
     *  ================================================================
     *   HIS FORMAT
     *  ================================================================
     * ================================================================ */

    /* --- HIS Options Dialog --- */

    private boolean showHISOptionsDialog(File file) {
        long fileSize = file.length();
        String sizeStr = formatFileSize(fileSize);

        int nFrames = 0;
        int width = 0, height = 0, type = 0;
        double firstTS = 0, lastTS = 0;
        try {
            RandomAccessFile raf = new RandomAccessFile(file, "r");
            try {
                HISFrameInfo info = readHISFrameHeader(raf);
                if (info != null) {
                    width = info.width;
                    height = info.height;
                    type = info.type;
                    HISScanResult scan = scanHISFrames(raf, file.length());
                    nFrames = scan.nFrames;
                    firstTS = scan.firstTimeStamp;
                    lastTS  = scan.lastTimeStamp;
                }
            } finally {
                raf.close();
            }
        } catch (IOException e) { /* show dialog without frame info */ }

        String typeStr = (type == TYPE_12BIT) ? "12-bit packed" : "16-bit";
        long memNeeded = (long) width * height * 2 * nFrames;
        String memStr = formatFileSize(memNeeded);

        GenericDialog gd = new GenericDialog("Open HIS File");
        gd.addMessage("File: " + file.getName());
        if (nFrames > 0) {
            gd.addMessage(width + " x " + height + " x " + nFrames
                + " frames (" + typeStr + ")");
            gd.addMessage("File size: " + sizeStr + "  /  Memory if loaded: " + memStr);
            if (nFrames > 1) {
                double span = lastTS - firstTS;
                double interval = span / (nFrames - 1);
                gd.addMessage("Time span: " + formatTime(span)
                    + "  /  Avg interval: " + formatTime(interval));
            }
        } else {
            gd.addMessage("File size: " + sizeStr);
        }
        gd.addCheckbox("Use Virtual Stack (memory-efficient)", false);
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        useVirtualStack = gd.getNextBoolean();
        return true;
    }

    /* --- HIS Frame Header --- */

    private static class HISFrameInfo {
        int commentLength;
        int width;
        int height;
        int type;
        int nImage1;
        int nImage2;
        int[] timeBytes = new int[8];
        double timeStamp;
    }

    private static HISFrameInfo readHISFrameHeader(RandomAccessFile raf) throws IOException {
        byte[] buf = new byte[HIS_HEADER_SIZE];
        try { raf.readFully(buf); } catch (EOFException e) { return null; }

        ByteBuffer hdr = ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN);
        if (hdr.get(0) != 'I' || hdr.get(1) != 'M') return null;

        HISFrameInfo info = new HISFrameInfo();
        info.commentLength = hdr.getShort(2) & 0xFFFF;
        info.width         = hdr.getShort(4) & 0xFFFF;
        info.height        = hdr.getShort(6) & 0xFFFF;
        info.type          = hdr.getShort(12) & 0xFFFF;
        info.nImage1       = hdr.getShort(14) & 0xFFFF;
        info.nImage2       = hdr.getShort(16) & 0xFFFF;

        /* Read timestamp as 8 individual bytes (per his_time.c) */
        byte[] tsBytes = new byte[8];
        for (int i = 0; i < 8; i++) {
            tsBytes[i] = buf[22 + i];
            info.timeBytes[i] = buf[22 + i] & 0xFF;
        }
        long bits = 0;
        for (int i = 7; i >= 0; i--) {
            bits = (bits << 8) | (tsBytes[i] & 0xFFL);
        }
        info.timeStamp = Double.longBitsToDouble(bits);
        /* Raw value is in milliseconds; convert to seconds */
        info.timeStamp /= 1000.0;

        return info;
    }

    /* --- HIS Scanner --- */

    private static class HISScanResult {
        int nFrames;
        double firstTimeStamp;
        double lastTimeStamp;
    }

    private HISScanResult scanHISFrames(RandomAccessFile raf, long fileLength) throws IOException {
        raf.seek(0);
        HISScanResult result = new HISScanResult();
        int count = 0;

        while (raf.getFilePointer() < fileLength) {
            HISFrameInfo info = readHISFrameHeader(raf);
            if (info == null) break;

            if (count == 0) result.firstTimeStamp = info.timeStamp;
            result.lastTimeStamp = info.timeStamp;

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

        result.nFrames = count;
        raf.seek(0);
        return result;
    }

    /* --- HIS Normal Reader --- */

    private ImagePlus readHIS(File file) throws IOException {
        IJ.showStatus("Reading HIS file: " + file.getName());

        RandomAccessFile raf = new RandomAccessFile(file, "r");
        try {
            HISFrameInfo firstInfo = readHISFrameHeader(raf);
            if (firstInfo == null) {
                IJ.error("Open HIS", "Not a valid HIS file.");
                return null;
            }

            int width  = firstInfo.width;
            int height = firstInfo.height;
            int type   = firstInfo.type;
            long nImages = firstInfo.nImage1 + 65536L * firstInfo.nImage2;

            if (width <= 0 || height <= 0) {
                IJ.error("Open HIS", "Invalid dimensions: " + width + " x " + height);
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
            IJ.log("  First frame timestamp: " + firstInfo.timeStamp
                + " (raw bytes: " + formatTimeBytes(firstInfo.timeBytes) + ")");

            double t0 = firstInfo.timeStamp;

            if (firstInfo.commentLength > 0) {
                raf.skipBytes(firstInfo.commentLength);
            }

            ImageStack stack = new ImageStack(width, height);
            int pixelCount = width * height;
            int frameIndex = 0;

            short[] pixels = readHISFrameData(raf, type, pixelCount);
            if (pixels != null) {
                frameIndex++;
                double dt = firstInfo.timeStamp - t0;
                stack.addSlice(formatSliceLabel(frameIndex, dt),
                    new ShortProcessor(width, height, pixels, null));
                IJ.showProgress(frameIndex, (int) Math.max(nImages, 1));
            }

            while (raf.getFilePointer() < raf.length()) {
                HISFrameInfo info = readHISFrameHeader(raf);
                if (info == null) {
                    IJ.log("  Warning: Frame header mismatch at frame "
                        + (frameIndex + 1) + ", stopping.");
                    break;
                }
                if (info.commentLength > 0) raf.skipBytes(info.commentLength);

                pixels = readHISFrameData(raf, info.type, info.width * info.height);
                if (pixels == null) break;

                frameIndex++;
                double dt = info.timeStamp - t0;
                stack.addSlice(formatSliceLabel(frameIndex, dt),
                    new ShortProcessor(info.width, info.height, pixels, null));
                IJ.showProgress(frameIndex, (int) Math.max(nImages, 1));
            }

            IJ.log("  Frames read: " + frameIndex);
            if (stack.getSize() == 0) {
                IJ.error("Open HIS", "No valid frames found.");
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

    /* --- HIS Virtual Stack Reader --- */

    private ImagePlus readHISVirtual(File file) throws IOException {
        IJ.showStatus("Scanning HIS file: " + file.getName());

        RandomAccessFile raf = new RandomAccessFile(file, "r");
        long fileLength = raf.length();

        ArrayList offsets = new ArrayList();
        ArrayList types = new ArrayList();
        ArrayList timestamps = new ArrayList();
        int width = 0, height = 0;
        int firstType = 0;
        int[] firstTimeBytes = null;

        try {
            int idx = 0;
            while (raf.getFilePointer() < fileLength) {
                HISFrameInfo info = readHISFrameHeader(raf);
                if (info == null) break;

                if (idx == 0) {
                    width = info.width;
                    height = info.height;
                    firstType = info.type;
                    firstTimeBytes = info.timeBytes;
                }

                if (info.type != TYPE_16BIT && info.type != TYPE_12BIT) break;
                if (info.width <= 0 || info.height <= 0) break;

                long dataOffset = raf.getFilePointer() + info.commentLength;
                offsets.add(new Long(dataOffset));
                types.add(new Integer(info.type));
                timestamps.add(new Double(info.timeStamp));

                long dataSize = (info.type == TYPE_16BIT)
                    ? (long) info.width * info.height * 2
                    : ((long) info.width * info.height * 3) / 2;
                long nextFrame = dataOffset + dataSize;
                if (nextFrame > fileLength) break;
                raf.seek(nextFrame);
                idx++;
                if (idx % 100 == 0) IJ.showStatus("Scanning HIS: " + idx + " frames...");
            }
        } finally {
            raf.close();
        }

        int nFrames = offsets.size();
        if (nFrames == 0) {
            IJ.error("Open HIS", "No valid frames found.");
            return null;
        }

        String typeStr = (firstType == TYPE_16BIT) ? "16-bit" : "12-bit packed";

        IJ.log("HIS file (Virtual Stack): " + file.getName());
        IJ.log("  Dimensions: " + width + " x " + height);
        IJ.log("  Data type: " + typeStr);
        IJ.log("  Frames found: " + nFrames);

        long[] offsetArray = new long[nFrames];
        int[] typeArray = new int[nFrames];
        double[] tsArray = new double[nFrames];
        for (int i = 0; i < nFrames; i++) {
            offsetArray[i] = ((Long) offsets.get(i)).longValue();
            typeArray[i] = ((Integer) types.get(i)).intValue();
            tsArray[i] = ((Double) timestamps.get(i)).doubleValue();
        }

        IJ.log("  First frame timestamp: " + tsArray[0]
            + (firstTimeBytes != null ?
                " (raw bytes: " + formatTimeBytes(firstTimeBytes) + ")" : ""));
        if (nFrames > 1) {
            IJ.log("  Last frame timestamp:  " + tsArray[nFrames - 1]);
            IJ.log("  Total time span: " + formatTime(tsArray[nFrames - 1] - tsArray[0]));
        }

        HISVirtualStack vstack = new HISVirtualStack(
            file, width, height, offsetArray, typeArray, tsArray);

        ImagePlus imp = new ImagePlus(file.getName() + " (Virtual)", vstack);
        imp.setProperty("Info",
            "Format: HIS (Hamamatsu Image Sequence) - Virtual Stack\n" +
            "Data type: " + typeStr + "\nFrames: " + nFrames + "\n" +
            "Source: " + file.getAbsolutePath());
        IJ.showStatus("HIS Virtual Stack: " + nFrames + " frames");
        return imp;
    }

    static class HISVirtualStack extends VirtualStack {
        private final File file;
        private final int width, height, nFrames;
        private final long[] dataOffsets;
        private final int[] dataTypes;
        private final double[] timeStamps;
        private final double t0;

        HISVirtualStack(File file, int width, int height,
                         long[] dataOffsets, int[] dataTypes, double[] timeStamps) {
            super(width, height, dataOffsets.length);
            this.file = file;
            this.width = width;  this.height = height;
            this.dataOffsets = dataOffsets;  this.dataTypes = dataTypes;
            this.timeStamps = timeStamps;
            this.t0 = (timeStamps.length > 0) ? timeStamps[0] : 0.0;
            this.nFrames = dataOffsets.length;
        }

        public ImageProcessor getProcessor(int n) {
            if (n < 1 || n > nFrames) {
                IJ.error("HIS Virtual Stack", "Frame index out of range: " + n);
                return new ShortProcessor(width, height);
            }
            int idx = n - 1;
            short[] pixels = null;
            try {
                RandomAccessFile raf = new RandomAccessFile(file, "r");
                try {
                    raf.seek(dataOffsets[idx]);
                    pixels = readHISFrameDataStatic(raf, dataTypes[idx], width * height);
                } finally { raf.close(); }
            } catch (IOException e) {
                IJ.log("HIS Virtual Stack: Error reading frame " + n + ": " + e.getMessage());
            }
            if (pixels == null) return new ShortProcessor(width, height);
            return new ShortProcessor(width, height, pixels, null);
        }

        public int getSize() { return nFrames; }

        public String getSliceLabel(int n) {
            if (n >= 1 && n <= nFrames) {
                return formatSliceLabel(n, timeStamps[n - 1] - t0);
            }
            return "Frame " + n;
        }
    }

    /* --- HIS Frame Data Readers --- */

    private short[] readHISFrameData(RandomAccessFile raf, int type, int pixelCount)
            throws IOException {
        return readHISFrameDataStatic(raf, type, pixelCount);
    }

    static short[] readHISFrameDataStatic(RandomAccessFile raf, int type, int pixelCount)
            throws IOException {

        if (type == TYPE_16BIT) {
            int byteCount = pixelCount * 2;
            byte[] raw = new byte[byteCount];
            try { raf.readFully(raw); } catch (EOFException e) { return null; }

            short[] pixels = new short[pixelCount];
            ByteBuffer bb = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);
            for (int i = 0; i < pixelCount; i++) pixels[i] = bb.getShort();
            return pixels;

        } else if (type == TYPE_12BIT) {
            int packedBytes = (pixelCount * 3) / 2;
            byte[] raw = new byte[packedBytes];
            try { raf.readFully(raw); } catch (EOFException e) { return null; }

            short[] pixels = new short[pixelCount];
            int pixIdx = 0;
            for (int i = 0; i < packedBytes; i += 3) {
                if (i + 2 >= packedBytes) break;
                int b0 = raw[i] & 0xFF, b1 = raw[i+1] & 0xFF, b2 = raw[i+2] & 0xFF;
                if (pixIdx < pixelCount) pixels[pixIdx++] = (short) (b0 * 16 + b1 / 16);
                if (pixIdx < pixelCount) pixels[pixIdx++] = (short) ((b1 % 16) * 256 + b2);
            }
            return pixels;
        }
        return null;
    }

    /* ================================================================
     * IMG format reader (single frame, HiPic format)
     * ================================================================ */

    private ImagePlus readIMG(File file) throws IOException {
        IJ.showStatus("Reading IMG file: " + file.getName());

        RandomAccessFile raf = new RandomAccessFile(file, "r");
        try {
            byte[] headerBuf = new byte[HIS_HEADER_SIZE];
            raf.readFully(headerBuf);
            ByteBuffer hdr = ByteBuffer.wrap(headerBuf).order(ByteOrder.LITTLE_ENDIAN);

            if (hdr.get(0) != 'I' || hdr.get(1) != 'M') {
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
                IJ.error("Open IMG", "Invalid dimensions: " + width + " x " + height);
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
            if (!comment.isEmpty()) IJ.log("  Comment: " + comment);

            int pixelCount = width * height;
            short[] pixels = readHISFrameData(raf, TYPE_16BIT, pixelCount);
            if (pixels == null) {
                IJ.error("Open IMG", "Error reading pixel data.");
                return null;
            }

            ShortProcessor sp = new ShortProcessor(width, height, pixels, null);
            ImagePlus imp = new ImagePlus(file.getName(), sp);
            imp.setProperty("Info",
                "Format: IMG (HiPic)\nType: " + type +
                "\nOffset: (" + xOffset + ", " + yOffset + ")" +
                "\nComment: " + comment + "\nSource: " + file.getAbsolutePath());
            IJ.showStatus("IMG file loaded: " + width + " x " + height);
            return imp;
        } finally {
            raf.close();
        }
    }


    /* ================================================================
     *  ================================================================
     *   KIF FORMAT
     *  ================================================================
     * ================================================================ */

    /* --- KIF File Header --- */

    private static class KIFFileInfo {
        int versionMajor, versionMinor;
        long nFrames;       /* 0 = unknown */
        int width, height;
        int bitsPerSample;
        int samplesPerPixel;
        int sampleFormat;   /* 1=uint, 2=int, 3=float */
        int colorSpace;     /* 1=gray, 2=RGB */
        int planarConfig;   /* 1=chunky, 2=planar */
        int compression;    /* 0=none */
        int fileCommentLength;
        long indexTableOffset;
        long customDataOffset;
    }

    private static KIFFileInfo readKIFFileHeader(RandomAccessFile raf) throws IOException {
        byte[] buf = new byte[KIF_FILE_HEADER_SIZE];
        raf.readFully(buf);

        /* Check magic ".KIF" */
        if (buf[0] != '.' || buf[1] != 'K' || buf[2] != 'I' || buf[3] != 'F') {
            return null;
        }

        ByteBuffer hdr = ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN);

        KIFFileInfo info = new KIFFileInfo();
        info.versionMajor     = hdr.get(4) & 0xFF;
        info.versionMinor     = hdr.get(5) & 0xFF;
        info.nFrames          = hdr.getInt(6) & 0xFFFFFFFFL;
        info.width            = hdr.getInt(10);
        info.height           = hdr.getInt(14);
        info.bitsPerSample    = hdr.getShort(18) & 0xFFFF;
        info.samplesPerPixel  = hdr.getShort(20) & 0xFFFF;
        info.sampleFormat     = hdr.getShort(22) & 0xFFFF;
        info.colorSpace       = hdr.getShort(24) & 0xFFFF;
        info.planarConfig     = hdr.getShort(26) & 0xFFFF;
        info.compression      = hdr.getShort(28) & 0xFFFF;
        info.fileCommentLength = hdr.getInt(46);
        info.indexTableOffset = hdr.getLong(50);
        info.customDataOffset = hdr.getLong(58);

        return info;
    }

    /* --- KIF Frame Header --- */

    private static class KIFFrameInfo {
        long ticks;             /* .NET DateTimeOffset ticks (100ns since 0001-01-01) */
        short utcOffsetMinutes; /* UTC offset in minutes */
        double timeSeconds;     /* Converted to seconds since epoch for relative calc */
        int commentLength;
        int compressedDataLength;
    }

    private static KIFFrameInfo readKIFFrameHeader(RandomAccessFile raf) throws IOException {
        byte[] buf = new byte[KIF_FRAME_HEADER_SIZE];
        try { raf.readFully(buf); } catch (EOFException e) { return null; }

        /* Check magic "FH" */
        if (buf[0] != 'F' || buf[1] != 'H') return null;

        ByteBuffer hdr = ByteBuffer.wrap(buf).order(ByteOrder.LITTLE_ENDIAN);

        KIFFrameInfo info = new KIFFrameInfo();
        /* DateTimeOffset: bytes 4-11 = Int64 ticks, bytes 12-13 = Int16 UTC offset */
        info.ticks = hdr.getLong(4);
        info.utcOffsetMinutes = hdr.getShort(12);
        info.commentLength = hdr.getInt(16);
        info.compressedDataLength = hdr.getInt(20);

        /* Convert ticks to seconds (since .NET epoch) for relative time calculation */
        info.timeSeconds = info.ticks / 10000000.0;

        return info;
    }

    /** Format a KIF DateTimeOffset as a human-readable string */
    private static String formatKIFTimestamp(long ticks, short utcOffsetMinutes) {
        /* Convert .NET ticks to Java millis */
        long javaMillis = (ticks - DOTNET_EPOCH_OFFSET) / 10000L;

        Calendar cal = Calendar.getInstance(TimeZone.getTimeZone("UTC"));
        cal.setTimeInMillis(javaMillis + utcOffsetMinutes * 60000L);

        int yr = cal.get(Calendar.YEAR);
        int mo = cal.get(Calendar.MONTH) + 1;
        int dy = cal.get(Calendar.DAY_OF_MONTH);
        int hr = cal.get(Calendar.HOUR_OF_DAY);
        int mi = cal.get(Calendar.MINUTE);
        int sc = cal.get(Calendar.SECOND);
        int ms = cal.get(Calendar.MILLISECOND);

        String sign = (utcOffsetMinutes >= 0) ? "+" : "-";
        int absOff = Math.abs(utcOffsetMinutes);
        int offH = absOff / 60;
        int offM = absOff % 60;

        return String.format("%04d-%02d-%02d %02d:%02d:%02d.%03d %s%02d:%02d",
            yr, mo, dy, hr, mi, sc, ms, sign, offH, offM);
    }

    /* --- KIF Pixel Format Description --- */

    private static String describeKIFPixelFormat(KIFFileInfo fi) {
        String fmt;
        if (fi.sampleFormat == KIF_SAMPLE_UNSIGNED) fmt = "uint";
        else if (fi.sampleFormat == KIF_SAMPLE_SIGNED) fmt = "int";
        else if (fi.sampleFormat == KIF_SAMPLE_FLOAT) fmt = "float";
        else fmt = "unknown";

        String color;
        if (fi.colorSpace == KIF_COLOR_RGB) color = "RGB";
        else color = "gray";

        return fi.bitsPerSample + "-bit " + fmt + " " + color
            + " (" + fi.samplesPerPixel + " spp)";
    }

    /** Calculate pixel data size in bytes for one frame */
    private static long kifFrameDataSize(KIFFileInfo fi) {
        return (long) fi.width * fi.height * fi.samplesPerPixel
            * (fi.bitsPerSample / 8);
    }

    /* --- KIF Options Dialog --- */

    private boolean showKIFOptionsDialog(File file) {
        long fileSize = file.length();
        String sizeStr = formatFileSize(fileSize);

        KIFFileInfo fi = null;
        int nFrames = 0;
        String firstTS = "", lastTS = "";
        double firstTimeSec = 0, lastTimeSec = 0;

        try {
            RandomAccessFile raf = new RandomAccessFile(file, "r");
            try {
                fi = readKIFFileHeader(raf);
                if (fi != null) {
                    /* Skip file comment */
                    if (fi.fileCommentLength > 0) raf.skipBytes(fi.fileCommentLength);

                    /* Scan frames to count and get timestamps */
                    long dataSize = kifFrameDataSize(fi);
                    long fileLen = raf.length();
                    int count = 0;

                    while (raf.getFilePointer() < fileLen) {
                        KIFFrameInfo frame = readKIFFrameHeader(raf);
                        if (frame == null) break;

                        if (count == 0) {
                            firstTS = formatKIFTimestamp(frame.ticks, frame.utcOffsetMinutes);
                            firstTimeSec = frame.timeSeconds;
                        }
                        lastTS = formatKIFTimestamp(frame.ticks, frame.utcOffsetMinutes);
                        lastTimeSec = frame.timeSeconds;

                        /* Skip comment + pixel data */
                        long skip = frame.commentLength + dataSize;
                        long newPos = raf.getFilePointer() + skip;
                        if (newPos > fileLen) break;
                        raf.seek(newPos);
                        count++;
                    }
                    nFrames = count;
                }
            } finally {
                raf.close();
            }
        } catch (IOException e) { /* show dialog without scan info */ }

        GenericDialog gd = new GenericDialog("Open KIF File");
        gd.addMessage("File: " + file.getName());
        if (fi != null && nFrames > 0) {
            gd.addMessage(fi.width + " x " + fi.height + " x " + nFrames
                + " frames");
            gd.addMessage("Pixel: " + describeKIFPixelFormat(fi)
                + "  /  KIF v" + fi.versionMajor + "." + fi.versionMinor);
            long memNeeded = kifFrameDataSize(fi) * nFrames;
            gd.addMessage("File size: " + sizeStr
                + "  /  Memory if loaded: " + formatFileSize(memNeeded));
            if (nFrames > 1) {
                double span = lastTimeSec - firstTimeSec;
                double interval = span / (nFrames - 1);
                gd.addMessage("Time span: " + formatTime(span)
                    + "  /  Avg interval: " + formatTime(interval));
            }
            if (!firstTS.isEmpty()) {
                gd.addMessage("First frame: " + firstTS);
            }
        } else if (fi != null) {
            gd.addMessage(fi.width + " x " + fi.height
                + "  " + describeKIFPixelFormat(fi));
            gd.addMessage("File size: " + sizeStr);
        } else {
            gd.addMessage("File size: " + sizeStr);
        }
        gd.addCheckbox("Use Virtual Stack (memory-efficient)", false);
        gd.showDialog();
        if (gd.wasCanceled()) return false;
        useVirtualStack = gd.getNextBoolean();
        return true;
    }

    /* --- KIF Normal Reader --- */

    private ImagePlus readKIF(File file) throws IOException {
        IJ.showStatus("Reading KIF file: " + file.getName());

        RandomAccessFile raf = new RandomAccessFile(file, "r");
        try {
            KIFFileInfo fi = readKIFFileHeader(raf);
            if (fi == null) {
                IJ.error("Open KIF", "Not a valid KIF file.");
                return null;
            }

            if (fi.width <= 0 || fi.height <= 0) {
                IJ.error("Open KIF", "Invalid dimensions: " + fi.width + " x " + fi.height);
                return null;
            }
            if (fi.compression != 0) {
                IJ.error("Open KIF", "Compressed KIF data is not supported.");
                return null;
            }

            IJ.log("KIF file: " + file.getName());
            IJ.log("  Version: " + fi.versionMajor + "." + fi.versionMinor);
            IJ.log("  Dimensions: " + fi.width + " x " + fi.height);
            IJ.log("  Pixel format: " + describeKIFPixelFormat(fi));
            IJ.log("  Declared frames: " + fi.nFrames);

            /* Skip file comment */
            if (fi.fileCommentLength > 0) raf.skipBytes(fi.fileCommentLength);

            long dataSize = kifFrameDataSize(fi);
            ImageStack stack = null;
            int frameIndex = 0;
            double t0 = 0;

            while (raf.getFilePointer() < raf.length()) {
                KIFFrameInfo frame = readKIFFrameHeader(raf);
                if (frame == null) break;

                if (frameIndex == 0) {
                    t0 = frame.timeSeconds;
                    IJ.log("  First frame: "
                        + formatKIFTimestamp(frame.ticks, frame.utcOffsetMinutes));
                }

                /* Skip frame comment */
                if (frame.commentLength > 0) raf.skipBytes(frame.commentLength);

                /* Read pixel data */
                ImageProcessor ip = readKIFPixelData(raf, fi);
                if (ip == null) break;

                if (stack == null) stack = new ImageStack(fi.width, fi.height);

                frameIndex++;
                double dt = frame.timeSeconds - t0;
                stack.addSlice(formatSliceLabel(frameIndex, dt), ip);
                IJ.showProgress(frameIndex, (int) Math.max(fi.nFrames, 1));
            }

            IJ.log("  Frames read: " + frameIndex);
            if (stack == null || stack.getSize() == 0) {
                IJ.error("Open KIF", "No valid frames found.");
                return null;
            }

            ImagePlus imp = new ImagePlus(file.getName(), stack);
            imp.setProperty("Info",
                "Format: KIF (Hamamatsu KAISHIN) v" + fi.versionMajor + "." + fi.versionMinor + "\n" +
                "Pixel: " + describeKIFPixelFormat(fi) + "\n" +
                "Frames: " + frameIndex + "\n" +
                "Source: " + file.getAbsolutePath());
            IJ.showStatus("KIF file loaded: " + frameIndex + " frames");
            return imp;
        } finally {
            raf.close();
        }
    }

    /* --- KIF Virtual Stack Reader --- */

    private ImagePlus readKIFVirtual(File file) throws IOException {
        IJ.showStatus("Scanning KIF file: " + file.getName());

        RandomAccessFile raf = new RandomAccessFile(file, "r");
        KIFFileInfo fi;
        try {
            fi = readKIFFileHeader(raf);
            if (fi == null) {
                IJ.error("Open KIF", "Not a valid KIF file.");
                return null;
            }
            if (fi.compression != 0) {
                IJ.error("Open KIF", "Compressed KIF data is not supported.");
                return null;
            }
        } catch (IOException e) {
            raf.close();
            throw e;
        }

        long dataSize = kifFrameDataSize(fi);
        long fileLength = raf.length();

        ArrayList offsets = new ArrayList();   /* Long: pixel data offsets */
        ArrayList timestamps = new ArrayList(); /* Double: timeSeconds */
        String firstTS = "";

        try {
            /* Skip file comment */
            if (fi.fileCommentLength > 0) raf.skipBytes(fi.fileCommentLength);

            int idx = 0;
            while (raf.getFilePointer() < fileLength) {
                KIFFrameInfo frame = readKIFFrameHeader(raf);
                if (frame == null) break;

                if (idx == 0) {
                    firstTS = formatKIFTimestamp(frame.ticks, frame.utcOffsetMinutes);
                }

                long pixelOffset = raf.getFilePointer() + frame.commentLength;
                offsets.add(new Long(pixelOffset));
                timestamps.add(new Double(frame.timeSeconds));

                long nextFrame = pixelOffset + dataSize;
                if (nextFrame > fileLength) break;
                raf.seek(nextFrame);
                idx++;
                if (idx % 100 == 0) IJ.showStatus("Scanning KIF: " + idx + " frames...");
            }
        } finally {
            raf.close();
        }

        int nFrames = offsets.size();
        if (nFrames == 0) {
            IJ.error("Open KIF", "No valid frames found.");
            return null;
        }

        IJ.log("KIF file (Virtual Stack): " + file.getName());
        IJ.log("  Version: " + fi.versionMajor + "." + fi.versionMinor);
        IJ.log("  Dimensions: " + fi.width + " x " + fi.height);
        IJ.log("  Pixel format: " + describeKIFPixelFormat(fi));
        IJ.log("  Frames found: " + nFrames);
        if (!firstTS.isEmpty()) IJ.log("  First frame: " + firstTS);

        long[] offsetArray = new long[nFrames];
        double[] tsArray = new double[nFrames];
        for (int i = 0; i < nFrames; i++) {
            offsetArray[i] = ((Long) offsets.get(i)).longValue();
            tsArray[i] = ((Double) timestamps.get(i)).doubleValue();
        }

        if (nFrames > 1) {
            IJ.log("  Total time span: " + formatTime(tsArray[nFrames-1] - tsArray[0]));
        }

        KIFVirtualStack vstack = new KIFVirtualStack(file, fi, offsetArray, tsArray);

        ImagePlus imp = new ImagePlus(file.getName() + " (Virtual)", vstack);
        imp.setProperty("Info",
            "Format: KIF (Hamamatsu KAISHIN) v" + fi.versionMajor + "." + fi.versionMinor
            + " - Virtual Stack\n" +
            "Pixel: " + describeKIFPixelFormat(fi) + "\n" +
            "Frames: " + nFrames + "\n" +
            "Source: " + file.getAbsolutePath());
        IJ.showStatus("KIF Virtual Stack: " + nFrames + " frames");
        return imp;
    }

    static class KIFVirtualStack extends VirtualStack {
        private final File file;
        private final int width, height, bitsPerSample, samplesPerPixel, sampleFormat, colorSpace;
        private final int nFrames;
        private final long[] dataOffsets;
        private final double[] timeStamps;
        private final double t0;

        KIFVirtualStack(File file, KIFFileInfo fi, long[] dataOffsets, double[] timeStamps) {
            super(fi.width, fi.height, dataOffsets.length);
            this.file = file;
            this.width = fi.width;  this.height = fi.height;
            this.bitsPerSample = fi.bitsPerSample;
            this.samplesPerPixel = fi.samplesPerPixel;
            this.sampleFormat = fi.sampleFormat;
            this.colorSpace = fi.colorSpace;
            this.dataOffsets = dataOffsets;
            this.timeStamps = timeStamps;
            this.t0 = (timeStamps.length > 0) ? timeStamps[0] : 0.0;
            this.nFrames = dataOffsets.length;
        }

        public ImageProcessor getProcessor(int n) {
            if (n < 1 || n > nFrames) {
                IJ.error("KIF Virtual Stack", "Frame index out of range: " + n);
                return createBlankProcessor();
            }
            try {
                RandomAccessFile raf = new RandomAccessFile(file, "r");
                try {
                    raf.seek(dataOffsets[n - 1]);
                    ImageProcessor ip = readKIFPixelDataStatic(raf,
                        width, height, bitsPerSample, samplesPerPixel,
                        sampleFormat, colorSpace);
                    if (ip != null) return ip;
                } finally { raf.close(); }
            } catch (IOException e) {
                IJ.log("KIF Virtual Stack: Error reading frame " + n + ": " + e.getMessage());
            }
            return createBlankProcessor();
        }

        private ImageProcessor createBlankProcessor() {
            if (colorSpace == KIF_COLOR_RGB && samplesPerPixel >= 3)
                return new ColorProcessor(width, height);
            if (sampleFormat == KIF_SAMPLE_FLOAT || bitsPerSample == 32)
                return new FloatProcessor(width, height);
            if (bitsPerSample <= 8)
                return new ByteProcessor(width, height);
            return new ShortProcessor(width, height);
        }

        public int getSize() { return nFrames; }

        public String getSliceLabel(int n) {
            if (n >= 1 && n <= nFrames) {
                return formatSliceLabel(n, timeStamps[n - 1] - t0);
            }
            return "Frame " + n;
        }
    }

    /* --- KIF Pixel Data Readers --- */

    private ImageProcessor readKIFPixelData(RandomAccessFile raf, KIFFileInfo fi)
            throws IOException {
        return readKIFPixelDataStatic(raf,
            fi.width, fi.height, fi.bitsPerSample, fi.samplesPerPixel,
            fi.sampleFormat, fi.colorSpace);
    }

    static ImageProcessor readKIFPixelDataStatic(
            RandomAccessFile raf,
            int width, int height, int bitsPerSample, int samplesPerPixel,
            int sampleFormat, int colorSpace) throws IOException {

        int pixelCount = width * height;
        int bytesPerSample = bitsPerSample / 8;
        int bytesPerFrame = pixelCount * samplesPerPixel * bytesPerSample;

        byte[] raw = new byte[bytesPerFrame];
        try { raf.readFully(raw); } catch (EOFException e) { return null; }

        ByteBuffer bb = ByteBuffer.wrap(raw).order(ByteOrder.LITTLE_ENDIAN);

        /* --- RGB 8-bit chunky --- */
        if (colorSpace == KIF_COLOR_RGB && samplesPerPixel >= 3
                && bitsPerSample == 8 && sampleFormat == KIF_SAMPLE_UNSIGNED) {
            int[] rgbPixels = new int[pixelCount];
            for (int i = 0; i < pixelCount; i++) {
                int r = raw[i * samplesPerPixel]     & 0xFF;
                int g = raw[i * samplesPerPixel + 1] & 0xFF;
                int b = raw[i * samplesPerPixel + 2] & 0xFF;
                rgbPixels[i] = (r << 16) | (g << 8) | b;
            }
            return new ColorProcessor(width, height, rgbPixels);
        }

        /* --- Grayscale 8-bit unsigned --- */
        if (bitsPerSample == 8 && sampleFormat == KIF_SAMPLE_UNSIGNED
                && samplesPerPixel == 1) {
            byte[] pixels = new byte[pixelCount];
            System.arraycopy(raw, 0, pixels, 0, pixelCount);
            return new ByteProcessor(width, height, pixels, null);
        }

        /* --- Grayscale 16-bit unsigned --- */
        if (bitsPerSample == 16 && sampleFormat == KIF_SAMPLE_UNSIGNED
                && samplesPerPixel == 1) {
            short[] pixels = new short[pixelCount];
            for (int i = 0; i < pixelCount; i++) pixels[i] = bb.getShort();
            return new ShortProcessor(width, height, pixels, null);
        }

        /* --- Grayscale 16-bit signed -> store as ShortProcessor --- */
        if (bitsPerSample == 16 && sampleFormat == KIF_SAMPLE_SIGNED
                && samplesPerPixel == 1) {
            short[] pixels = new short[pixelCount];
            for (int i = 0; i < pixelCount; i++) pixels[i] = bb.getShort();
            /* ImageJ ShortProcessor stores unsigned 0-65535;
             * shift signed range -32768..32767 to 0..65535 */
            for (int i = 0; i < pixelCount; i++) {
                pixels[i] = (short) (pixels[i] + 32768);
            }
            ShortProcessor sp = new ShortProcessor(width, height, pixels, null);
            sp.setMinAndMax(0, 65535);
            return sp;
        }

        /* --- Grayscale 32-bit float --- */
        if (bitsPerSample == 32 && sampleFormat == KIF_SAMPLE_FLOAT
                && samplesPerPixel == 1) {
            float[] pixels = new float[pixelCount];
            for (int i = 0; i < pixelCount; i++) pixels[i] = bb.getFloat();
            return new FloatProcessor(width, height, pixels, null);
        }

        /* --- Grayscale 32-bit unsigned/signed -> convert to float --- */
        if (bitsPerSample == 32 && samplesPerPixel == 1) {
            float[] pixels = new float[pixelCount];
            if (sampleFormat == KIF_SAMPLE_UNSIGNED) {
                for (int i = 0; i < pixelCount; i++) {
                    pixels[i] = (float) (bb.getInt() & 0xFFFFFFFFL);
                }
            } else {
                for (int i = 0; i < pixelCount; i++) {
                    pixels[i] = (float) bb.getInt();
                }
            }
            return new FloatProcessor(width, height, pixels, null);
        }

        /* --- Fallback: unsupported format --- */
        IJ.log("  Unsupported KIF pixel format: " + bitsPerSample + "bps, "
            + samplesPerPixel + "spp, fmt=" + sampleFormat + ", cs=" + colorSpace);
        return null;
    }


    /* ================================================================
     * Common formatting helpers
     * ================================================================ */

    static String formatSliceLabel(int frameNumber, double dt) {
        return "Frame " + frameNumber + " / t=" + formatTime(dt);
    }

    static String formatTime(double seconds) {
        double abs = Math.abs(seconds);
        if (abs < 0.001 && abs > 0.0) {
            return String.format("%.1fus", seconds * 1e6);
        } else if (abs < 1.0) {
            return String.format("%.3fms", seconds * 1e3);
        } else if (abs < 60.0) {
            return String.format("%.3fs", seconds);
        } else if (abs < 3600.0) {
            int min = (int) (seconds / 60);
            double sec = seconds - min * 60;
            return String.format("%dm%.1fs", min, sec);
        } else {
            int hr = (int) (seconds / 3600);
            int min = (int) ((seconds - hr * 3600) / 60);
            double sec = seconds - hr * 3600 - min * 60;
            return String.format("%dh%dm%.1fs", hr, min, sec);
        }
    }

    static String formatTimeBytes(int[] b) {
        if (b == null || b.length < 8) return "N/A";
        return String.format("%02X %02X %02X %02X %02X %02X %02X %02X",
            b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]);
    }

    static String formatFileSize(long bytes) {
        if (bytes > 1024L * 1024 * 1024) {
            return String.format("%.1f GB", bytes / (1024.0 * 1024 * 1024));
        } else {
            return String.format("%.1f MB", bytes / (1024.0 * 1024));
        }
    }
}
