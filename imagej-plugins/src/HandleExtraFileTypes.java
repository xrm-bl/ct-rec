/*
 * HandleExtraFileTypes.java
 *
 * This plugin enables ImageJ/Fiji to open HIS and IMG files via
 * drag-and-drop, File > Open, or double-click.
 *
 * ImageJ's Opener calls this class (via IJ.runPlugIn) when it encounters
 * a file type it does not recognize. The key mechanism is:
 *
 *   1. HandleExtraFileTypes extends ImagePlus (not just implements PlugIn)
 *   2. run() reads the file and calls setStack() on itself (this)
 *   3. Opener receives this object as an ImagePlus with valid image data
 *
 * If this object has width > 0 after run(), Opener treats it as a
 * successfully loaded image. If width == 0, Opener shows the
 * "Format not supported" error.
 *
 * IMPORTANT:
 *   - For vanilla ImageJ: place this .class in plugins/ as a standalone file
 *     (NOT inside a jar). It must be in the default package (no package statement).
 *   - For Fiji: Fiji ships its own HandleExtraFileTypes in the IO_.jar plugin.
 *     You need to either:
 *       (a) Edit Fiji's HandleExtraFileTypes.java source to add HIS/IMG entries, or
 *       (b) Remove IO_.jar and use this file instead (not recommended), or
 *       (c) Use the menu-based approach (Plugins > SP8CT > IO > Open HIS/IMG...)
 *
 * Author: Generated for SPring-8 BL47XU
 * No external dependencies other than ij.jar
 */

import ij.*;
import ij.plugin.*;
import java.io.*;

public class HandleExtraFileTypes extends ImagePlus implements PlugIn {

    static final int IMAGE_OPENED = -1;
    static final int PLUGIN_NOT_FOUND = -2;

    /**
     * Called from io/Opener.java when ImageJ cannot identify a file type.
     *
     * @param path  full path to the file being opened
     */
    @Override
    public void run(String path) {
        if (path == null || path.equals("")) return;

        File theFile = new File(path);
        String fileName = theFile.getName();
        String lowerName = fileName.toLowerCase();

        /* Only handle .his, .img, .kif, and .fpi extensions */
        if (!lowerName.endsWith(".his") && !lowerName.endsWith(".img")
                && !lowerName.endsWith(".kif") && !lowerName.endsWith(".fpi")) {
            return;  /* width remains 0 -> Opener knows we didn't handle it */
        }

        /* Delegate to Open_HIS_IMG to do the actual reading */
        Open_HIS_IMG opener = new Open_HIS_IMG();
        ImagePlus imp = opener.openFile(path);

        if (imp == null) {
            /* Reading failed or user cancelled.
             * width remains 0 -> Opener may show "Format not supported".
             * This is acceptable for user-initiated cancel.
             */
            return;
        }

        /* Transfer the image data into this HandleExtraFileTypes object.
         * This is the standard mechanism: Opener casts the return of
         * IJ.runPlugIn() to ImagePlus and checks getWidth() > 0.
         */
        ImageStack stack = imp.getStack();
        String title = imp.getTitle();
        if (title == null || title.isEmpty()) title = fileName;

        setStack(title, stack);

        /* Copy calibration info */
        setCalibration(imp.getCalibration());

        /* Copy "Show Info" property */
        if (imp.getProperty("Info") != null) {
            setProperty("Info", imp.getProperty("Info"));
        }

        /* Copy FileInfo if available */
        if (imp.getOriginalFileInfo() != null) {
            setFileInfo(imp.getOriginalFileInfo());
        }
    }
}
