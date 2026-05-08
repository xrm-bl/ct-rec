#!/bin/bash
# ============================================================
#  build.sh - SP8CT ImageJ plugin build script (Linux/Mac)
#
#  Usage:  ./build.sh [path/to/ij.jar]
#
#  Run from the src directory where all .java files reside.
# ============================================================

set -e

JARNAME="SP8CT_Plugins.jar"
BUILDDIR="build"

# --- Locate ij.jar ---
if [ -n "$1" ]; then
    IJJAR="$1"
elif [ -f "../ij.jar" ]; then
    IJJAR="../ij.jar"
elif [ -f "../../ij.jar" ]; then
    IJJAR="../../ij.jar"
elif [ -n "$IMAGEJ_DIR" ] && [ -f "$IMAGEJ_DIR/ij.jar" ]; then
    IJJAR="$IMAGEJ_DIR/ij.jar"
else
    echo "ERROR: ij.jar not found."
    echo "Usage: ./build.sh path/to/ij.jar"
    exit 1
fi

echo "Using ij.jar: $IJJAR"

# --- Clean and create build directory ---
rm -rf "$BUILDDIR"
mkdir -p "$BUILDDIR"

# --- Compile all Java files ---
echo "Compiling..."
javac -cp "$IJJAR" -d "$BUILDDIR" -encoding UTF-8 -source 8 -target 8 \
    *.java

echo "Compilation OK."

# --- Copy plugins.config into build ---
cp plugins.config "$BUILDDIR/"

# --- Create JAR ---
# Include all .class EXCEPT HandleExtraFileTypes* (including $N inner classes)
echo "Creating $JARNAME..."
cd "$BUILDDIR"
ls *.class | grep -v '^HandleExtraFileTypes' > _jarfiles.txt
echo plugins.config >> _jarfiles.txt
jar cf "../$JARNAME" @_jarfiles.txt
rm _jarfiles.txt
cd ..

echo ""
echo "============================================================"
echo " Build complete: $JARNAME"
echo ""
echo " Install:"
echo "   1. Copy $JARNAME to ImageJ/plugins/"
echo "   2. Copy $BUILDDIR/HandleExtraFileTypes.class"
echo "      to ImageJ/plugins/ (for D&D support)"
echo "   3. Restart ImageJ"
echo "============================================================"