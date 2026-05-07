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
echo "Creating $JARNAME..."
cd "$BUILDDIR"
jar cf "../$JARNAME" .
cd ..

echo ""
echo "============================================================"
echo " Build complete: $JARNAME"
echo ""
echo " Install:"
echo "   Copy $JARNAME to ImageJ/plugins/ or Fiji.app/plugins/"
echo "   and restart ImageJ/Fiji."
echo "============================================================"
