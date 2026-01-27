#!/bin/bash

# ==== Argument check ====
if [ -z "$2" ]; then
    echo "Usage: $0 <number_of_frames> <number_of_measurements>"
    echo "Example: $0 2101 4"
    exit 1
fi

N="$2"
BASEDIR="$(pwd)"

# Generate split data (same behavior as: spl %1 %2 > aaa.bat && call aaa.bat)
spl "$1" "$2" > aaa.sh
sed -i 's/copy/cp/g' aaa.sh
chmod +x aaa.sh
./aaa.sh

# Create output directory
mkdir -p rc-check

# ==== Loop from 001 to N ====
for ((i=1; i<=N; i++)); do
    DIR=$(printf "%03d" "$i")
    echo "==== Processing $DIR ===="

    if [ ! -d "$DIR/raw" ]; then
        echo "Error: $DIR/raw not found. Skipping."
    else
        cd "$DIR/raw" || continue

        # Execute conv.bat equivalent (assuming conv.bat is executable or replaced)
        sh ./conv.bat

        # Trial reconstruction
        ct_rec_g_c 120

        # Append rotation center info
        pid rec00120.tif >> "$BASEDIR/center.log"

        # Copy preview image
        cp rec00120.tif "$BASEDIR/rc-check/${DIR}.tif"

        cd "$BASEDIR"
    fi
done

echo "Completed."
