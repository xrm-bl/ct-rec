#!/bin/bash

# ==== Argument check ====
if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_center.log>"
    echo "Example: $0 center.log"
    exit 1
fi

LOGFILE="$1"
BASEDIR="$(pwd)"
LINE_NUM=0

# ==== Read center.log line by line ====
# This avoids issues with missing newline at EOF
while IFS=$'\t '" " read -r COL1 COL2 COL3 COL4 COL5 COL6 || [ -n "$COL1" ]; do
    # Skip empty lines
    if [ -z "$COL1" ]; then
        continue
    fi

    LINE_NUM=$((LINE_NUM + 1))
    DIR=$(printf "%03d" "$LINE_NUM")

    echo "==== Processing ${DIR} (center=${COL2}) ===="

    if [ ! -d "$DIR" ]; then
        echo "Error: ${DIR} not found. Skipping."
        continue
    fi

    pushd "$DIR" >/dev/null

    # Create directories (ignore errors if they already exist)
    mkdir -p rec ro_xy ro_zx

    # Reconstruction commands
    hp_tg_g_c raw 5.64 "$COL2" 0 rec
    tif_f2i 8 rec ro_xy -0.5 3.0
    si_rar ro_xy - +x +z +y ro_zx

    popd >/dev/null

done < "$LOGFILE"

echo "Completed."
