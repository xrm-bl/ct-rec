#!/bin/bash
# his_zrec.sh - one slice of one scan of a Z scan, kept in rc-check/
#
#   his_zrec.sh scan layer [center]
#     scan   : scan number in a.his (1..M, see his_info)
#     layer  : layer to reconstruct
#     center : rotation axis position (optional; his_rec estimates it
#              from that layer when omitted)
#
# Run in the directory holding a.his, conv.bat and output.log.  With
# RHIS_SCAN=scan his_rec reconstructs the layer and the image is moved to
# rc-check/SSS_LLLLL_C.tif (scan as %03d, layer as %05d, C = the axis
# actually used, as the value to pass to his_srec: the ImageDescription
# read back with pid, minus 0.5 for a 360-degree offset scan where
# ofct_rec records the image centre dcnt+0.5; trailing zeros removed,
# e.g. rc-check/003_01024_1023.5.tif), so that trial reconstructions of
# several layers and axis positions can be compared side by side.  The
# program can be changed with HIS_REC (default his_rec_g_c, or his_rec_t_c
# when no GPU build is on the PATH).
set -u
if [ $# -lt 2 ] || [ $# -gt 3 ]; then
    echo "usage: $(basename "$0") scan layer [center]"
    echo "  scan   : scan number in a.his (1..M, see his_info)"
    echo "  layer  : layer to reconstruct"
    echo "  center : rotation axis position (optional; estimated when omitted)"
    echo "  run where a.his, conv.bat and output.log are; writes rc-check/SSS_LLLLL_C.tif"
    exit 1
fi
SCAN="$1"; LAYER="$2"; CENTER="${3:-}"
if [ -z "${HIS_REC:-}" ]; then
    if command -v his_rec_g_c > /dev/null 2>&1; then HIS_REC=his_rec_g_c; else HIS_REC=his_rec_t_c; fi
fi
if [ ! -f conv.bat ]; then
    echo "error: conv.bat not found here (run beside a.his, conv.bat and output.log)" >&2
    exit 1
fi
KIND=$(his_info 2>/dev/null | awk '{print $7}')
[ -z "$KIND" ] && { echo "error: his_info failed" >&2; exit 1; }
SSS=$(printf "%03d" "$SCAN")
LLLLL=$(printf "%05d" "$LAYER")
REC="rec$LLLLL.tif"
mkdir -p rc-check

echo "==== scan $SSS, layer $LLLLL ($KIND-degree scan) ===="
if ! RHIS_SCAN=$SCAN "$HIS_REC" "$LAYER" ${CENTER:+"$CENTER"}; then
    echo "error: $HIS_REC failed" >&2
    exit 1
fi
if [ ! -f "$REC" ]; then
    echo "error: $REC was not written" >&2
    exit 1
fi

# axis actually used: column 2 of the ImageDescription (pid), -0.5 for an offset scan, trailing zeros trimmed
C=$(pid "$REC" 2>/dev/null | awk -v kind="$KIND" '{ v = $2; if (kind == 360) v = v - 0.5; printf "%.3f", v }')
[ -z "$C" ] && C="${CENTER:-auto}"
C=$(echo "$C" | sed -e 's/\(\.[0-9]*[1-9]\)0*$/\1/' -e 's/\.0*$//')

OUT="rc-check/${SSS}_${LLLLL}_${C}.tif"
mv -f "$REC" "$OUT"
echo "$OUT"
