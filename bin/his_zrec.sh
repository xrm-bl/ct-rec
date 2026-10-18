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
# actually used, read back from the ImageDescription with pid; trailing
# zeros removed, e.g. rc-check/003_01024_1023.5.tif), so that trial
# reconstructions of several layers and axis positions can be compared
# side by side.  The program can be changed with HIS_REC (default
# his_rec_g_c, or his_rec_t_c when no GPU build is on the PATH); a
# 360-degree offset scan is handled by his_rec itself.
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
SSS=$(printf "%03d" "$SCAN")
LLLLL=$(printf "%05d" "$LAYER")
REC="rec$LLLLL.tif"
mkdir -p rc-check

echo "==== scan $SSS, layer $LLLLL ===="
if ! RHIS_SCAN=$SCAN "$HIS_REC" "$LAYER" ${CENTER:+"$CENTER"}; then
    echo "error: $HIS_REC failed" >&2
    exit 1
fi
if [ ! -f "$REC" ]; then
    echo "error: $REC was not written" >&2
    exit 1
fi

# axis actually used: column 2 of the ImageDescription (pid), trailing zeros trimmed
C=$(pid "$REC" 2>/dev/null | awk '{print $2}')
[ -z "$C" ] && C="${CENTER:-auto}"
C=$(echo "$C" | sed -e 's/\(\.[0-9]*[1-9]\)0*$/\1/' -e 's/\.0*$//')

OUT="rc-check/${SSS}_${LLLLL}_${C}.tif"
mv -f "$REC" "$OUT"
echo "$OUT"
