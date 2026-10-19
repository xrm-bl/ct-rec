#!/bin/bash
# rc-check-his.sh - rotation-axis check of a Z scan directly from a.his
#
#   rc-check-his.sh [layer]
#     layer : layer reconstructed for the check (default: the centre row
#             of the frame, Ny/2, as his_info reports it)
#
# Run in the directory holding a.his, conv.bat and output.log.  his_info
# reads the his header, conv.bat and output.log and reports the number of
# scans M in the file (frames / frames per scan), the centre layer and the
# scan kind (180 / 360).  For every scan k the environment variable
# RHIS_SCAN=k selects the frames of that scan (no spl / his_spl_E, no
# conv.bat run, no .img files), his_rec reconstructs the layer with its
# own axis estimate, the ImageDescription (pid) goes to center.log and the
# image to rc-check/00k.tif for a visual check.
#
# center.log holds, in column 2, the value to pass to his_srec: for a
# 180-degree scan that is the ImageDescription value itself (ct_rec and
# hp_tg share the convention); for a 360-degree offset scan the
# ImageDescription of his_rec (ofct_rec) carries the image centre
# dcnt+0.5 while his_srec (ofct_srec) wants the detector axis position
# Rc = dcnt, so 0.5 is subtracted here.  gen-all-his.sh reads center.log.
# The program can be changed with HIS_REC (default his_rec_g_c, or
# his_rec_t_c when no GPU build is on the PATH).
set -u
if [ $# -gt 1 ] || [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    echo "usage: $(basename "$0") [layer]"
    echo "  layer : layer used for the check (default: centre row of the frame)"
    echo "  run where a.his, conv.bat and output.log are; the number of scans and the"
    echo "  scan kind come from his_info; writes center.log (column 2 = value for"
    echo "  his_srec) and rc-check/00k.tif"
    exit 1
fi
if [ -z "${HIS_REC:-}" ]; then
    if command -v his_rec_g_c > /dev/null 2>&1; then HIS_REC=his_rec_g_c; else HIS_REC=his_rec_t_c; fi
fi
if [ ! -f conv.bat ]; then
    echo "error: conv.bat not found here (run beside a.his, conv.bat and output.log)" >&2
    exit 1
fi
INFO=$(his_info) || { echo "error: his_info failed" >&2; exit 1; }
read -r NX NY NF L N LAYER0 KIND <<< "$INFO"
LAYER="${1:-$LAYER0}"
echo "$NX x $NY, $NF frames: $N scan(s) of $L frames; $KIND-degree scan; check layer $LAYER"
REC=$(printf "rec%05d.tif" "$LAYER")
rm -f center.log
mkdir -p rc-check

for ((i = 1; i <= N; i++)); do
    DIR=$(printf "%03d" "$i")
    echo "==== scan $DIR / $N ===="
    if ! RHIS_SCAN=$i "$HIS_REC" "$LAYER"; then
        echo "error: $HIS_REC failed on scan $i" >&2
        exit 1
    fi
    # ImageDescription: Dr centre views RA0 min max; column 2 -> the his_srec value
    pid "$REC" | awk -v kind="$KIND" '{ v = $2; if (kind == 360) v = v - 0.5; printf "%s %.3f %s %s %s %s\n", $1, v, $3, $4, $5, $6 }' >> center.log
    echo "axis for his_srec: $(tail -1 center.log | awk '{print $2}')"
    mv -f "$REC" "rc-check/$DIR.tif"
done
echo "done: center.log ($N lines) and rc-check/*.tif written; next: gen-all-his.sh center.log Dr"
