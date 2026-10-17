#!/bin/bash
# rc-check-his.sh - rotation-axis check of a Z scan directly from a.his
#
#   rc-check-his.sh [layer]
#     layer : layer reconstructed for the check (default: the centre row
#             of the frame, Ny/2, as his_info reports it)
#
# Run in the directory holding a.his, conv.bat and output.log.  his_info
# reads the his header and conv.bat and reports the number of scans M in
# the file (frames / frames per scan) and the centre layer.  For every
# scan k the environment variable RHIS_SCAN=k selects the frames of that
# scan (no spl / his_spl_E, no conv.bat run, no .img files), his_rec
# reconstructs the layer with its own axis estimate, pid appends the
# ImageDescription (column 2 = axis) to center.log and the image goes to
# rc-check/00k.tif for a visual check.  center.log is the input of
# gen-all-his.sh.  The program can be changed with HIS_REC (default
# his_rec_g_c, or his_rec_t_c when no GPU build is on the PATH); his_rec
# picks ct_rec or ofct_rec from the angles in output.log, so 180-degree
# and 360-degree offset Z scans both work.
set -u
if [ $# -gt 1 ] || [ "${1:-}" = "-h" ] || [ "${1:-}" = "--help" ]; then
    echo "usage: $(basename "$0") [layer]"
    echo "  layer : layer used for the check (default: centre row of the frame)"
    echo "  run where a.his, conv.bat and output.log are; the number of scans comes"
    echo "  from his_info; writes center.log and rc-check/00k.tif"
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
read -r NX NY NF L N LAYER0 <<< "$INFO"
LAYER="${1:-$LAYER0}"
echo "$NX x $NY, $NF frames: $N scan(s) of $L frames; check layer $LAYER"
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
    pid "$REC" >> center.log
    mv -f "$REC" "rc-check/$DIR.tif"
done
echo "done: center.log ($N lines) and rc-check/*.tif written; next: gen-all-his.sh center.log Dr"
