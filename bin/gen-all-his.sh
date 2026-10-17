#!/bin/bash
# gen-all-his.sh - whole-volume reconstruction of a Z scan directly from a.his
#
#   gen-all-his.sh center.log Dr [LACmin LACmax]
#     center.log : written by rc-check-his.sh (one line per scan, column 2 = axis)
#     Dr         : pixel size
#     LACmin/max : range of the 8-bit conversion (default -0.5 3.0)
#
# Run in the directory holding a.his, conv.bat and output.log (the same
# place as rc-check-his.sh).  Line k of center.log belongs to scan k: with
# RHIS_SCAN=k his_srec reconstructs that scan into 00k/rec and tif_f2i
# converts it to 8-bit in 00k/ro_xy, as gen-all.sh does after spl +
# conv.bat + hp_tg (the ro_zx reslice of gen-all.sh is not made).  The
# directories 00k/rec are what rec_stk expects.  The program can be
# changed with HIS_SREC (default his_srec_g_c, or his_srec_t_c when no GPU
# build is on the PATH); a 360-degree offset scan is handled by his_srec.
set -u
if [ -z "${2:-}" ]; then
    echo "usage: $(basename "$0") center.log Dr [LACmin LACmax]"
    echo "  center.log : from rc-check-his.sh, Dr : pixel size, LACmin LACmax : 8-bit range (default -0.5 3.0)"
    echo "  run where a.his, conv.bat and output.log are; writes 00k/rec and 00k/ro_xy"
    exit 1
fi
LOGFILE="$1"
DR="$2"
LMIN="${3:--0.5}"
LMAX="${4:-3.0}"
if [ -z "${HIS_SREC:-}" ]; then
    if command -v his_srec_g_c > /dev/null 2>&1; then HIS_SREC=his_srec_g_c; else HIS_SREC=his_srec_t_c; fi
fi
if [ ! -f "$LOGFILE" ]; then
    echo "error: $LOGFILE not found" >&2
    exit 1
fi
LINE_NUM=0

# center.log line by line (the || [ -n ] clause keeps a last line without newline)
while IFS=$' \t' read -r COL1 COL2 _REST || [ -n "$COL1" ]; do
    [ -z "$COL1" ] && continue
    LINE_NUM=$((LINE_NUM + 1))
    DIR=$(printf "%03d" "$LINE_NUM")
    echo "==== scan $DIR (center $COL2) ===="
    mkdir -p "$DIR/rec" "$DIR/ro_xy"
    if ! RHIS_SCAN=$LINE_NUM "$HIS_SREC" . "$DR" "$COL2" 0 "$DIR/rec"; then
        echo "error: $HIS_SREC failed on scan $LINE_NUM" >&2
        exit 1
    fi
    ( cd "$DIR" && tif_f2i 8 rec ro_xy "$LMIN" "$LMAX" )
done < "$LOGFILE"
echo "done: $LINE_NUM scans; next: rec_stk $LINE_NUM start end > ddd.sh"
