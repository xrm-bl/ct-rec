#!/bin/bash
# ----------------------------------------------------------------------
# ct_rec_loop.sh  (bash port of ct_rec_loop.bat)
#
# Reproduce hp_tg's whole-volume reconstruction by calling the per-slice
# program ct_rec once per layer, running several jobs in parallel.
# Arguments follow hp_tg; the LAST argument is the number of parallel
# jobs (optional, default 1).
#
#   ct_rec_loop.sh  indir  Dr  RC  RA0  outdir  [Njobs]
#   ct_rec_loop.sh  indir  Dr  L1 C1 L2 C2  RA0  outdir  [Njobs]
#
#   indir  : directory holding dark.img or dark.tif, q*.img / q*.tif,
#            and output.log (ct_rec auto-detects img/tif from the dark file)
#   Dr     : pixel size [um]
#   RC     : rotation center (fixed for every layer)
#   L1 C1  : center C1 at layer L1   (linear-center form)
#   L2 C2  : center C2 at layer L2 -> center(z)=C1+(z-L1)(C2-C1)/(L2-L1)
#   RA0    : offset angle [deg]
#   outdir : directory to collect rec*.tif into
#   Njobs  : number of concurrent jobs (default 1; tune to GPU memory)
#
# The executable can be overridden by the CT_EXE environment variable:
#   CT_EXE=~/bin/ct_rec_g_r ct_rec_loop.sh raw 1.0 4630.5 0.0 rec 20
#
# Run this from the directory ONE LEVEL ABOVE indir.
# ----------------------------------------------------------------------
set -u

# ct_rec executable (on PATH, or override with CT_EXE=...)
CT_EXE=${CT_EXE:-ct_rec_g_c}

usage() {
    echo "usage: $(basename "$0") indir Dr RC RA0 outdir [Njobs]"        >&2
    echo "       $(basename "$0") indir Dr L1 C1 L2 C2 RA0 outdir [Njobs]" >&2
    exit 1
}

# ---- parse arguments (hp_tg style + optional trailing Njobs) ----
if   [ $# -ge 8 ]; then
    MODE=LINEAR
    INDIR=$1; DR=$2; Z1=$3; C1=$4; Z2=$5; C2=$6; RA0=$7; OUTDIR=$8; NJ=${9:-1}
elif [ $# -ge 5 ]; then
    MODE=FIXED
    INDIR=$1; DR=$2; C1=$3; RA0=$4; OUTDIR=$5; NJ=${6:-1}
else
    usage
fi
case $NJ in (*[!0-9]*|'') NJ=1;; esac
[ "$NJ" -ge 1 ] || NJ=1

# absolute paths (workers cd into INDIR, so OUTDIR/log must be absolute)
INDIR_ABS=$(cd "$INDIR" 2>/dev/null && pwd) ||
    { echo "no such directory: $INDIR" >&2; exit 1; }
mkdir -p "$OUTDIR" || exit 1
OUTDIR_ABS=$(cd "$OUTDIR" && pwd)
LOGDIR_ABS=$OUTDIR_ABS/log
mkdir -p "$LOGDIR_ABS"

# ---- dark file: img preferred, tif fallback (as in ct_rec itself) ----
if   [ -f "$INDIR_ABS/dark.img" ]; then DARKFILE=dark.img
elif [ -f "$INDIR_ABS/dark.tif" ]; then DARKFILE=dark.tif
else echo "no dark.img / dark.tif in $INDIR_ABS" >&2; exit 1; fi

# ---- binary readers (od is POSIX; le=1 for little-endian) ----
rd_u16() { # file le offset
    local f=$1 le=$2 o=$3
    if [ "$le" = 1 ]; then od -An -tu2 -j"$o" -N2 "$f" | tr -d ' '
    else set -- $(od -An -tu1 -j"$o" -N2 "$f"); echo $(( $1*256 + $2 )); fi
}
rd_u32() { # file le offset
    local f=$1 le=$2 o=$3
    if [ "$le" = 1 ]; then od -An -tu4 -j"$o" -N4 "$f" | tr -d ' '
    else set -- $(od -An -tu1 -j"$o" -N4 "$f")
         echo $(( (($1*256 + $2)*256 + $3)*256 + $4 )); fi
}

# height = ImageLength (tag 257) from the first IFD of a TIFF (II and MM)
tif_height() {
    local f=$1 le=0 ifd n i e tag typ
    [ "$(od -An -tx1 -N2 "$f" | tr -d ' \n')" = "4949" ] && le=1
    ifd=$(rd_u32 "$f" $le 4); n=$(rd_u16 "$f" $le "$ifd")
    for ((i=0; i<n; i++)); do
        e=$((ifd + 2 + 12*i))
        tag=$(rd_u16 "$f" $le $e)
        if [ "$tag" -eq 257 ]; then
            typ=$(rd_u16 "$f" $le $((e+2)))
            if [ "$typ" -eq 3 ]; then rd_u16 "$f" $le $((e+8))
            else                      rd_u32 "$f" $le $((e+8)); fi
            return 0
        fi
    done
    return 1
}

# ---- fixed-center form reconstructs every layer (0..height-1) ----
if [ "$MODE" = FIXED ]; then
    if [ "$DARKFILE" = dark.tif ]; then
        HEIGHT=$(tif_height "$INDIR_ABS/dark.tif")
    else
        # 16-bit value at byte offset 6 of the HiPic .img header
        HEIGHT=$(rd_u16 "$INDIR_ABS/dark.img" 1 6)
    fi
    if [ -z "${HEIGHT:-}" ] || [ "$HEIGHT" -le 0 ] 2>/dev/null; then
        echo "cannot read height from $DARKFILE" >&2; exit 1
    fi
    Z1=0; Z2=$((HEIGHT-1)); C2=$C1
fi

echo "ct_rec_loop: $MODE  layers $Z1..$Z2  Njobs=$NJ  indir=\"$INDIR_ABS\"  out=\"$OUTDIR_ABS\""

# ---- worker: run ct_rec for one layer in INDIR, then move the result ----
run_one() { # z zp center
    ( cd "$INDIR_ABS" &&
      "$CT_EXE" "$1" "$3" "$DR" "$RA0" > "$LOGDIR_ABS/rec$2.log" 2>&1
      mv -f "rec$2.tif" "$OUTDIR_ABS/" 2>/dev/null )
}

while read -r z zp c; do
    while (( $(jobs -rp | wc -l) >= NJ )); do wait -n; done
    run_one "$z" "$zp" "$c" &
done < <(awk -v z1="$Z1" -v z2="$Z2" -v c1="$C1" -v c2="$C2" -v m="$MODE" 'BEGIN{
    for (z=z1; z<=z2; z++) {
        c = (m=="LINEAR" && z2!=z1) ? c1+(z-z1)*(c2-c1)/(z2-z1) : c1
        printf "%d %05d %.4f\n", z, z, c
    }}')

wait
echo "done -> $OUTDIR_ABS/"
