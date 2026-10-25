#!/usr/bin/env bash
# make_n2i_portable.sh -- build a portable n2i bundle on an internet-connected Linux PC.
#
# The bundle contains a relocatable Python (python-build-standalone), PyTorch (CUDA build),
# numpy, tifffile and n2i_v3.py. Copy n2i_portable.tar.gz to the target machine, extract it
# anywhere and run ./n2i_portable/n2i -- nothing is installed there.
#
# Requirements: build PC = Linux x86_64 with internet, curl, python3 (any version, only to
# read the GitHub release list), tar. Target = Linux x86_64, glibc >= 2.28, NVIDIA driver.
#
# Usage:
#   bash make_n2i_portable.sh                     # CUDA 12.8 build (driver "CUDA Version" >= 12.8)
#   bash make_n2i_portable.sh --cuda-tag cu126    # other CUDA build
#   bash make_n2i_portable.sh --python 3.11       # other Python minor version
#   bash make_n2i_portable.sh --out mydir         # output directory name

set -u
CUDA_TAG="cu128"
PYMINOR="3.12"
OUT="n2i_portable"
ARCH="x86_64-unknown-linux-gnu"

info() { printf '\033[36m[n2i]\033[0m %s\n' "$*"; }
fail() { printf '\033[31m[n2i] ERROR:\033[0m %s\n' "$*"; exit 1; }

while [ $# -gt 0 ]; do
    case "$1" in
        --cuda-tag) CUDA_TAG="$2"; shift 2 ;;
        --python)   PYMINOR="$2"; shift 2 ;;
        --out)      OUT="$2"; shift 2 ;;
        -h|--help)  sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)          fail "unknown option: $1" ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
[ -f "$SCRIPT_DIR/n2i_v3.py" ] || fail "n2i_v3.py not found next to this script."
[ "$(uname -m)" = "x86_64" ] || fail "this script builds an x86_64 bundle; run it on an x86_64 PC."
command -v curl >/dev/null || fail "curl not found (sudo apt install curl)."
command -v python3 >/dev/null || fail "python3 not found (needed only to read the release list)."
[ -e "$OUT" ] && fail "$OUT already exists; remove it or use --out."

# ------------------------------------------------------------------ 1. standalone Python
info "Looking up the latest python-build-standalone release (Python $PYMINOR) ..."
URL=$(curl -fsSL https://api.github.com/repos/astral-sh/python-build-standalone/releases/latest | \
      python3 -c "
import json, re, sys
pat = re.compile(r'cpython-' + re.escape(sys.argv[1]) + r'\.\d+\+\d+-' + re.escape(sys.argv[2]) + r'-install_only\.tar\.gz$')
for a in json.load(sys.stdin)['assets']:
    if pat.search(a['name']):
        print(a['browser_download_url']); break
" "$PYMINOR" "$ARCH")
[ -n "$URL" ] || fail "no Python $PYMINOR build found in the latest release (network / proxy?)."
info "Downloading $(basename "$URL") ..."
mkdir -p "$OUT"
curl -fL "$URL" | tar -xz -C "$OUT" || fail "download/extract of Python failed."
PY="$OUT/python/bin/python3"
[ -x "$PY" ] || fail "unexpected archive layout (python/bin/python3 missing)."
"$PY" -m pip install --upgrade pip --quiet || fail "pip upgrade failed."

# ------------------------------------------------------------------ 2. packages
if [ "${N2I_SKIP_TORCH:-0}" != "1" ]; then
    info "Installing PyTorch ($CUDA_TAG) -- several GB ..."
    "$PY" -m pip install --no-cache-dir torch --index-url "https://download.pytorch.org/whl/$CUDA_TAG" \
        || fail "PyTorch installation failed (is $CUDA_TAG available?)."
fi
info "Installing numpy, tifffile ..."
"$PY" -m pip install --no-cache-dir numpy tifffile --quiet || fail "numpy/tifffile installation failed."

# ------------------------------------------------------------------ 3. script + launcher
cp "$SCRIPT_DIR/n2i_v3.py" "$OUT/n2i_v3.py"
cat > "$OUT/n2i" <<'EOF'
#!/usr/bin/env bash
# portable n2i launcher: uses only the Python inside this folder
DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
exec "$DIR/python/bin/python3" "$DIR/n2i_v3.py" "$@"
EOF
chmod +x "$OUT/n2i"

cat > "$OUT/check.sh" <<'EOF'
#!/usr/bin/env bash
# run on the target machine (on a GPU node) to check the bundle
DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
export PYTHONNOUSERSITE=1
unset PYTHONPATH PYTHONHOME
echo "glibc: $(ldd --version 2>/dev/null | head -n1)"
exec "$DIR/python/bin/python3" - <<'PY'
import sys
import numpy, tifffile, torch
print("python", sys.version.split()[0], "| torch", torch.__version__, "| CUDA build", torch.version.cuda,
      "| numpy", numpy.__version__, "| tifffile", tifffile.__version__)
if not torch.cuda.is_available():
    sys.exit("CUDA not available (GPU node? driver too old for this CUDA build?)")
x = torch.randn(2, 1, 64, 64, device="cuda")
torch.nn.functional.conv2d(x, torch.randn(8, 1, 3, 3, device="cuda")).sum().item()
print("GPU OK:", torch.cuda.get_device_name(0),
      "| training precision (auto):", "bf16" if torch.cuda.get_device_capability(0)[0] >= 8 else "fp16")
PY
EOF
chmod +x "$OUT/check.sh"

# ------------------------------------------------------------------ 4. test + archive
info "Testing imports ..."
if [ "${N2I_SKIP_TORCH:-0}" != "1" ]; then
    "$PY" -c "import torch, numpy, tifffile; print('torch', torch.__version__, 'CUDA', torch.version.cuda)" \
        || fail "import test failed."
fi
"$OUT/n2i" --help >/dev/null 2>&1 || [ "${N2I_SKIP_TORCH:-0}" = "1" ] || fail "launcher test failed."

info "Creating $OUT.tar.gz ..."
if command -v pigz >/dev/null; then
    tar -cf - "$OUT" | pigz > "$OUT.tar.gz"
else
    tar -czf "$OUT.tar.gz" "$OUT"
fi
info "Done: $OUT.tar.gz ($(du -h "$OUT.tar.gz" | cut -f1))"
echo "  On the target machine:"
echo "    tar -xzf $OUT.tar.gz"
echo "    ./$OUT/check.sh              # on a GPU node"
echo "    ./$OUT/n2i train            # in a data directory with rec/0, rec/1 ..."
echo "    ./$OUT/n2i infer"
