#!/usr/bin/env bash
# setup_n2i.sh -- set up a Python venv with CUDA-enabled PyTorch for n2i_v3.py on Ubuntu.
#
#   1. checks the NVIDIA driver (nvidia-smi) and the highest CUDA version it supports
#   2. checks python3 (>= 3.10, 64-bit); optionally installs python3-venv with apt
#   3. creates <dir>/venv
#   4. installs PyTorch, newest CUDA build the driver supports first, falling back to
#      older builds if the GPU test fails (e.g. older GPU generations)
#   5. installs numpy + tifffile, copies n2i_v3.py, writes the launcher <dir>/n2i
# Safe to re-run: a working PyTorch is detected and kept.
#
# Usage:
#   bash setup_n2i.sh                      # install to ~/n2i
#   bash setup_n2i.sh --dir /data/n2i      # other location
#   bash setup_n2i.sh --cuda-tag cu126     # force a PyTorch CUDA build
#   bash setup_n2i.sh --install-deps       # apt install python3-venv etc. (uses sudo)
#   bash setup_n2i.sh --cpu                # CPU-only PyTorch
#   bash setup_n2i.sh --reinstall          # reinstall PyTorch
#   bash setup_n2i.sh --link               # also link the launcher to ~/.local/bin/n2i

set -u

INSTALL_DIR="$HOME/n2i"
CUDA_TAG=""
PYTHON="${PYTHON:-python3}"
CPU=0; REINSTALL=0; INSTALL_DEPS=0; LINK=0
# PyTorch wheel indexes, newest first. Edit when PyTorch adds/drops CUDA builds.
TORCH_TAGS=(cu132 cu130 cu128 cu126 cu124 cu121 cu118)

info() { printf '\033[36m[n2i]\033[0m %s\n' "$*"; }
warn() { printf '\033[33m[n2i]\033[0m %s\n' "$*"; }
fail() { printf '\033[31m[n2i] ERROR:\033[0m %s\n' "$*"; exit 1; }

usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; }

while [ $# -gt 0 ]; do
    case "$1" in
        --dir)          INSTALL_DIR="$2"; shift 2 ;;
        --cuda-tag)     CUDA_TAG="$2"; shift 2 ;;
        --python)       PYTHON="$2"; shift 2 ;;
        --cpu)          CPU=1; shift ;;
        --reinstall)    REINSTALL=1; shift ;;
        --install-deps) INSTALL_DEPS=1; shift ;;
        --link)         LINK=1; shift ;;
        -h|--help)      usage; exit 0 ;;
        *)              warn "unknown option: $1"; usage; exit 1 ;;
    esac
done

[ "$(id -u)" -eq 0 ] && warn "Running as root: the venv will be owned by root. A normal user is recommended."
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
mkdir -p "$INSTALL_DIR" || fail "cannot create $INSTALL_DIR"
INSTALL_DIR="$(cd "$INSTALL_DIR" && pwd)"

# ------------------------------------------------------------------ 1. GPU / driver
DRIVER_CUDA=0
if [ "$CPU" -eq 0 ]; then
    command -v nvidia-smi >/dev/null 2>&1 || \
        fail "nvidia-smi not found. Install the NVIDIA driver (e.g. 'sudo ubuntu-drivers install', then reboot), or use --cpu."
    ver=$(nvidia-smi 2>/dev/null | grep -oP 'CUDA Version:\s*\K[0-9]+\.[0-9]+' | head -n1)
    [ -n "$ver" ] || fail "could not read the CUDA version from nvidia-smi (driver loaded? try rebooting)."
    DRIVER_CUDA=$(( ${ver%%.*} * 10 + ${ver##*.} ))
    info "GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader | paste -sd ';' -)"
    info "Driver supports CUDA up to $ver"
fi

# ------------------------------------------------------------------ 2. Python
if [ "$INSTALL_DEPS" -eq 1 ]; then
    info "Installing python3, python3-venv, python3-pip with apt (sudo) ..."
    sudo apt-get update && sudo apt-get install -y python3 python3-venv python3-pip || fail "apt install failed."
fi
command -v "$PYTHON" >/dev/null 2>&1 || fail "$PYTHON not found. Install it ('sudo apt install python3') or use --install-deps."
PYVER=$("$PYTHON" -c 'import sys,struct; ok=sys.version_info>=(3,10) and struct.calcsize("P")==8; print("%d.%d"%sys.version_info[:2]); sys.exit(0 if ok else 3)') \
    || fail "64-bit Python >= 3.10 required (found $PYVER)."
info "Python: $(command -v "$PYTHON") ($PYVER)"

# ------------------------------------------------------------------ 3. venv
VENV="$INSTALL_DIR/venv"
VPY="$VENV/bin/python"
if [ ! -x "$VPY" ]; then
    info "Creating venv: $VENV"
    if ! "$PYTHON" -m venv "$VENV"; then
        rm -rf "$VENV"
        fail "venv creation failed. Install 'python3-venv' (or 'python$PYVER-venv'): sudo apt install python$PYVER-venv, or rerun with --install-deps."
    fi
else
    info "Using existing venv: $VENV"
fi
"$VPY" -m pip install --upgrade pip --quiet || fail "pip upgrade failed (network / proxy?)."

# ------------------------------------------------------------------ 4. PyTorch
CHECK="$INSTALL_DIR/_torch_check.py"
cat > "$CHECK" <<'EOF'
import sys
try:
    import torch
except Exception as e:
    print("torch import failed:", e); sys.exit(1)
print("torch", torch.__version__, "| built for CUDA", torch.version.cuda)
if "--cpu" in sys.argv:
    sys.exit(0)
if not torch.cuda.is_available():
    print("CUDA not available"); sys.exit(2)
try:
    x = torch.randn(2, 1, 64, 64, device="cuda")
    w = torch.randn(8, 1, 3, 3, device="cuda")
    y = torch.nn.functional.conv2d(x, w).sum()
    torch.cuda.synchronize()
except Exception as e:
    print("GPU kernel test failed:", e); sys.exit(3)
cap = torch.cuda.get_device_capability(0)
print("GPU OK:", torch.cuda.get_device_name(0), "| sm_%d%d" % cap,
      "| training precision (auto):", "bf16" if torch.cuda.get_device_capability(0)[0] >= 8 else "fp16")
EOF
CHK_ARGS=("$CHECK"); [ "$CPU" -eq 1 ] && CHK_ARGS+=("--cpu")

TORCH_OK=0
if [ "$REINSTALL" -eq 0 ] && "$VPY" "${CHK_ARGS[@]}" 2>/dev/null; then
    TORCH_OK=1
    info "Working PyTorch already installed; keeping it (use --reinstall to change)."
fi

if [ "$TORCH_OK" -eq 0 ]; then
    if [ "$CPU" -eq 1 ]; then CANDS=(cpu)
    elif [ -n "$CUDA_TAG" ]; then CANDS=("$CUDA_TAG")
    else
        CANDS=()
        for t in "${TORCH_TAGS[@]}"; do
            [ "${t#cu}" -le "$DRIVER_CUDA" ] && CANDS+=("$t")
        done
    fi
    [ "${#CANDS[@]}" -gt 0 ] || fail "driver CUDA $ver is older than every build in the list. Update the NVIDIA driver."
    for t in "${CANDS[@]}"; do
        "$VPY" -m pip uninstall -y torch >/dev/null 2>&1
        info "Installing PyTorch ($t) -- this downloads ~2-3 GB ..."
        if ! "$VPY" -m pip install torch --index-url "https://download.pytorch.org/whl/$t"; then
            warn "No usable wheel for $t (or download failed); trying the next one."
            continue
        fi
        if "$VPY" "${CHK_ARGS[@]}"; then
            TORCH_OK=1; info "PyTorch $t works."; break
        fi
        warn "PyTorch $t installed but the GPU test failed; trying an older CUDA build."
    done
    [ "$TORCH_OK" -eq 1 ] || fail "no PyTorch build passed the test. Check the driver / GPU, or try --cuda-tag."
fi
rm -f "$CHECK"

# ------------------------------------------------------------------ 5. packages, script, launcher
info "Installing numpy, tifffile ..."
"$VPY" -m pip install --upgrade numpy tifffile --quiet || fail "numpy/tifffile installation failed."

SRC="$SCRIPT_DIR/n2i_v3.py"
DST="$INSTALL_DIR/n2i_v3.py"
if [ -f "$SRC" ]; then
    [ "$(readlink -f "$SRC")" != "$(readlink -f "$DST")" ] && cp -f "$SRC" "$DST"
    info "n2i_v3.py -> $DST"
elif [ ! -f "$DST" ]; then
    warn "n2i_v3.py not found next to this script. Copy it to $INSTALL_DIR manually."
    ls "$SCRIPT_DIR"/n2i*.py 2>/dev/null | sed 's/^/        similar file: /'
fi

LAUNCHER="$INSTALL_DIR/n2i"
cat > "$LAUNCHER" <<'EOF'
#!/usr/bin/env bash
DIR="$(cd "$(dirname "$(readlink -f "$0")")" && pwd)"
exec "$DIR/venv/bin/python" "$DIR/n2i_v3.py" "$@"
EOF
chmod +x "$LAUNCHER"

if [ "$LINK" -eq 1 ]; then
    mkdir -p "$HOME/.local/bin"
    ln -sf "$LAUNCHER" "$HOME/.local/bin/n2i"
    info "Linked: ~/.local/bin/n2i (open a new terminal if 'n2i' is not found)"
fi

echo
info "Setup complete."
echo "  Launcher : $LAUNCHER"
echo "  Example (in a data directory with rec/0, rec/1 ...):"
echo "    $LAUNCHER train"
echo "    $LAUNCHER infer"
