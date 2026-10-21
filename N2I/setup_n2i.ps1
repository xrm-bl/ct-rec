<#
.SYNOPSIS
  Sets up a Python venv with CUDA-enabled PyTorch for n2i_v2.py (Noise2Inverse) on Windows.

.DESCRIPTION
  1. Checks the NVIDIA driver (nvidia-smi) and the highest CUDA version it supports.
  2. Finds Python (py launcher / PATH); optionally installs it with winget.
  3. Creates <InstallDir>\venv.
  4. Installs PyTorch, trying the newest CUDA build the driver supports first and
     falling back to older builds if the GPU test fails (e.g. older GPU generations).
  5. Installs numpy and tifffile, copies n2i_v2.py, and writes the launcher n2i.bat.
  Safe to re-run: an existing working installation is detected and kept.

.EXAMPLE
  powershell -ExecutionPolicy Bypass -File setup_n2i.ps1
  powershell -ExecutionPolicy Bypass -File setup_n2i.ps1 -InstallDir D:\n2i
  powershell -ExecutionPolicy Bypass -File setup_n2i.ps1 -CudaTag cu126
  powershell -ExecutionPolicy Bypass -File setup_n2i.ps1 -InstallPython
  powershell -ExecutionPolicy Bypass -File setup_n2i.ps1 -Cpu
#>
param(
    [string]$InstallDir    = "C:\n2i",
    [string]$PythonVersion = "3.12",
    [string]$CudaTag       = "",      # force a PyTorch wheel tag, e.g. cu130 / cu126
    [switch]$InstallPython,           # install Python with winget if it is not found
    [switch]$Cpu,                     # CPU-only PyTorch (no NVIDIA GPU)
    [switch]$Reinstall                # reinstall PyTorch even if a working one exists
)

$ErrorActionPreference = "Continue"
function Info($m) { Write-Host "[n2i] $m" -ForegroundColor Cyan }
function Warn($m) { Write-Host "[n2i] $m" -ForegroundColor Yellow }
function Fail($m) { Write-Host "[n2i] ERROR: $m" -ForegroundColor Red; exit 1 }

# PyTorch wheel indexes, newest first. Edit this list when PyTorch adds/drops CUDA builds.
$TorchTags = @("cu132", "cu130", "cu128", "cu126", "cu124", "cu121", "cu118")

# ------------------------------------------------------------------ 1. GPU / driver
$driverCuda = 0
if (-not $Cpu) {
    if (-not (Get-Command nvidia-smi -ErrorAction SilentlyContinue)) {
        Fail "nvidia-smi not found. Install the latest NVIDIA driver, or rerun with -Cpu."
    }
    $smiOut = (& nvidia-smi) | Out-String
    if ($smiOut -match "CUDA Version:\s*([0-9]+)\.([0-9]+)") {
        $driverCudaStr = "$($Matches[1]).$($Matches[2])"
        $driverCuda = [int]$Matches[1] * 10 + [int]$Matches[2]
    } else {
        Fail "Could not read the CUDA version from nvidia-smi."
    }
    $gpus = (& nvidia-smi --query-gpu=name,memory.total --format=csv,noheader) -join "; "
    Info "GPU: $gpus"
    Info "Driver supports CUDA up to $driverCudaStr"
}

# ------------------------------------------------------------------ 2. Python
function Test-Python($exe, $pre) {
    $code = "import sys,struct; ok = sys.version_info >= (3,10) and struct.calcsize('P') == 8; print(sys.version.split()[0]); sys.exit(0 if ok else 3)"
    $v = & $exe @pre -c $code 2>$null
    if ($LASTEXITCODE -eq 0) { return $v }
    return $null
}

function Find-Python {
    if (Get-Command py -ErrorAction SilentlyContinue) {
        foreach ($sel in @("-$PythonVersion", "-3")) {
            $v = Test-Python "py" @($sel)
            if ($v) { return @{ Exe = "py"; Pre = @($sel); Ver = $v } }
        }
    }
    $cmd = Get-Command python -ErrorAction SilentlyContinue
    if ($cmd -and $cmd.Source -notmatch "WindowsApps") {   # skip the Microsoft Store stub
        $v = Test-Python $cmd.Source @()
        if ($v) { return @{ Exe = $cmd.Source; Pre = @(); Ver = $v } }
    }
    return $null
}

$py = Find-Python
if (-not $py -and $InstallPython) {
    if (-not (Get-Command winget -ErrorAction SilentlyContinue)) {
        Fail "winget not available. Install Python $PythonVersion (64-bit) from python.org manually."
    }
    Info "Installing Python $PythonVersion with winget ..."
    & winget install -e --id "Python.Python.$PythonVersion" --scope user `
        --accept-package-agreements --accept-source-agreements
    $env:Path = [Environment]::GetEnvironmentVariable("Path", "Machine") + ";" +
                [Environment]::GetEnvironmentVariable("Path", "User")
    $py = Find-Python
}
if (-not $py) {
    Fail "64-bit Python >= 3.10 not found. Install Python $PythonVersion from python.org, or rerun with -InstallPython."
}
Info "Python: $($py.Exe) $($py.Pre -join ' ')  (version $($py.Ver))"

# ------------------------------------------------------------------ 3. venv
if ($InstallDir.Length -gt 60) {
    Warn "InstallDir is long; PyTorch has deep paths. A short path such as C:\n2i is safer."
}
New-Item -ItemType Directory -Force -Path $InstallDir | Out-Null
$venv = Join-Path $InstallDir "venv"
$vpy  = Join-Path $venv "Scripts\python.exe"
if (-not (Test-Path $vpy)) {
    Info "Creating venv: $venv"
    & $py.Exe @($py.Pre) -m venv $venv
    if ($LASTEXITCODE -ne 0 -or -not (Test-Path $vpy)) { Fail "venv creation failed." }
} else {
    Info "Using existing venv: $venv"
}
& $vpy -m pip install --upgrade pip --quiet
if ($LASTEXITCODE -ne 0) { Fail "pip upgrade failed (network / proxy?)." }

# ------------------------------------------------------------------ 4. PyTorch
$check = Join-Path $InstallDir "_torch_check.py"
@'
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
      "| bf16 (--amp):", "yes" if torch.cuda.is_bf16_supported() else "no")
'@ | Set-Content -Path $check -Encoding ASCII

$chkArgs = @($check)
if ($Cpu) { $chkArgs += "--cpu" }

$torchOk = $false
if (-not $Reinstall) {
    & $vpy @chkArgs 2>$null
    if ($LASTEXITCODE -eq 0) { $torchOk = $true; Info "Working PyTorch already installed; keeping it (use -Reinstall to change)." }
}

if (-not $torchOk) {
    if ($Cpu)          { $cands = @("cpu") }
    elseif ($CudaTag)  { $cands = @($CudaTag) }
    else               { $cands = @($TorchTags | Where-Object { [int]$_.Substring(2) -le $driverCuda }) }
    if ($cands.Count -eq 0) {
        Fail "Driver CUDA $driverCudaStr is older than every build in the list. Update the NVIDIA driver."
    }
    foreach ($t in $cands) {
        & $vpy -m pip uninstall -y torch 2>$null | Out-Null
        Info "Installing PyTorch ($t) -- this downloads ~2-3 GB ..."
        & $vpy -m pip install torch --index-url "https://download.pytorch.org/whl/$t"
        if ($LASTEXITCODE -ne 0) { Warn "No usable wheel for $t (or download failed); trying the next one."; continue }
        & $vpy @chkArgs
        if ($LASTEXITCODE -eq 0) { $torchOk = $true; Info "PyTorch $t works."; break }
        Warn "PyTorch $t installed but the GPU test failed; trying an older CUDA build."
    }
    if (-not $torchOk) { Fail "No PyTorch build passed the test. Check the driver / GPU, or try -CudaTag." }
}
Remove-Item $check -ErrorAction SilentlyContinue

# ------------------------------------------------------------------ 5. other packages, script, launcher
Info "Installing numpy, tifffile ..."
& $vpy -m pip install --upgrade numpy tifffile --quiet
if ($LASTEXITCODE -ne 0) { Fail "numpy/tifffile installation failed." }

$src = Join-Path $PSScriptRoot "n2i_v2.py"
$dst = Join-Path $InstallDir "n2i_v2.py"
if (Test-Path $src) {
    if ((Resolve-Path $src).Path -ne [IO.Path]::GetFullPath($dst)) { Copy-Item $src $dst -Force }
    Info "n2i_v2.py -> $dst"
} elseif (-not (Test-Path $dst)) {
    Warn "n2i_v2.py not found next to this script. Copy it to $InstallDir manually."
}

$bat = Join-Path $InstallDir "n2i.bat"
@'
@echo off
"%~dp0venv\Scripts\python.exe" "%~dp0n2i_v2.py" %*
'@ | Set-Content -Path $bat -Encoding ASCII

Write-Host ""
Info "Setup complete."
Write-Host "  Launcher : $bat"
Write-Host "  Example  :"
Write-Host "    $bat train --mode 25d --splits D:\ct\rec0 D:\ct\rec1 D:\ct\rec2 D:\ct\rec3 --out D:\ct\n2i_25d.pt --workers 4"
Write-Host "    $bat infer --splits D:\ct\rec0 D:\ct\rec1 D:\ct\rec2 D:\ct\rec3 --model D:\ct\n2i_25d.pt --outdir D:\ct\denoised --fp16"