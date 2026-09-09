# Installing on a Windows PC and running the first reconstruction

2026-08-27

The steps from downloading the release zip to reconstructing one slice, for the
case where the bundled Windows x64 binaries (`exe/`) are used as they are. To
build from source instead, see "Binaries and building" in README.md.

---

## 1. Prerequisites

OS / hardware: 64-bit Windows 10 or 11. Every executable in `exe/` is native
x64.

CPU builds (programs ending in `_t`): the Microsoft Visual C++ 2015-2022
Redistributable (x64). The C runtime is linked statically, so normally nothing
is needed, but the 33 programs that use OpenMP (`ct_rec_t_*`, `hp_tg_t_*` and
so on) require `vcomp140.dll`, which the redistributable provides.

GPU builds (programs ending in `_g`): three things are needed.

| Requirement | Comes from | Applies to |
|---|---|---|
| NVIDIA driver (`nvcuda.dll`) | the driver installation | all 28 GPU programs |
| `cufft64_12.dll` | CUDA Toolkit (12.x / 13.x) | the 18 programs `ct_rec_g_*`, `hp_tg_g_*`, `ofct_rec_g_*`, `ofct_srec_g_*`, `p_rec_g_*`, `rec2rec_g_*` |
| a Turing (sm_75) or newer GPU | - | the bundled executables are built with CUDA Toolkit 13.2 |

The driver alone is not enough: the 18 programs that use cuFFT need the CUDA
Toolkit installed. A PC without a GPU can still do everything through the CPU
builds (`_t`).

> Only `ct_rec_tif_g_c` / `_r` / `_s` are left over from an older CUDA 11 build
> and require `cudart64_110.dll` and `cufft64_10.dll`. They were merged into
> `ct_rec` in ver 2.1, so use `ct_rec_g_*` instead.

Viewing the images (optional): the output is 32-bit float TIFF, so ImageJ or
Fiji is convenient. Putting the bundled `imagej-plugins/SP8CT_Plugins.jar` into
`plugins/` also lets you open HIS/IMG directly. See imagej-plugins/readme.txt.

## 2. Download and unpack

1. Get `ct-rec_v2.4.zip` from https://github.com/xrm-bl/ct-rec/releases/latest
2. Unblock the zip before extracting it: right-click the zip -> Properties ->
   at the bottom of the General tab tick "Security: Unblock" -> OK. Otherwise
   the extracted executables keep the "downloaded from the Web" mark and may be
   stopped by SmartScreen.
3. Extract it. Everything is inside a single folder, `ct-rec_v2.4\`.

## 3. Where to put it

Use a path with no spaces and no non-ASCII characters. The batch files move
around with pushd / popd, so avoid deep paths and folders under a sync client
(OneDrive and the like).

```
C:\ct-rec\ct-rec_v2.4\exe\      <- executables and batch files
C:\ct-rec\ct-rec_v2.4\src\      <- sources (ignore them if you do not build)
C:\ct-rec\ct-rec_v2.4\readme.txt
```

To update, extract the new zip into a separate folder and repoint PATH at it;
that way versions can coexist.

## 4. Putting the exe folder on PATH

Doing it through the GUI is the reliable way.

1. `Win + R` -> `rundll32 sysdm.cpl,EditEnvironmentVariables` -> Enter
2. In the upper pane (user variables), select `Path` and click Edit
3. Click New, add `C:\ct-rec\ct-rec_v2.4\exe`, then OK
4. Close and reopen any command prompt (running processes do not pick it up)

`setx PATH "%PATH%;C:\ct-rec\ct-rec_v2.4\exe"` also works, but setx truncates at
1024 characters and can damage an existing PATH. The GUI is recommended.

To leave PATH alone, put an `env.bat` like this in the working folder and run it
at the start of a session.

```bat
@echo off
set PATH=C:\ct-rec\ct-rec_v2.4\exe;%PATH%
set KERNEL_SIZE=11
set OMP_NUM_THREADS=16
cmd /k
```

## 5. Checking that it works

Open a command prompt and run a program with no arguments. The usage message
means the installation is good.

```
> ct_rec_t_c
parameter was wrong!!!
usage : ...\ct_rec_t_c.exe layer (center) (pixel size) (offsetangle)
default pixel size 1.0um
```

Then check a GPU build.

```
> ct_rec_g_c
```

The same usage message means the GPU dependencies are in place. If the window
closes immediately, or Windows reports that `cufft64_12.dll` is missing, the
CUDA Toolkit is not installed. `nvidia-smi` is a good way to confirm the driver
and the GPU generation as well.

## 6. Environment variables

Only the ones worth touching are listed here; the full list is in
20260723_CT_env_vars.md.

| Variable | Default | Notes |
|---|---|---|
| `KERNEL_SIZE` | 5 | strength of the ring removal; 0 or 1 disables it; around 11-17 for a sample with no prominent structure |
| `OMP_NUM_THREADS` | 40 | threads for the ring removal; it is 40 when unset, so lower it to about the number of logical cores on a small PC |
| `CBP_THREADS` | logical cores - 1 | threads for the CPU back projection; normally leave it alone |
| `CUDA_GPU` | 0 | which GPU to use when there are several |

`set KERNEL_SIZE=11` lasts for that command prompt only; for a permanent setting
add it as a user variable in the same GUI as step 4.

## 7. The first reconstruction

Preparing the data: the projections `q0001.img …` (or `q0001.tif …`) and the
dark frame `dark.img` (or `dark.tif`) must be in the same directory. The input
format is detected from the extension of the dark file. If the data came from
the beamline as `a.his` + `conv.bat` + `output.log`, run `conv.bat` in that
directory to split and convert it.

One slice:

```
> cd D:\data\001\raw
> ct_rec_g_c 120
```

`120` is the layer (height) to reconstruct. Omitting the rotation-axis position
makes the program estimate it. The result is `rec00120.tif`, a 32-bit float
TIFF. To give the axis and the pixel size explicitly, continue as
`ct_rec_g_c 120 1024.5 5.64`.

Checking: `pid rec00120.tif` prints the pixel size, rotation-axis position,
number of projections and min/max values embedded in the TIFF tags. Open the
image in ImageJ as well and look for rings and cupping.

The whole volume:

```
> mkdir rec
> hp_tg_g_c raw 5.64 1024.5 0 rec
```

The arguments are, in order: the directory holding the projections, the pixel
size (um), the rotation-axis position, the rotation-angle offset, and the output
directory. The output directory must exist beforehand.

Converting to 8/16-bit:

```
> mkdir ro_xy
> tif_f2i 8 rec ro_xy -0.5 3.0
```

## 8. Batch processing of continuously acquired data

The procedure for processing several measurements at once is in
rec-all_memo.md: `rc-check.bat 2101 4` for trial reconstructions, visual
inspection of the images in `rc-check\`, then `gen-all center.log` for the real
run (`act_rc-check.bat` / `act_gen-all.bat` for automatic CT).

`gen-all.bat` has the pixel size and similar values written into it, so copy it
into the working folder and edit the copy rather than the one in the exe folder;
that way an update does not overwrite your edits.

## 9. Troubleshooting

| Symptom | Cause and remedy |
|---|---|
| `'ct_rec_g_c' is not recognized as an internal or external command` | PATH is not set, or the command prompt was not reopened |
| only the GPU builds exit without printing anything | `cufft64_12.dll` is missing; install the CUDA Toolkit |
| `VCOMP140.DLL was not found` | install the VC++ 2015-2022 Redistributable (x64) |
| `CUDA error sort_filter_g.cu:266: an illegal memory access` | fixed in ver 2.4; check that you are not running an old executable |
| low-transmission warnings | the sample is too thick or the exposure too short; see 20260806_low_transmission_guard.md |
| complaint that the output directory does not exist | the hp_tg family does not create it; mkdir first |
