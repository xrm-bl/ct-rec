@echo off
setlocal enabledelayedexpansion
rem ----------------------------------------------------------------------
rem ct_rec_loop.bat
rem
rem Reproduce hp_tg's whole-volume reconstruction by calling the per-slice
rem program ct_rec once per layer, running several jobs in parallel.
rem Arguments follow hp_tg; the LAST argument is the number of parallel
rem jobs (optional, default 1).
rem
rem   ct_rec_loop.bat  indir  Dr  RC  RA0  outdir  [Njobs]
rem   ct_rec_loop.bat  indir  Dr  L1 C1 L2 C2  RA0  outdir  [Njobs]
rem
rem   indir  : directory holding dark.img, q*.img, output.log
rem   Dr     : pixel size [um]
rem   RC     : rotation center (fixed for every layer)
rem   L1 C1  : center C1 at layer L1   (linear-center form)
rem   L2 C2  : center C2 at layer L2 -> center(z)=C1+(z-L1)(C2-C1)/(L2-L1)
rem   RA0    : offset angle [deg]
rem   outdir : directory to collect rec*.tif into
rem   Njobs  : number of concurrent jobs (default 1; tune to GPU memory)
rem
rem Run this from the directory ONE LEVEL ABOVE indir.
rem ----------------------------------------------------------------------

rem ct_rec executable name (assumed to be on PATH; a full path also works)
set "CT_EXE=ct_rec_g_c.exe"

rem ---- parse arguments (hp_tg style + optional trailing Njobs) ----
if not "%~8"=="" goto :linear
if not "%~5"=="" goto :fixed
goto :usage

:fixed
set "MODE=FIXED"
set "INDIR=%~1"
set "DR=%~2"
set "C1=%~3"
set "RA0=%~4"
set "OUTDIR=%~5"
set "NJ=%~6"
goto :parsed

:linear
set "MODE=LINEAR"
set "INDIR=%~1"
set "DR=%~2"
set "Z1=%~3"
set "C1=%~4"
set "Z2=%~5"
set "C2=%~6"
set "RA0=%~7"
set "OUTDIR=%~8"
set "NJ=%~9"
goto :parsed

:usage
echo usage: %~nx0 indir Dr RC RA0 outdir [Njobs]
echo        %~nx0 indir Dr L1 C1 L2 C2 RA0 outdir [Njobs]
exit /b 1

:parsed
if "%NJ%"=="" set "NJ=1"

rem image name used by the parallel-job counter (tasklist)
for %%E in ("%CT_EXE%") do set "CT_IMG=%%~nxE"

rem absolute paths (workers cd into INDIR, so OUTDIR/log must be absolute)
for %%I in ("%INDIR%")  do set "INDIR_ABS=%%~fI"
for %%O in ("%OUTDIR%") do set "OUTDIR_ABS=%%~fO"
set "LOGDIR_ABS=%OUTDIR_ABS%\log"

if not exist "%INDIR_ABS%\dark.img" ( echo no dark.img in "%INDIR_ABS%" & exit /b 1 )
if not exist "%OUTDIR_ABS%" mkdir "%OUTDIR_ABS%"
if not exist "%LOGDIR_ABS%" mkdir "%LOGDIR_ABS%"

rem fixed-center form reconstructs every layer (0..height-1);
rem height = 16-bit value at byte offset 6 of the HiPic .img header
rem (kept out of an if(...) block: the PowerShell parens would break it)
if /i not "%MODE%"=="FIXED" goto :afterheight
set "HEIGHT="
for /f "usebackq" %%H in (`powershell -NoProfile -Command "$f=[IO.File]::OpenRead('%INDIR_ABS%\dark.img');$b=New-Object byte[] 8;[void]$f.Read($b,0,8);$f.Close();[BitConverter]::ToUInt16($b,6)"`) do set "HEIGHT=%%H"
if not defined HEIGHT (
    echo cannot read height from dark.img
    exit /b 1
)
set /a Z1=0
set /a Z2=HEIGHT-1
set "C2=%C1%"
:afterheight

rem --- build the job list "z zp center" with a single PowerShell call ---
set "JOBLIST=%TEMP%\ctjobs_%RANDOM%%RANDOM%.txt"
powershell -NoProfile -Command "$z1=%Z1%;$z2=%Z2%;$c1=[double]'%C1%';$c2=[double]'%C2%';$m='%MODE%';for($z=$z1;$z -le $z2;$z++){if($m -eq 'LINEAR' -and $z2 -ne $z1){$c=$c1+($z-$z1)*($c2-$c1)/($z2-$z1)}else{$c=$c1};'{0} {1:00000} {2:0.0000}' -f $z,$z,$c}" > "%JOBLIST%"
if not exist "%JOBLIST%" ( echo failed to build job list & exit /b 1 )

rem --- worker: run ct_rec for one layer in INDIR, then move the result ---
set "WORKER=%TEMP%\ctjob_%RANDOM%%RANDOM%.bat"
> "%WORKER%" echo @echo off
>>"%WORKER%" echo cd /d "%%INDIR_ABS%%"
>>"%WORKER%" echo "%%CT_EXE%%" %%~1 %%~3 %%DR%% %%RA0%% ^> "%%LOGDIR_ABS%%\rec%%~2.log" 2^>^&1
>>"%WORKER%" echo move /y "rec%%~2.tif" "%%OUTDIR_ABS%%\" ^>nul 2^>^&1

echo ct_rec_loop: %MODE%  layers %Z1%..%Z2%  Njobs=%NJ%  indir="%INDIR_ABS%"  out="%OUTDIR_ABS%"

for /f "usebackq tokens=1,2,3" %%a in ("%JOBLIST%") do (
    call :throttle
    start "" /b cmd /c call "%WORKER%" %%a %%b %%c
)

call :waitall
del "%JOBLIST%" "%WORKER%" >nul 2>&1
echo done -^> "%OUTDIR_ABS%\"
endlocal
exit /b 0

rem --- block until fewer than NJ ct_rec processes are running ---
:throttle
for /f %%n in ('tasklist /nh /fi "imagename eq %CT_IMG%" 2^>nul ^| find /c /i "%CT_IMG%"') do set "RUN=%%n"
if !RUN! geq %NJ% (
    timeout /t 1 /nobreak >nul
    goto :throttle
)
exit /b

rem --- block until all ct_rec processes have finished ---
:waitall
for /f %%n in ('tasklist /nh /fi "imagename eq %CT_IMG%" 2^>nul ^| find /c /i "%CT_IMG%"') do set "RUN=%%n"
if !RUN! gtr 0 (
    timeout /t 1 /nobreak >nul
    goto :waitall
)
exit /b
