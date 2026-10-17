@echo off
REM rc-check-his.bat - rotation-axis check of a Z scan directly from a.his
REM
REM   rc-check-his {layer}
REM     layer : layer reconstructed for the check (default: the centre row
REM             of the frame, Ny/2, as his_info reports it)
REM
REM Run in the directory holding a.his, conv.bat and output.log.  his_info
REM reads the his header and conv.bat and reports the number of scans M in
REM the file (frames / frames per scan) and the centre layer.  For every
REM scan k the environment variable RHIS_SCAN=k selects the frames of that
REM scan (no spl / his_spl_E, no conv.bat run, no .img files), his_rec
REM reconstructs the layer with its own axis estimate, pid appends the
REM ImageDescription (column 2 = axis) to center.log and the image goes to
REM rc-check\00k.tif for a visual check.  center.log is the input of
REM gen-all-his.bat.  The program can be changed with HIS_REC (default
REM his_rec_g_c); his_rec picks ct_rec or ofct_rec from the angles in
REM output.log, so 180-degree and 360-degree offset Z scans both work.
if "%~1"=="-h" goto usage
if "%~1"=="/?" goto usage
if not "%~2"=="" goto usage
setlocal enabledelayedexpansion
if "%HIS_REC%"=="" set "HIS_REC=his_rec_g_c"
if not exist conv.bat (
    echo error: conv.bat not found here ^(run beside a.his, conv.bat and output.log^)
    exit /b 1
)
set "N="
for /f "tokens=1-6" %%a in ('his_info') do (
    set "NX=%%a"
    set "NY=%%b"
    set "NF=%%c"
    set "L=%%d"
    set "N=%%e"
    set "LAYER0=%%f"
)
if "%N%"=="" (
    echo error: his_info failed
    exit /b 1
)
set "LAYER=%~1"
if "%LAYER%"=="" set "LAYER=%LAYER0%"
echo %NX% x %NY%, %NF% frames: %N% scan^(s^) of %L% frames; check layer %LAYER%
set "REC=00000%LAYER%"
set "REC=rec!REC:~-5!.tif"
if exist center.log del center.log
if not exist rc-check mkdir rc-check

for /l %%i in (1,1,%N%) do (
    set "DIR=00%%i"
    set "DIR=!DIR:~-3!"
    echo ==== scan !DIR! / %N% ====
    set "RHIS_SCAN=%%i"
    %HIS_REC% %LAYER%
    if errorlevel 1 (
        echo error: %HIS_REC% failed on scan %%i
        exit /b 1
    )
    pid %REC% >> center.log
    move /y %REC% rc-check\!DIR!.tif > nul
)
echo done: center.log ^(%N% lines^) and rc-check\*.tif written; next: gen-all-his center.log Dr
endlocal
exit /b 0

:usage
echo usage: %~nx0 {layer}
echo   layer : layer used for the check ^(default: centre row of the frame^)
echo   run where a.his, conv.bat and output.log are; the number of scans comes
echo   from his_info; writes center.log and rc-check\00k.tif
exit /b 1
