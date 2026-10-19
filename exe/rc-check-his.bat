@echo off
REM rc-check-his.bat - rotation-axis check of a Z scan directly from a.his
REM
REM   rc-check-his {layer}
REM     layer : layer reconstructed for the check (default: the centre row
REM             of the frame, Ny/2, as his_info reports it)
REM
REM Run in the directory holding a.his, conv.bat and output.log.  his_info
REM reads the his header, conv.bat and output.log and reports the number of
REM scans M in the file (frames / frames per scan), the centre layer and the
REM scan kind (180 / 360).  For every scan k the environment variable
REM RHIS_SCAN=k selects the frames of that scan (no spl / his_spl_E, no
REM conv.bat run, no .img files), his_rec reconstructs the layer with its
REM own axis estimate, the ImageDescription (pid) goes to center.log and the
REM image to rc-check\00k.tif for a visual check.
REM
REM center.log holds, in column 2, the value to pass to his_srec: for a
REM 180-degree scan that is the ImageDescription value itself (ct_rec and
REM hp_tg share the convention); for a 360-degree offset scan the
REM ImageDescription of his_rec (ofct_rec) carries the image centre
REM dcnt+0.5 while his_srec (ofct_srec) wants the detector axis position
REM Rc = dcnt, so 0.5 is subtracted here.  gen-all-his.bat reads center.log.
REM The program can be changed with HIS_REC (default his_rec_g_c).
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
for /f "tokens=1-7" %%a in ('his_info') do (
    set "NX=%%a"
    set "NY=%%b"
    set "NF=%%c"
    set "L=%%d"
    set "N=%%e"
    set "LAYER0=%%f"
    set "KIND=%%g"
)
if "%N%"=="" (
    echo error: his_info failed
    exit /b 1
)
set "LAYER=%~1"
if "%LAYER%"=="" set "LAYER=%LAYER0%"
echo %NX% x %NY%, %NF% frames: %N% scan^(s^) of %L% frames; %KIND%-degree scan; check layer %LAYER%
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
    REM ImageDescription: Dr centre views RA0 min max; column 2 -> the his_srec value
    set "LINE="
    for /f "tokens=1-6" %%a in ('pid %REC%') do (
        set "D1=%%a"
        set "D2=%%b"
        set "D3=%%c"
        set "D4=%%d"
        set "D5=%%e"
        set "D6=%%f"
    )
    if "%KIND%"=="360" (
        for /f %%v in ('powershell -NoProfile -Command "[math]::Round(!D2! - 0.5, 3)"') do set "D2=%%v"
    )
    echo !D1! !D2! !D3! !D4! !D5! !D6!>> center.log
    echo axis for his_srec: !D2!
    move /y %REC% rc-check\!DIR!.tif > nul
)
echo done: center.log ^(%N% lines^) and rc-check\*.tif written; next: gen-all-his center.log Dr
endlocal
exit /b 0

:usage
echo usage: %~nx0 {layer}
echo   layer : layer used for the check ^(default: centre row of the frame^)
echo   run where a.his, conv.bat and output.log are; the number of scans and the
echo   scan kind come from his_info; writes center.log ^(column 2 = value for
echo   his_srec^) and rc-check\00k.tif
exit /b 1
