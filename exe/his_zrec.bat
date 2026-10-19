@echo off
REM his_zrec.bat - one slice of one scan of a Z scan, kept in rc-check\
REM
REM   his_zrec scan layer {center}
REM     scan   : scan number in a.his (1..M, see his_info)
REM     layer  : layer to reconstruct
REM     center : rotation axis position (optional; his_rec estimates it
REM              from that layer when omitted)
REM
REM Run in the directory holding a.his, conv.bat and output.log.  With
REM RHIS_SCAN=scan his_rec reconstructs the layer and the image is moved to
REM rc-check\SSS_LLLLL_C.tif (scan as %03d, layer as %05d, C = the axis
REM actually used, as the value to pass to his_srec: the ImageDescription
REM read back with pid, minus 0.5 for a 360-degree offset scan where
REM ofct_rec records the image centre dcnt+0.5; trailing zeros removed,
REM e.g. rc-check\003_01024_1023.5.tif), so that trial reconstructions of
REM several layers and axis positions can be compared side by side.  The
REM program can be changed with HIS_REC (default his_rec_g_c).
if "%~2"=="" goto usage
if not "%~4"=="" goto usage
setlocal enabledelayedexpansion
set "SCAN=%~1"
set "LAYER=%~2"
set "CENTER=%~3"
if "%HIS_REC%"=="" set "HIS_REC=his_rec_g_c"
if not exist conv.bat (
    echo error: conv.bat not found here ^(run beside a.his, conv.bat and output.log^)
    exit /b 1
)
set "KIND="
for /f "tokens=7" %%a in ('his_info') do set "KIND=%%a"
if "%KIND%"=="" (
    echo error: his_info failed
    exit /b 1
)
set "SSS=00%SCAN%"
set "SSS=%SSS:~-3%"
set "LLLLL=0000%LAYER%"
set "LLLLL=%LLLLL:~-5%"
set "REC=rec%LLLLL%.tif"
if not exist rc-check mkdir rc-check

set "RHIS_SCAN=%SCAN%"
echo ==== scan %SSS%, layer %LLLLL% ^(%KIND%-degree scan^) ====
%HIS_REC% %LAYER% %CENTER%
if errorlevel 1 (
    echo error: %HIS_REC% failed
    exit /b 1
)
if not exist %REC% (
    echo error: %REC% was not written
    exit /b 1
)

REM axis actually used: column 2 of the ImageDescription (pid); -0.5 for an offset scan; trailing zeros trimmed
set "C="
for /f "tokens=2" %%a in ('pid %REC%') do set "C=%%a"
if "%C%"=="" set "C=%CENTER%"
if "%C%"=="" set "C=auto"
if "%KIND%"=="360" if not "%C%"=="auto" (
    for /f %%v in ('powershell -NoProfile -Command "[math]::Round(!C! - 0.5, 3)"') do set "C=%%v"
)
echo !C! | find "." > nul
if errorlevel 1 goto named
:trim
if "!C:~-1!"=="0" (
    set "C=!C:~0,-1!"
    goto trim
)
if "!C:~-1!"=="." set "C=!C:~0,-1!"
:named

set "OUT=rc-check\%SSS%_%LLLLL%_%C%.tif"
move /y %REC% "%OUT%" > nul
echo %OUT%
endlocal
exit /b 0

:usage
echo usage: %~nx0 scan layer {center}
echo   scan   : scan number in a.his ^(1..M, see his_info^)
echo   layer  : layer to reconstruct
echo   center : rotation axis position ^(optional; estimated when omitted^)
echo   run where a.his, conv.bat and output.log are; writes rc-check\SSS_LLLLL_C.tif
exit /b 1
