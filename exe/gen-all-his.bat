@echo off
REM gen-all-his.bat - whole-volume reconstruction of a Z scan directly from a.his
REM
REM   gen-all-his center.log Dr {LACmin LACmax}
REM     center.log : written by rc-check-his (one line per scan, column 2 = axis)
REM     Dr         : pixel size
REM     LACmin/max : range of the 8-bit conversion (default -0.5 3.0)
REM
REM Run in the directory holding a.his, conv.bat and output.log (the same
REM place as rc-check-his).  Line k of center.log belongs to scan k: with
REM RHIS_SCAN=k his_srec reconstructs that scan into 00k\rec and tif_f2i
REM converts it to 8-bit in 00k\ro_xy, as gen-all.bat does after spl +
REM conv.bat + hp_tg (the ro_zx reslice of gen-all.bat is not made).  The
REM directories 00k\rec are what rec_stk expects.  The program can be
REM changed with HIS_SREC (default his_srec_g_c); a 360-degree offset scan
REM is handled by his_srec itself (ofct_srec).
if "%~2"=="" (
    echo usage: %~nx0 center.log Dr {LACmin LACmax}
    echo   center.log : from rc-check-his, Dr : pixel size, LACmin LACmax : 8-bit range ^(default -0.5 3.0^)
    echo   run where a.his, conv.bat and output.log are; writes 00k\rec and 00k\ro_xy
    exit /b 1
)
setlocal enabledelayedexpansion
set "LOGFILE=%~1"
set "DR=%~2"
set "LMIN=%~3"
set "LMAX=%~4"
if "%LMIN%"=="" set "LMIN=-0.5"
if "%LMAX%"=="" set "LMAX=3.0"
if "%HIS_SREC%"=="" set "HIS_SREC=his_srec_g_c"
if not exist "%LOGFILE%" (
    echo error: %LOGFILE% not found
    exit /b 1
)
set "LINE_NUM=0"

REM ==== center.log line by line (echo. guards a missing final newline) ====
for /f "usebackq tokens=1,2,3,4,5,6 delims=	 " %%a in (`type "%LOGFILE%" ^& echo.`) do (
    if not "%%a"=="" (
        set /a LINE_NUM+=1
        set "DIR=00!LINE_NUM!"
        set "DIR=!DIR:~-3!"
        echo ==== scan !DIR! ^(center %%b^) ====
        set "RHIS_SCAN=!LINE_NUM!"
        if not exist !DIR! mkdir !DIR!
        if not exist !DIR!\rec mkdir !DIR!\rec
        if not exist !DIR!\ro_xy mkdir !DIR!\ro_xy
        %HIS_SREC% . %DR% %%b 0 !DIR!\rec
        if errorlevel 1 (
            echo error: %HIS_SREC% failed on scan !LINE_NUM!
            exit /b 1
        )
        pushd !DIR!
        tif_f2i 8 rec ro_xy %LMIN% %LMAX%
        popd
    )
)
echo done: !LINE_NUM! scans; next: rec_stk !LINE_NUM! start end ^> ddd.bat
endlocal
