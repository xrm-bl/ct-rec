@echo off
rem # 33-169
if "%~1"=="" (
    echo %0 stack_num
    exit /b
)

setlocal enabledelayedexpansion

set stepnum=%1

mkdir rc-check

for /l %%a in (1,1,%stepnum%) do (
  set num=000%%a
  set num=!num:~-3!
  
  cd !num!
rem    copy ..\output.log .
rem     copy ..\conv.bat .
rem    mkdir raw
rem     move *.tif raw
rem     move output.log raw
rem     move conv.bat raw
    cd raw
    call conv.bat
    tf_rec_g_r 70
    pid rec00070.tif >> ..\..\center.log
    copy rec00070.tif ..\..\rc-check\!num!.tif
    
    cd ..\..
)

rem move *.tif rc-check

endlocal



