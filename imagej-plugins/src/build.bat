@echo off
setlocal

REM ============================================================
REM  build.bat - SP8CT ImageJ plugin build script (Windows)
REM
REM  Usage:  build.bat [path\to\ij.jar]
REM
REM  Run from the src directory where all .java files reside.
REM  If ij.jar path is not specified, it tries common locations.
REM ============================================================

set JARNAME=SP8CT_Plugins.jar
set BUILDDIR=build

REM --- Locate ij.jar ---
if not "%~1"=="" (
    set IJJAR=%~1
) else if exist "..\ij.jar" (
    set IJJAR=..\ij.jar
) else if exist "..\..\ij.jar" (
    set IJJAR=..\..\ij.jar
) else if defined IMAGEJ_DIR (
    set IJJAR=%IMAGEJ_DIR%\ij.jar
) else (
    echo ERROR: ij.jar not found.
    echo Usage: build.bat path\to\ij.jar
    exit /b 1
)

echo Using ij.jar: %IJJAR%

REM --- Clean and create build directory ---
if exist %BUILDDIR% rmdir /s /q %BUILDDIR%
mkdir %BUILDDIR%

REM --- Compile all Java files ---
echo Compiling...
javac -cp "%IJJAR%" -d %BUILDDIR% -encoding UTF-8 -source 8 -target 8 ^
    Open_HIS_IMG.java ^
    HandleExtraFileTypes.java ^
    Gaussian_Filter_2D.java ^
    Gaussian_Filter_3D.java ^
    Gaussian_Filter_3D_Ext.java ^
    Median_Filter_2D.java ^
    Median_Filter_3D_Ext.java ^
    Bilateral_Filter_2D.java ^
    Bilateral_Filter_3D_Ext.java ^
    NLM_Filter_2D.java ^
    NLM_Filter_3D_Ext.java ^
    TV_Denoise_2D.java ^
    TV_Denoise_3D_Ext.java ^
    Wavelet_Denoise_2D.java ^
    Wavelet_Denoise_3D_Ext.java ^
    Anisotropic_Diffusion_2D.java ^
    Anisotropic_Diffusion_3D_Ext.java ^
    BM3D_Filter_2D.java ^
    BM4D_Filter_3D_Ext.java ^
    Stack_Crop_3D.java

if errorlevel 1 (
    echo.
    echo Compilation FAILED.
    exit /b 1
)

echo Compilation OK.

REM --- Copy plugins.config into build ---
copy plugins.config %BUILDDIR%\plugins.config >nul

REM --- Create JAR ---
echo Creating %JARNAME%...
cd %BUILDDIR%
jar cf ..\%JARNAME% .
cd ..

echo.
echo ============================================================
echo  Build complete: %JARNAME%
echo.
echo  Install:
echo    Copy %JARNAME% to ImageJ/plugins/ or Fiji.app/plugins/
echo    and restart ImageJ/Fiji.
echo.
echo  D^&D support (optional):
echo    Also copy HandleExtraFileTypes.class from
echo    %BUILDDIR%\ to your plugins/ folder (ImageJ only).
echo    For Fiji, see DragAndDrop_Setup.txt.
echo ============================================================

endlocal
