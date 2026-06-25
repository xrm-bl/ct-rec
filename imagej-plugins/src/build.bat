@echo off
setlocal

REM ============================================================
REM  build.bat - SP8CT ImageJ plugin build script (Windows)
REM
REM  Usage:  build.bat [path\to\ij.jar]
REM
REM  Run from the src directory where all .java files reside.
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
javac -cp "%IJJAR%" -d %BUILDDIR% -encoding UTF-8 -source 8 -target 8 *.java

if errorlevel 1 (
    echo.
    echo Compilation FAILED.
    exit /b 1
)

echo Compilation OK.

REM --- Copy plugins.config into build ---
copy plugins.config %BUILDDIR%\plugins.config >nul

REM --- Create JAR ---
REM   Include all .class files EXCEPT HandleExtraFileTypes*
REM   (HandleExtraFileTypes must be placed as a standalone .class)
echo Creating %JARNAME%...
cd %BUILDDIR%

REM Build file list excluding HandleExtraFileTypes
(for %%f in (*.class) do (
    echo %%f | findstr /i /v "HandleExtraFileTypes" >nul && echo %%f
)) > _jarfiles.txt
echo plugins.config >> _jarfiles.txt

jar cf ..\%JARNAME% @_jarfiles.txt
del _jarfiles.txt
cd ..

echo.
echo ============================================================
echo  Build complete: %JARNAME%
echo.
echo  Install:
echo    1. Copy %JARNAME% to ImageJ\plugins\
echo    2. Copy %BUILDDIR%\HandleExtraFileTypes.class
echo       to ImageJ\plugins\ (for D^&D support)
echo    3. Restart ImageJ
echo ============================================================

endlocal