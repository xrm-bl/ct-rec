@echo off
rem Runs setup_n2i.ps1 without changing the system execution policy.
rem Options are passed through, e.g.:  setup_n2i.bat -InstallDir D:\n2i -InstallPython
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0setup_n2i.ps1" %*
pause
