@echo off
REM ==== 引数チェック ====
if "%~2"=="" (
    echo 使用方法: %0 撮影枚数 測定回数
    echo 例: %0 2101 4
    exit /b
)

setlocal enabledelayedexpansion
set "N=%2"
set "BASEDIR=%CD%"

spl  %1 %2 > aaa.bat
call aaa.bat
mkdir rc-check

REM ==== 001 ～ N までループ ====
for /l %%i in (1,1,%N%) do (
    set "DIR=00%%i"
    set "DIR=!DIR:~-3!"
    echo ==== Processing !DIR! ====
    
    if not exist "!DIR!\raw" (
        echo エラー: !DIR!\raw が見つかりません。スキップします。
    ) else (
        pushd "!DIR!\raw"
        call conv.bat
        ct_rec_g_c 120
        pid rec00120.tif >> "!BASEDIR!\center.log"
        copy rec00120.tif "!BASEDIR!\rc-check\!DIR!.tif"
        popd
    )
)

echo 完了しました。
endlocal
