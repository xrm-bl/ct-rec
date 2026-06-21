@echo off
REM ==== 引数チェック ====
if "%~1"=="" (
    echo 使用方法: %0 center.logのパス
    echo 例: %0 center.log
    exit /b
)

setlocal enabledelayedexpansion
set "LOGFILE=%~1"
set "BASEDIR=%CD%"
set "LINE_NUM=0"

REM ==== center.log を1行ずつ読み込み（末尾改行問題を回避） ====
for /f "usebackq tokens=1,2,3,4,5,6 delims=	 " %%a in (`type "%LOGFILE%" ^& echo.`) do (
    if not "%%a"=="" (
        set /a LINE_NUM+=1
        set "DIR=00!LINE_NUM!"
        set "DIR=!DIR:~-3!"
        set "COL2=%%b"
        
        echo "==== Processing !DIR! (center=!COL2!) ===="
        
        if not exist "!DIR!" (
            echo エラー: !DIR! が見つかりません。スキップします。
        ) else (
            pushd "!DIR!"
            mkdir rec ro_xy ro_zx 2>nul
            hp_tg_g_c raw 5.64 !COL2! 0 rec
            tif_f2i 8 rec ro_xy -0.5 3.0
            si_rar.exe ro_xy - +x +z +y ro_zx
            popd
        )
    )
)

echo 完了しました。
endlocal