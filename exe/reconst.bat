@echo off
rem %1=layer		// not negligible
rem %2=center		// 省略すると一応重心計算をする

if "%1" =="" goto err

@echo on
ct_rec_g_c %1 %2 %3 %4
@echo off
rem ct_sino %1 
rem sino_conv %2
rem ct_cbp %4
rem rec2tif 16 %3
rem rec2img
goto fin

:err
echo Usage is: reconst layer (center) (pixel size) (offset angle)

:fin
