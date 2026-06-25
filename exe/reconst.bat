@echo off
rem %1=layer		// not negligible
rem %2=center		// È—ª‚·‚é‚Æˆê‰dSŒvZ‚ğ‚·‚é

if "%1" =="" goto err

@echo on
ct_rec_g_c %1 %2 %3 %4
rem ct_sino %1 
rem sino_conv %2
rem ct_cbp %4
rem rec2tif 16 %3
@echo off
rem rec2img
goto fin

:err
echo Usage is: reconst layer (center) (pixel size) (offset angle)

:fin
