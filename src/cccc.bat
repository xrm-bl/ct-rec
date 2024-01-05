

set   ICX=icl /DWINDOWS /O3 /Qparallel /Qprec-div
set   CC2=cl /DWINDOWS /O2
set   CCX=cl /DWINDOWS /Ox
set   CBP=cbp_thread_int.c
set SIF_F=sif_f_fast.c
set   CUDART="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2\lib\x64\cudart.lib"
set    CUFFT="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2\lib\x64\cufft.lib"
set CUDAINCL="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.2\include"


nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Chesler
%CC2% /Fect_rec_tif_g_c.exe  ct_rec_tif.c error.c libtiff.lib %CUFFT% %CUDART% cbp.obj

