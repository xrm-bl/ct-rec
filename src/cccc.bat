

set   ICX=icl /DWINDOWS /O3 /Qparallel /Qprec-div
set   CC2=cl /DWINDOWS /O2 /D_USE_MATH_DEFINES
set   CCX=cl /DWINDOWS /Ox /D_USE_MATH_DEFINES
set   CBP=cbp_thread_int.c
set SIF_F=sif_f_fast.c
set   CUDAVER=v13.2
set   CUDART="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\%CUDAVER%\lib\x64\cudart.lib"
set    CUFFT="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\%CUDAVER%\lib\x64\cufft.lib"
set CUDAINCL="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\%CUDAVER%\include"
rem CUDA 13.x minimum: Turing (sm_75).  Change to sm_89 (Ada) or sm_100 (Blackwell) as needed.
set CUDA_ARCH=-arch=all-major
set NVCC=nvcc -O3 -I%CUDAINCL% %CUDA_ARCH% -use_fast_math -Xcompiler "/wd 4819"


rem GPU ring removal object (pure CUDA C); built once, linked with /DUSE_GPU.
%NVCC% sort_filter_g.cu -c

%NVCC% cbp.cu -DFloat=float -c -DFilter=Chesler
%CC2% /openmp /DUSE_GPU /Fetf_rec_g_c.exe  tf_rec.c error.c sort_filter_g.obj libtiff.lib %CUFFT% %CUDART% cbp.obj

