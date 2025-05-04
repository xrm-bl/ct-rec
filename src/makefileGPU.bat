
set   ICX=icl /DWINDOWS /O3 /Qparallel /Qprec-div
set   CC2=cl /DWINDOWS /O2 /D_USE_MATH_DEFINES
set   CCX=cl /DWINDOWS /Ox /D_USE_MATH_DEFINES
set   CBP=cbp_thread_int.c
set   SIF_F=sif_f_fast.c
set   CUDART="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\lib\x64\cudart.lib"
set    CUFFT="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\lib\x64\cufft.lib"
set CUDAINCL="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\include"


rem normal CT reconstruction
rem GPU
nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Ramachandran
%CC2% /openmp /Fect_rec_g_r.exe ct_rec.c error.c sort_filter_omp.c %SIF_F% %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -DFloat=float -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Ramachandran
%CC2% /openmp /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fehp_tg_g_r.exe hp_tg_ku.c error.c rhp.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Ramachandran
%CC2% /openmp /Fetf_rec_g_r.exe tf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Shepp
%CC2% /openmp /Fect_rec_g_s.exe ct_rec.c error.c sort_filter_omp.c %SIF_F% %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -DFloat=float -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Shepp
%CC2% /openmp /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fehp_tg_g_s.exe hp_tg_ku.c error.c rhp.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Shepp
%CC2% /openmp /Fetf_rec_g_s.exe  tf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Chesler
%CC2% /openmp /Fect_rec_g_c.exe  ct_rec.c error.c sort_filter_omp.c %SIF_F% %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -DFloat=float -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Chesler
%CC2% /openmp /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fehp_tg_g_c.exe hp_tg_ku.c error.c rhp.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Chesler
%CC2% /openmp /Fetf_rec_g_c.exe  tf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

rem reconstruction for offset CT
rem GPU

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Ramachandran
%CC2% /openmp /Feofct_srec_g_r.exe ofct_srec.c error.c rhp.c rl.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Shepp
%CC2% /openmp /Feofct_srec_g_s.exe ofct_srec.c error.c rhp.c rl.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Chesler
%CC2% /openmp /Feofct_srec_g_c.exe ofct_srec.c error.c rhp.c rl.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Ramachandran
%CC2% /openmp /Feotf_rec_g_r.exe otf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Shepp
%CC2% /openmp /Feotf_rec_g_s.exe  otf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Chesler
%CC2% /openmp /Feotf_rec_g_c.exe  otf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


rem p image CT reconstruction
rem GPU

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Ramachandran
%CC2% /openmp /Fep_rec_g_r.exe /DFOM=float /DFloat=float p_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Shepp
%CC2% /openmp /Fep_rec_g_s.exe /DFOM=float /DFloat=float p_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

nvcc -O3 -I%CUDAINCL% cbp.cu -c -use_fast_math -Xcompiler "/wd 4819" -DFilter=Chesler
%CC2% /openmp /Fep_rec_g_c.exe /DFOM=float /DFloat=float p_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


move *.exe ..\exe
del *.obj

