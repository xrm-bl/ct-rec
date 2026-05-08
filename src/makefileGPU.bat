set   ICX=icl /DWINDOWS /O3 /Qparallel /Qprec-div
set   CC2=cl /DWINDOWS /O2 /D_USE_MATH_DEFINES
set   CCX=cl /DWINDOWS /Ox /D_USE_MATH_DEFINES
set   CBP=cbp_thread_int.c
set   SIF_F=sif_f_fast.c
set   CUDAVER=v13.2
set   CUDART="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\%CUDAVER%\lib\x64\cudart.lib"
set    CUFFT="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\%CUDAVER%\lib\x64\cufft.lib"
set CUDAINCL="C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\%CUDAVER%\include"
rem CUDA 13.x minimum: Turing (sm_75).  Change to sm_89 (Ada) or sm_100 (Blackwell) as needed.
set CUDA_ARCH=-arch=sm_75
set NVCC=nvcc -O3 -I%CUDAINCL% %CUDA_ARCH% -use_fast_math -Xcompiler "/wd 4819"


rem normal CT reconstruction
rem GPU
%NVCC% cbp.cu -DFloat=float -c -DFilter=Ramachandran
%CC2% /openmp /Fect_rec_g_r.exe ct_rec.c error.c sort_filter_omp.c %SIF_F% %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Ramachandran
%CC2% /openmp /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fehp_tg_g_r.exe hp_tg_ku.c error.c rhp.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Ramachandran
%CC2% /openmp /Fetf_rec_g_r.exe tf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Ramachandran
%CC2% /openmp /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fetf_tg_g_r.exe hp_tg_ku.c error.c rtf.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


%NVCC% cbp.cu -DFloat=float -c -DFilter=Shepp
%CC2% /openmp /Fect_rec_g_s.exe ct_rec.c error.c sort_filter_omp.c %SIF_F% %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Shepp
%CC2% /openmp /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fehp_tg_g_s.exe hp_tg_ku.c error.c rhp.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Shepp
%CC2% /openmp /Fetf_rec_g_s.exe  tf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Shepp
%CC2% /openmp /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fetf_tg_g_s.exe hp_tg_ku.c error.c rtf.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


%NVCC% cbp.cu -DFloat=float -c -DFilter=Chesler
%CC2% /openmp /Fect_rec_g_c.exe  ct_rec.c error.c sort_filter_omp.c %SIF_F% %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Chesler
%CC2% /openmp /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fehp_tg_g_c.exe hp_tg_ku.c error.c rhp.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Chesler
%CC2% /openmp /Fetf_rec_g_c.exe  tf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Chesler
%CC2% /openmp /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fetf_tg_g_c.exe hp_tg_ku.c error.c rtf.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


rem reconstruction for offset CT
rem GPU

%NVCC% cbp.cu -DFloat=float -c -DFilter=Ramachandran
%CC2% /openmp /Feofct_srec_g_r.exe ofct_srec.c error.c rhp.c rl.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Shepp
%CC2% /openmp /Feofct_srec_g_s.exe ofct_srec.c error.c rhp.c rl.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Chesler
%CC2% /openmp /Feofct_srec_g_c.exe ofct_srec.c error.c rhp.c rl.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


%NVCC% cbp.cu -DFloat=float -c -DFilter=Ramachandran
%CC2% /openmp /Feotf_rec_g_r.exe otf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Shepp
%CC2% /openmp /Feotf_rec_g_s.exe  otf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Chesler
%CC2% /openmp /Feotf_rec_g_c.exe  otf_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


%NVCC% cbp.cu -DFloat=float -c -DFilter=Ramachandran
%CC2% /openmp /Feoftf_srec_g_r.exe ofct_srec.c error.c rtf.c rl.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Shepp
%CC2% /openmp /Feoftf_srec_g_s.exe ofct_srec.c error.c rtf.c rl.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Chesler
%CC2% /openmp /Feoftf_srec_g_c.exe ofct_srec.c error.c rtf.c rl.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj


rem p image CT reconstruction
rem GPU

%NVCC% cbp.cu -DFloat=float -c -DFilter=Ramachandran
%CC2% /openmp /Fep_rec_g_r.exe /DFOM=float /DFloat=float p_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Shepp
%CC2% /openmp /Fep_rec_g_s.exe /DFOM=float /DFloat=float p_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

%NVCC% cbp.cu -DFloat=float -c -DFilter=Chesler
%CC2% /openmp /Fep_rec_g_c.exe /DFOM=float /DFloat=float p_rec.c error.c sort_filter_omp.c libtiff.lib %CUFFT% %CUDART% cbp.obj

rem filters made by ClaudAI
%NVCC% -o tif_blf_g.exe tif_blf_g.cu libtiff.lib
%NVCC% -o tif_gsf_g.exe tif_gsf_g.cu libtiff.lib
%NVCC% -o tif_mdf_g.exe tif_mdf_g.cu libtiff.lib
%NVCC% -o tif_nlm_g.exe tif_nlm_g.cu libtiff.lib
%NVCC% -o tif_tvd_g.exe tif_tvd_g.cu libtiff.lib
%NVCC% -o tif_wvd_g.exe tif_wvd_g.cu libtiff.lib
%NVCC% -o tif_adf_g.exe tif_adf_g.cu libtiff.lib
%NVCC% -o tif_bm4d_g.exe tif_bm4d_g.cu libtiff.lib

move *.exe ..\exe
del *.obj tif_*_g.exp  tif_*_g.lib
