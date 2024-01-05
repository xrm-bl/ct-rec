
set   ICI=icl /DWINDOWS /O3 /D_USE_MATH_DEFINES /Qparallel /Qprec-div /Qstd=c99 /Dpopen=_popen /Dpclose=_pclose /D_WIN64 /D_AMD64_ 
set   ICX=icl /DWINDOWS /O3  /Qparallel /Qprec-div
set   CC2=cl /DWINDOWS /O2 /D_USE_MATH_DEFINES
set   CCX=cl /DWINDOWS /Ox /D_USE_MATH_DEFINES
set   CBP=cbp_thread_int.c
set   CBPd=cbp_thread_nai.c
set SIF_F=sif_f_fast.c


rem normal CT reconstruction
rem CPU

%CC2% /Fetf_rec_t_r.exe tf_rec.c error.c libtiff.lib %CBP% /DFilter=Ramachandran
%CC2% /Fetf_rec_t_s.exe tf_rec.c error.c libtiff.lib %CBP% /DFilter=Shepp
%CC2% /Fetf_rec_t_c.exe tf_rec.c error.c libtiff.lib %CBP% /DFilter=Chesler

%CC2% /Fect_rec_t_r.exe ct_rec.c error.c %SIF_F% %CBP% /DFilter=Ramachandran
%CC2% /Fect_rec_t_s.exe ct_rec.c error.c %SIF_F% %CBP% /DFilter=Shepp
%CC2% /Fect_rec_t_c.exe ct_rec.c error.c %SIF_F% %CBP% /DFilter=Chesler

%CC2% /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fehp_tg_t_r.exe hp_tg_ku.c error.c rhp.c libtiff.lib %CBP% /DFilter=Ramachandran
%CC2% /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fehp_tg_t_s.exe hp_tg_ku.c error.c rhp.c libtiff.lib %CBP% /DFilter=Shepp
%CC2% /DONLY_CT_VIEWS /DFOM=float /DFloat=float /Fehp_tg_t_c.exe hp_tg_ku.c error.c rhp.c libtiff.lib %CBP% /DFilter=Chesler

rem reconstruction for offset CT
rem CPU

%CC2% /Feofct_srec_t_r.exe /DFloat=float /DONLY_CT_VIEWS error.c rhp.c rl.c %CBP% libtiff.lib ofct_srec.c /DFilter=Ramachandran
%CC2% /Feofct_srec_t_s.exe /DFloat=float /DONLY_CT_VIEWS error.c rhp.c rl.c %CBP% libtiff.lib ofct_srec.c /DFilter=Shepp
%CC2% /Feofct_srec_t_c.exe /DFloat=float /DONLY_CT_VIEWS error.c rhp.c rl.c %CBP% libtiff.lib ofct_srec.c /DFilter=Chesler


rem guess rotation center for offset CT
%CC2% /Feofct_xy.exe /DONLY_CT_VIEWS error.c rhp.c msd.c %SIF_F% oct_xy.c

rem p image CT reconstruction
rem CPU

%CC2% /Fep_rec_t_r.exe p_rec.c error.c libtiff.lib %CBP% /DFloat=float /DFilter=Ramachandran
%CC2% /Fep_rec_t_s.exe p_rec.c error.c libtiff.lib %CBP% /DFloat=float /DFilter=Shepp
%CC2% /Fep_rec_t_c.exe p_rec.c error.c libtiff.lib %CBP% /DFloat=float /DFilter=Chesler

rem CT reconstruction from tiff 32bit sinogram
rem CPU

%CC2% /Fesf_rec_t_r.exe sf_rec.c error.c libtiff.lib %CBPd% /DFloat=double /DFilter=Ramachandran
%CC2% /Fesf_rec_t_s.exe sf_rec.c error.c libtiff.lib %CBPd% /DFloat=double /DFilter=Shepp
%CC2% /Fesf_rec_t_c.exe sf_rec.c error.c libtiff.lib %CBPd% /DFloat=double /DFilter=Chesler

rem normalize
%CC2% /Fetif_f2i.exe tif_f2i.c rif_f.c libtiff.lib

rem rec crop
%CC2% /Ferec_crop.exe rec_crop.c rif_f.c libtiff.lib

rem sinogram
%CC2% /Fesinog.exe sinog.c %SIF_F%

rem print image discription
%CC2% /Fepid.exe pid.c

rem stop watch
%CC2% stop_watch.c

rem img average
%CC2% /Feimg_ave.exe img_ave.c

rem spl
%CC2% /Fespl.exe spl.c
%CC2% /Fehis_spl_K.exe his_spl_K.c
%CC2% /Fehis_spl_E.exe his_spl_E.c
%CC2% /Fehis_spl_tif.exe his_spl_tif.c libtiff.lib

rem rec_stk
%CC2% /Ferec_stk.exe rec_stk.c
%CC2% /Fechk-rc.exe chk-rc.c
%CC2% /Feset-rc.exe set-rc.c

rem his2img
%CC2% /Fehis2img.exe his2img.c

rem ct_prj_f
%CC2% /Fect_prj_f.exe ct_prj_f.c libtiff.lib

rem ct_sub_f
%CC2% /Fect_sub_f.exe ct_sub_f.c libtiff.lib

rem ict_prj_fc
%CC2% /Feict_prj_fc.exe ict_prj_fc.c libtiff.lib

rem ct_prj_fc
%CC2% /Fect_prj_fc.exe ct_prj_fc.c libtiff.lib

rem tf_prj_fc
%CC2% /Fetf_prj_f.exe rif.c tf_prj_f.c libtiff.lib

rem tif2hst
%CC2% /Fetif2hst.exe tif2hst.c libtiff.lib

rem 3D gaussian filter
%CCX% /Ferec_gf.exe fft.c rec_gf.c rif_f.c libtiff.lib

rem 3D Rectangle rotation
rem si_rar.exe ro - +y +z +x ro_yz
rem si_rar.exe ro - +z +x +y ro_zx
%CC2% error.c rif.c csi.c rsi.c sif.c si_rar.c /Fesi_rar.exe

rem 3D binning
rem si_sir ro - 2 ro_2x2x2
%CC2% error.c rif.c csi.c rsi.c sif.c si_sir.c /Fesi_sir.exe
%CC2% error.c rif.c csi.c rsi.c si_sir_a.c libtiff.lib /Fesi_sir_a.exe
%CC2% error.c rif.c csi.c rsi.c si_sir_s.c libtiff.lib /Fesi_sir_s.exe

rem 3D gaussian filter for 8 or 16bit tiff
rem   si_gf  orgDir  nameFile  radius  {bias}  newDir
%CC2% error.c fft.c csi.c rif.c sif.c si_gf.c /Fesi_gf.exe

rem hp2DO
rem %ICX% hp2DO.c sif.c "C:\Program Files (x86)\IntelSWTools\compilers_and_libraries_2019.5.281\windows\compiler\lib\intel64_win\libiomp5md.lib"
%CCX% /openmp hp2DO.c sif.c 

rem tif_ave
%CC2% /Fetif_ave.exe tif_ave.c rif.c libtiff.lib

rem his_ave
%CC2% /Fehis_ave.exe his_ave.c sif.c

rem q3_ave, q4_ave
%CC2% /Feq3_ave.exe q3_ave.c
%CC2% /Feq4_ave.exe q4_ave.c

rem act_spl
rem %CC2% /Feact_spl.exe act_spl.c

rem gf_sd fd
%CC2% /Fegf_sd.exe error.c rif.c sif.c gf_sd.c
%CC2% /Fegf_fd.exe error.c fft.c rif.c sif.c gf_fd.c

rem his2tif6
%CC2% /Fehis2tif6.exe his2tif6.c libtiff.lib

move *.exe ..\exe
del *.obj

