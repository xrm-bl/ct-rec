# ct-rec

CT reconstruction and processing suite for the SPring-8 synchrotron radiation
CT systems (parallel-beam absorption CT). Based on Nakano's software.

Raw projections (Hamamatsu HiPic `.img` or 16-bit TIFF) are normalized with
dark and I0 frames, ring artifacts are removed by a sorting-based filter
(Vo et al., Opt. Express 26, 2018), and slices are reconstructed by filtered
back projection into 32-bit float TIFF, then optionally normalized to
8/16-bit. Most reconstruction programs come in GPU (CUDA) and CPU (OpenMP)
variants and are configured through environment variables.

```
his/tif split -> axis search -> reconstruction (ring removal + FBP) -> 32bit tiff -> 8/16bit
```

## Main programs

| Category | Programs | Notes |
|---|---|---|
| Single-slice reconstruction (180°) | `ct_rec_G_F` | rotation-axis RMSD search built in |
| Whole-volume reconstruction | `hp_tg_G_F`, `ct_rec_loop` (`.bat`/`.sh`) | memory-adaptive chunking |
| Reconstruction from p-images | `p_rec_G_F`, `sf_rec_t_F` | |
| Offset CT (360°, half-beam) | `ofct_DO(_g)`, `ofct_rec_G_F`, `ofct_srec_G_F` | axis finder + reconstruction |
| Re-projection / re-reconstruction | `rec2rec_F`, `rec2rec_g_F` | |
| Sinogram / projection tools | `sinog`, `of_sinog`, `ct_prj_f`, `ict_prj_fc` | |
| 8/16-bit normalization | `tif_f2i` | min/max recorded in TIFF tags |
| 2D/3D filters for 32-bit volumes | `tif_mgf(_g)`, `tif_mdf(_g)`, `tif_gsf(_g)`, `tif_blf(_g)`, `tif_adf_g`, `tif_nlm_g`, `tif_tvd_g`, `tif_wvd_g`, `tif_bm4d_g`, `rec_gf` | median, gaussian, bilateral, anisotropic diffusion, NLM, TV, wavelet, BM4D |
| Acquisition-data utilities | `his_spl_*`, `act_spl(2)`, `his2img`, `his2tif6`, `*_ave` | split / convert / average |
| Volume utilities | `rec_crop`, `rec_sir` (3D binning), `rec_stk`, `tif2hst` | |
| ImageJ plugins | `imagej-plugins/SP8CT_Plugins.jar` | HIS/IMG opener, 2D/3D filters, stack crop, cylinder unwrap |

Suffix convention: `G` = `g` (GPU/CUDA) or `t` (CPU/OpenMP); `F` = FBP filter,
`r` (Ramachandran+Hann), `s` (Shepp-Logan), `c` (Chesler).
Example: `ct_rec_g_c` = GPU, Chesler.

## Robustness features

- Ring removal on GPU or CPU in every reconstruction program.
- Truncation (cupping) correction in the shared CBP layer when the sample
  overfills the field of view (`PAD_THRESH`, auto by default).
- Per-pixel low-transmission guard: opaque pixels are floored instead of
  producing Inf/NaN (`CT_REC_BLACK_THRESH` / `CT_REC_BLACK_FRAC`); clipping
  statistics are reported after each run. See
  [20260806_low_transmission_guard.md](20260806_low_transmission_guard.md).

## Binaries and building

Prebuilt Windows x64 binaries are in [exe/](exe/). To build from source:

```
cd src
makefileCPU.bat        (MSVC, OpenMP)
makefileGPU.bat        (MSVC + CUDA Toolkit 13.2, requires sm_75/Turing or newer)
```

libtiff 4.6.0 is bundled. Linux builds of selected tools: [bin/gen-all.sh](bin/gen-all.sh).

## Documentation

| File | Content |
|---|---|
| [readme.txt](readme.txt) | full manual, Japanese (current: ver 2.4) |
| [readme_e.txt](readme_e.txt) | full manual, English |
| [20260723_CT_env_vars.md](20260723_CT_env_vars.md) | all environment variables with defaults |
| [20260504_filter_readme.md](20260504_filter_readme.md) | filter tools |
| [20260806_low_transmission_guard.md](20260806_low_transmission_guard.md) | technical note: low-transmission guard |
| [imagej-plugins/readme.txt](imagej-plugins/readme.txt) | ImageJ plugin build/install |

Bug reports and requests: contact the author (Uesugi).
