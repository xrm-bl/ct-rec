#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>  // uint32_t, uint16_tのための定義を追加
#ifdef _WIN32
#include "msdirent.h"
#else
#include <dirent.h>
#endif
#include "tiffio.h"  // TIFFライブラリをインクルード
#include "rtf.h"

#define EPS 1e-9  /* for time interval */

static void Error(char *dir, char *name, char *msg)
{
    if (*dir != '\0') (void)fprintf(stderr, "%s/", dir);
    
    (void)fprintf(stderr, "%s : %s\n", name, msg); exit(1);
}

static int StrCmp(char *s1, char *s2)
{
    int c;
    
    while ((c = (*s1 >= 'A' && *s1 <= 'Z') ? *s1 | 32 : *s1) == (*s2)) {
        if (c == '\0') return 0;
        
        ++s1; ++s2;
    }
    return 1;
}

#define LEN 2048

static TIFF* OpenTiff(char *dir, char *name, int *Nx, int *Ny)
{
    char path[LEN];
    TIFF *tif;
    
    (void)sprintf(path, "%s/%s", dir, name);
    if ((tif = TIFFOpen(path, "r")) == NULL) Error(dir, name, "file not open.");
    
    // TIFF画像のサイズを取得
    uint32_t width, height;
    
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);
    
    *Nx = (int)width;
    *Ny = (int)height;
    
    // 16bitモノクロチェック
    uint16_t bps, spp, photo;
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
    TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &spp);
    TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photo);
    
    if (bps != 16 || spp != 1 || (photo != PHOTOMETRIC_MINISBLACK && photo != PHOTOMETRIC_MINISWHITE))
        Error(dir, name, "not a 16-bit monochrome TIFF image.");
    
    return tif;
}

#define WORD unsigned short

static void ReadTiff(char *dir, char *name, int Nx, int Ny, WORD *W)
{
    int x, y;
    TIFF *tif = OpenTiff(dir, name, &x, &y);
    char *dre = "data read error.";
    
    if (Nx != x || Ny != y) Error(dir, name, "image size not match.");
    
    // 画像データを読み込む
    tsize_t scanline_size = TIFFScanlineSize(tif);
    tdata_t buf = _TIFFmalloc(scanline_size);
    
    if (buf == NULL) Error(dir, name, "no memory for TIFF buffer.");
    
    for (y = 0; y < Ny; y++) {
        if (TIFFReadScanline(tif, buf, y, 0) < 0) Error(dir, name, dre);
        
        // 1行分のデータをコピー
        memcpy(W + y * Nx, buf, Nx * sizeof(WORD));
    }
    
    _TIFFfree(buf);
    TIFFClose(tif);
}

static int Compare(OutputLog *OL1, OutputLog *OL2)
{
    int it = OL2->it - OL1->it;
    
    return (it) ? it : (OL1->it == 0) ? (OL1->c > OL2->c) ? 1 : (OL1->c < OL2->c) ? -1 : 0 :
                       (OL1->it == 1) ? (OL1->a > OL2->a) ? 1 : (OL1->a < OL2->a) ? -1 : 0 : 0;
}

void InitReadHiPic(char *dir, HiPic *hp)
{
    char *env, q_img[LEN], str[LEN],
         output_log[LEN] = "",
         dark_img[LEN] = "",
         *fnf = "file not found.",
         *nmfdd = "no memory for directory data.",
         *bsn = "bad sequence number.",
         *dsn = "duplicated sequence number.";
    DIR *Dir;
    struct dirent *sd;
    int l, q, it, x, y, i;
    FILE *file;
    double c, a, c1, c2;
    OutputLog *OL;
    
    // 環境変数の読み込み
    if ((env = getenv("RHP_O")) != NULL) (void)strcpy(output_log, env);
    if ((env = getenv("RHP_D")) != NULL) (void)strcpy(dark_img, env);
    
    // TIFFファイルのパターンを設定
    (void)sprintf(q_img,
                  "%c[0-9]*[0-9].tif",
                  ((env = getenv("RHP_Q")) != NULL && *env != '\0' &&
                   (*env | 32) >= 'a' && (*env | 32) <= 'z') ? *env | 32 : 'q');
    
    hp->Nq = 0;
    
    if ((Dir = opendir(dir)) == NULL) Error("", dir, "directory not open.");
    
#define NAME sd->d_name
    
    while ((sd = readdir(Dir)) != NULL)
        if (*output_log == '\0' && !StrCmp(NAME, "output.log"))
            (void)sprintf(output_log, "%s/%s", dir, NAME);
        else if (*dark_img == '\0' && !StrCmp(NAME, "dark.tif"))
            (void)strcpy(dark_img, NAME);
        else if ((NAME[0] | 32) == q_img[0] &&
                 (l = (int)strlen(NAME)) > 5 &&
                 !StrCmp(NAME + l - 4, ".tif") &&
                 sscanf(NAME + 1, "%d", &q) == 1 &&
                 q >= hp->Nq)
            hp->Nq = q + 1;
    
    if (*output_log == '\0') Error(dir, "output.log", fnf);
    if (*dark_img == '\0') Error(dir, "dark.tif", fnf);
    if (hp->Nq == 0) Error(dir, q_img, fnf);
    
#define MALLOC(type, noe)    (type *)malloc(sizeof(type) * (noe))
    
    if ((hp->dir = strdup(dir)) == NULL ||
        (hp->q_img = MALLOC(char *, hp->Nq)) == NULL) Error("", dir, nmfdd);
    
    rewinddir(Dir);
    
    for (q = 0; q < hp->Nq; q++) hp->q_img[q] = NULL;
    
    while ((sd = readdir(Dir)) != NULL)
        if ((NAME[0] | 32) == q_img[0] &&
            (l = (int)strlen(NAME)) > 5 &&
            !StrCmp(NAME + l - 4, ".tif") &&
            sscanf(NAME + 1, "%d", &q) == 1 &&
            q >= 0) {
            if (q >= hp->Nq) Error(dir, NAME, bsn);
            
            if (hp->q_img[q] != NULL) Error(dir, NAME, dsn);
            
            if ((hp->q_img[q] = strdup(NAME)) == NULL) Error(dir, NAME, nmfdd);
        }
    
    (void)closedir(Dir);
    
    if ((hp->OL = MALLOC(OutputLog, hp->Nq + 1)) == NULL)
        Error("", output_log, "no memory for log data.");
    
    for (q = 0; q <= hp->Nq; q++) hp->OL[q].it = (-1);
    
    // 暗画像のサイズを取得
    TIFF *tif = OpenTiff(dir, dark_img, &(hp->Nx), &(hp->Ny));
    TIFFClose(tif);
    
    hp->Ni = hp->Nt = 0;
    
    if (*output_log == '-')
        file = stdin;
    else
        if ((file = fopen(output_log, "r")) == NULL)
            Error("", output_log, "file not open.");
    
    while (fgets(str, LEN, file) != NULL)
        if ((sscanf(str, "%d %lf %lf %d%n", &q, &c, &a, &it, &l) == 4)) {
            str[l] = '\0';
            
            if (q < 0 || q >= hp->Nq || hp->q_img[q] == NULL) Error("", str, bsn);
            
            if (it != 0 && it != 1) Error("", str, "unacceptable log data.");
            
            // 画像サイズチェック
            tif = OpenTiff(dir, hp->q_img[q], &x, &y);
            TIFFClose(tif);
            
            if (hp->Nx != x || hp->Ny != y)
                Error(dir, hp->q_img[q], "image size not match.");
            
            if (hp->OL[q].it != (-1)) Error("", str, dsn);
            
            if (it == 0) ++(hp->Ni); else ++(hp->Nt);
            
            if (hp->Ni + hp->Nt == 1) c1 = c2 = c; else if (c1 > c) c1 = c;
                                       else if (c2 < c) c2 = c;
            
            hp->OL[q].q = q; hp->OL[q].c = c; hp->OL[q].a = a; hp->OL[q].it = it;
        }
        else
#ifndef ONLY_CT_VIEWS
            if (file == stdin)
#endif
            break;
    
    if (file != stdin) (void)fclose(file);
    
    if (hp->Ni == 0) Error(dir, q_img, "no file referred as I0-image.");
    if (hp->Nt == 0) Error(dir, q_img, "no file referred as I-image.");
    
    if ((hp->D = MALLOC(WORD, hp->Ny * hp->Nx)) == NULL ||
        (hp->I = MALLOC(WORD *, hp->Ni + (hp->Ni == 1))) == NULL ||
        (*(hp->I) = MALLOC(WORD, hp->Ni * hp->Ny * hp->Nx)) == NULL ||
        (hp->T = MALLOC(FOM *, hp->Ny)) == NULL ||
        (*(hp->T) = MALLOC(FOM, hp->Ny * hp->Nx)) == NULL)
        Error("", dir, "no memory for image data.");
    
    ReadTiff(dir, dark_img, hp->Nx, hp->Ny, hp->D);
    
    qsort(hp->OL, hp->Nq + 1, sizeof(OutputLog), (int (*)())Compare);
    
    OL = hp->OL + hp->Nt;
    for (i = 0; i < hp->Ni; i++) {
        if (i > 0) {
            if (OL[i - 1].c + EPS >= OL[i].c) {
                fprintf(stderr, "%d\t%lf\t%lf\n", i, OL[i - 1].c + EPS, OL[i].c);
                // Error("", output_log, "bad time interval.");
            }
            
            hp->I[i] = hp->I[i - 1] + (size_t)hp->Ny * hp->Nx;
        }
        ReadTiff(dir, hp->q_img[OL[i].q], hp->Nx, hp->Ny, hp->I[i]);
    }
    if (hp->Ni == 1) {
        hp->Ni = 2; OL[1].q = OL[0].q;
                   OL[0].c = c1 - EPS; OL[1].c = c2 + EPS;
                   OL[1].a = OL[0].a; OL[1].it = OL[0].it; hp->I[1] = hp->I[0];
    }
    for (y = 1; y < hp->Ny; y++) hp->T[y] = hp->T[y - 1] + hp->Nx;
}

void ReadHiPic(HiPic *hp, int t)
{
    char str[LEN];
    OutputLog *ol, *OL;
    int i, y, x;
    WORD *D, *I1, *I2, *T;
    double r1, r2, I_D, T_D;
    
    if (t < 0 || t >= hp->Nt) {
        (void)sprintf(str, "%d", t); Error("", str, "bad sequence number.");
    }
    ol = hp->OL + t;
    ReadTiff(hp->dir, hp->q_img[ol->q], hp->Nx, hp->Ny, (WORD *)*(hp->T));
    
    OL = hp->OL + hp->Nt;
    for (i = hp->Ni - 2; i > 0; i--) if (ol->c >= OL[i].c) break;
    
    D = hp->D + (size_t)hp->Ny * hp->Nx;
    I1 = hp->I[i] + (size_t)hp->Ny * hp->Nx;
    I2 = hp->I[i + 1] + (size_t)hp->Ny * hp->Nx;
    T = (WORD *)*(hp->T) + (size_t)hp->Ny * hp->Nx;
    r1 = 1.0 - (r2 = (ol->c - OL[i].c) / (OL[i + 1].c - OL[i].c));
    
    for (y = hp->Ny - 1; y >= 0; y--)
    for (x = hp->Nx - 1; x >= 0; x--) {
        --D; --I1; --I2; --T;
        hp->T[y][x] = ((I_D = r1 * (double)*I1 + r2 * (double)*I2 - (double)*D) > 0.0 &&
                      (T_D = (double)*T - (double)*D) > 0.0) ? T_D / I_D : ERROR_VALUE;
    }
}

void TermReadHiPic(HiPic *hp)
{
    int q;
    
    free(*(hp->T)); free(hp->T);
    free(*(hp->I)); free(hp->I); free(hp->D);
    
    free(hp->OL);
    
    for (q = 0; q < hp->Nq; q++) if (hp->q_img[q] != NULL) free(hp->q_img[q]);
    
    free(hp->q_img); free(hp->dir);
}