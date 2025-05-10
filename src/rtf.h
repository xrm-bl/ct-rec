#ifndef FOM
#define FOM double
#endif

// TIFFライブラリに関する依存性をインクルード
#ifndef TIFF_SUPPORT
#define TIFF_SUPPORT
#endif

typedef struct {
    int    q;
    double c,
           a;
    int    it;
} OutputLog;
        
typedef struct {
    char        *dir,
                **q_img;
    int         Nq,
                Nx,
                Ny,
                Ni,
                Nt;
    OutputLog   *OL;
    unsigned short *D,
                **I;
    FOM         **T;
} HiPic;

extern void InitReadHiPic(char *dir,
                         HiPic *hipic),
           ReadHiPic(HiPic *hipic,
                    int t),
           TermReadHiPic(HiPic *hipic);

#ifndef ERROR_VALUE
#define ERROR_VALUE 0.0
#endif