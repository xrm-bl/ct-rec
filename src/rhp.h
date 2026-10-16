#ifndef	FOM
#define FOM	double
#endif

/* Two copies of the reader with different FOM in one program (his_srec:
   hp_tg_ku.c wants FOM=float, ofct_srec.c FOM=double).  A translation
   unit compiled with -DRHP_NAME_SUFFIX=_f sees and defines
   InitReadHiPic_f etc., so both variants can be linked side by side. */
#ifdef	RHP_NAME_SUFFIX
#define RHP_CAT_(a,b)	a##b
#define RHP_CAT(a,b)	RHP_CAT_(a,b)
#define InitReadHiPic	RHP_CAT(InitReadHiPic,RHP_NAME_SUFFIX)
#define ReadHiPic	RHP_CAT(ReadHiPic,RHP_NAME_SUFFIX)
#define ReadHiPicBand	RHP_CAT(ReadHiPicBand,RHP_NAME_SUFFIX)
#define TermReadHiPic	RHP_CAT(TermReadHiPic,RHP_NAME_SUFFIX)
#endif

typedef struct {
		int	q;
		double	c,
			a;
		int	it;
	} OutputLog;

typedef struct {
		char		*dir,
				**q_img;
		int		Nq,
				Nx,
				Ny,
				Ni,
				Nt;
		OutputLog	*OL;
		unsigned short	*D,
				**I;
		FOM		**T;
	} HiPic;

extern void	InitReadHiPic(	char	*dir,
				HiPic	*hipic),
		ReadHiPic(	HiPic	*hipic,
				int	t),
		/* 再入(スレッド安全)版: 投影 t を raw (Ny*Nx) に読み、補正した
		   透過率を行 y1..y2 について dst[y-y1][x] へ書く。異なる t に
		   対して並列に呼べる (rhp_c.c 参照)。 */
		ReadHiPicBand(	HiPic	*hipic,
				int	t,
				int	y1,
				int	y2,
				FOM	**dst,
				unsigned short	*raw),
		TermReadHiPic(	HiPic	*hipic);

#ifndef	ERROR_VALUE
#define ERROR_VALUE	0.0
#endif
