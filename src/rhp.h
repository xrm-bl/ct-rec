
#ifndef	FOM
#define FOM	double
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
