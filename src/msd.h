
#ifndef	FOM
#define FOM	double
#endif

typedef	struct {
		double	r,i;
	} Complex;

typedef struct {
		int	Ngx,Ngy,Nhx,Nhy,Dx1,Dy1,Dx2,Dy2,Ogx,Ogy,Ohx,Ohy,Lx,Ly;
		Complex	*Ex,*Ey,**Q,*U,*V;
		FOM	*D,*C;
	} MSD;

extern int	InitMSD(MSD *msd,int Ngx,int Ngy,int Nhx,int Nhy,
				 int Dx1,int Dy1,int Dx2,int Dy2);
extern void	CalcMSD(MSD *msd,FOM **G,FOM **H),
		TermMSD(MSD *msd);
