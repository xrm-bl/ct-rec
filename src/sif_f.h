
#ifndef	FOM
#define FOM	double
#endif

extern void	StoreImageFile_Float(char *path,
				     int  Nx,
				     int  Ny,
				     FOM  **cell,
				     char *desc);

extern FOM	SIF_F_min,SIF_F_max;

#ifndef	SIF_F_LEN
#define SIF_F_LEN	2048
#endif

#ifndef	SIF_F_FORM
#define SIF_F_FORM	"%lf\t%lf"
#endif

extern char	SIF_F_form[],SIF_F_desc[];
