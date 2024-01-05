
#ifdef	CSI_BS
typedef	struct {
		char	*path;
		double	base,step;
	} PBS;

extern PBS	*CheckSliceImages(
#else
extern char	**CheckSliceImages(
#endif
/*
		char	*directory,
		char	*nameFile,
		int	*Nx,
		int	*Ny,
		int	*Nz,
		int	*BPS
*/
		);
