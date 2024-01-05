
#ifndef	FOM
#define FOM	double
#endif

extern int	ReadImageFile_Float(char *path,
				    int  *Nx,
				    int  *Ny,
				    FOM  ***cell,
				    char **desc);
