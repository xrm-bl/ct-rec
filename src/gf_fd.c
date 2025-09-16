
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"
#include "cell.h"
#include "rif.h"
#include "sif.h"

extern void	Error(char *msg);

int	main(int argc,char **argv)
{
	int	Nx,Ny,BPS,Lx,Ly,y,x;
	Cell	**cell,c;
	char	*desc;
	double	R,d,PI2R2,x2,g,bias=0.0;
	Complex	*Ex,*Ey,**C,*F;

	if (argc!=4 && argc!=5)
	    Error("usage : gf_fd orgTIFF radius {bias} newTIFF");

	ReadImageFile(argv[1],&Nx,&Ny,&BPS,&cell,&desc);

	if ((R=atof(argv[2]))<0.0) Error("bad Gaussian radius.");

	if (argc==5) bias=atof(argv[3]);

	if ((Lx=SetUpFFT(Nx,&Ex))==0 ||
	    (Ly=SetUpFFT(Ny,&Ey))==0 ||
	    (C=(Complex **)malloc(sizeof(Complex *)*Ly))==NULL ||
	    (*C=(Complex *)malloc(sizeof(Complex)*Ly*Lx))==NULL ||
	    (F=(Complex *)malloc(sizeof(Complex)*Ly))==NULL)
	    Error("memory allocation error.");

	d=1.0/((double)Lx*(double)Ly);
	for (y=0; y<Ny; y++) {
	    if (y>0) C[y]=C[y-1]+Lx;

	    for (x=0; x<Nx; x++) {
		C[y][x].r=((double)cell[y][x]-bias)*d; C[y][x].i=0.0;
	    }
	    for (; x<Lx; x++) C[y][x].r=C[y][x].i=0.0;

	    FFT(-1,Lx,Ex,C[y]);
	}
	for (; y<Ly; y++) {
	    C[y]=C[y-1]+Lx; for (x=0; x<Lx; x++) C[y][x].r=C[y][x].i=0.0;
	}
	PI2R2=M_PI*M_PI*R*R;
	for (x=0; x<Lx; x++) {
	    for (y=0; y<Ly; y++) F[y]=C[y][x]; FFT(-1,Ly,Ey,F);

	    d=0.5-fabs(0.5-(double)x/(double)Lx); x2=d*d;
	    for (y=0; y<Ly; y++) {
		d=0.5-fabs(0.5-(double)y/(double)Ly);
		g=exp(-PI2R2*(x2+d*d)); F[y].r*=g; F[y].i*=g;
	    }
	    FFT( 1,Ly,Ey,F); for (y=0; y<Ly; y++) C[y][x]=F[y];
        }
	c=(Cell)((1<<BPS)-1);
	for (y=0; y<Ny; y++) {
	    FFT( 1,Lx,Ex,C[y]);

	    for (x=0; x<Nx; x++)
		cell[y][x]=((d=C[y][x].r+bias)<0.0)?0:
			   (d>(double)c)?c:
			   (Cell)(d+0.5);
	}
	StoreImageFile(argv[argc-1],Nx,Ny,BPS,cell,desc);

	free(F); free(*C); free(C); free(Ey); free(Ex);
	for (y=0; y<Ny; y++) free(cell[y]); free(cell);
	if (desc!=NULL) free(desc);

		// append to log file
	FILE		*flog;
	int		i;
	if ((flog = fopen("cmd-hst.log", "a")) == NULL) {
		return(-1);
	}
	for (i = 0; i<argc; ++i) fprintf(flog, "%s ", argv[i]);
	fprintf(flog, "\n");
	fclose(flog);

	return 0;
}
