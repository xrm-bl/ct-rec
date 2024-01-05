
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"
#include "csi.h"
#include "cell.h"
#include "rif.h"
#include "sif.h"

#ifndef	COMPLEX
typedef	struct {
		float	r,i;
	} COMPLEX;
#endif

extern void	Error(char *msg);

#define LEN	1024

int	main(int argc,char **argv)
{
	char	**path,**desc,str[LEN];
	int	Nx,Ny,Nz,BPS,Lx,Ly,Lz,z,nx,ny,y,x,s;
	double	R,d,PI2R2,y2,x2_y2,g,bias=0.0;
	Complex	*Ex,*Ey,*Ez,*F;
	COMPLEX	***C;
	Cell	**cell,c;

	if (argc!=5 && argc!=6)
	    Error("usage : si_gf orgDir nameFile radius {bias} newDir");

	path=CheckSliceImages(argv[1],argv[2],&Nx,&Ny,&Nz,&BPS);

	if ((R=atof(argv[3]))<0.0) Error("bad Gaussian radius.");

	if (argc==6) bias=atof(argv[4]);

	if ((desc=(char **)malloc(sizeof(char *)*Nz))==NULL ||
	    (Lx=SetUpFFT(Nx,&Ex))==0 ||
	    (Ly=SetUpFFT(Ny,&Ey))==0 ||
	    (Lz=SetUpFFT(Nz,&Ez))==0 ||
	(F=(Complex *)malloc(sizeof(Complex)*((Lx>Ly)?(Lx>Lz)?Lx:Lz:
						      (Ly>Lz)?Ly:Lz)))==NULL ||
	    (C=(COMPLEX ***)malloc(sizeof(COMPLEX **)*Nz))==NULL ||
	    (*C=(COMPLEX **)malloc(sizeof(COMPLEX *)*Nz*Ly))==NULL ||
	    (**C=(COMPLEX *)malloc(sizeof(COMPLEX)*Nz*Ly*Lx))==NULL)
	    Error("memory allocation error.");

	d=1.0/((double)Lx*(double)Ly*(double)Lz);
	for (z=0; z<Nz; z++) {
	    if (z>0) C[z]=C[z-1]+Ly;

	    ReadImageFile(path[z],&nx,&ny,NULL,&cell,desc+z);

	    for (y=0; y<ny; y++) {
		if (z>0 || y>0) C[z][y]=C[z][y-1]+Lx;

		for (x=0; x<nx; x++) {
		    F[x].r=((double)cell[y][x]-bias)*d; F[x].i=0.0;
		}
		for (; x<Lx; x++) F[x].r=F[x].i=0.0;

		FFT(-1,Lx,Ex,F);
		for (x=0; x<Lx; x++) {
		    C[z][y][x].r=F[x].r; C[z][y][x].i=F[x].i;
		}
		free(cell[y]);
	    }
	    free(cell);

	    for (; y<Ly; y++) C[z][y]=C[z][y-1]+Lx;

	    for (x=0; x<Lx; x++) {
		for (y=0; y<ny; y++) {
		    F[y].r=C[z][y][x].r; F[y].i=C[z][y][x].i;
		}
		for (; y<Ly; y++) F[y].r=F[y].i=0.0;

		FFT(-1,Ly,Ey,F);
		for (y=0; y<Ly; y++) {
		    C[z][y][x].r=F[y].r; C[z][y][x].i=F[y].i;
		}
	    }
	}
	PI2R2=M_PI*M_PI*R*R;
	for (y=0; y<Ly; y++) {
	    d=0.5-fabs(0.5-(double)y/(double)Ly); y2=d*d;
	    for (x=0; x<Lx; x++) {
		for (z=0; z<Nz; z++) {
		    F[z].r=C[z][y][x].r; F[z].i=C[z][y][x].i;
		}
		for (; z<Lz; z++) F[z].r=F[z].i=0.0;

		FFT(-1,Lz,Ez,F);

		d=0.5-fabs(0.5-(double)x/(double)Lx); x2_y2=d*d+y2;
		for (z=0; z<Lz; z++) {
		    d=0.5-fabs(0.5-(double)z/(double)Lz);
		    g=exp(-PI2R2*(x2_y2+d*d)); F[z].r*=g; F[z].i*=g;
		}
		FFT(1,Lz,Ez,F);
		for (z=0; z<Nz; z++) {
		    C[z][y][x].r=F[z].r; C[z][y][x].i=F[z].i;
		}
	    }
        }
	cell=(Cell **)*C;
	c=(Cell)((1<<BPS)-1);
	s=0; while (argv[1][s++]!='\0') ;

	for (z=0; z<Nz; z++) {
	    for (x=0; x<Lx; x++) {
		for (y=0; y<Ly; y++) {
		    F[y].r=C[z][y][x].r; F[y].i=C[z][y][x].i;
		}
		FFT(1,Ly,Ey,F);
		for (y=0; y<Ly; y++) {
		    C[z][y][x].r=F[y].r; C[z][y][x].i=F[y].i;
		}
	    }
	    for (y=0; y<Ny; y++) {
		for (x=0; x<Lx; x++) {
		    F[x].r=C[z][y][x].r; F[x].i=C[z][y][x].i;
		}
		FFT(1,Lx,Ex,F);
		for (x=0; x<Nx; x++)
		    cell[y][x]=((d=F[x].r+bias)<0.0)?0:
			       (d>(double)c)?c:
			       (Cell)(d+0.5);
	    }
	    (void)sprintf(str,"%s/%s",argv[argc-1],path[z]+s);
	    StoreImageFile(str,Nx,Ny,BPS,cell,desc[z]);

	    if (desc[z]!=NULL) free(desc[z]);

	    free(path[z]);
	}
	free(**C); free(*C); free(C); free(F); free(Ez); free(Ey); free(Ex);
	free(desc); free(path);

	// append to log file
	FILE		*f;
	int			i;
	if((f = fopen("cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(f,"%s ",argv[i]);
	fprintf(f,"\n");
	fclose(f);

	return 0;
}
