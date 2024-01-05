
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cell.h"
#include "rif.h"
#include "sif.h"

extern void	Error(char *msg);

int	main(int argc,char **argv)
{
	int	Nx,Ny,BPS,Ox,Oy,y,x,y0,y1,x0,x1,v,u;
	Cell	**C;
	char	*desc;
	double	R,**F,**G,*f,R2,y2,*g,h,H=0.0;

	if (argc!=4 && argc!=5)
	    Error("usage : gf_sd orgTIFF radius {bias} newTIFF");

	ReadImageFile(argv[1],&Nx,&Ny,&BPS,&C,&desc);

	if ((R=atof(argv[2]))<=0.0) Error("bad Gaussian radius.");

	if (argc==5) H=atof(argv[3]);

	Ox=Oy=(int)ceil(sqrt(M_LN2*(double)(BPS+1))*R); if (Ox>Nx) Ox=Nx;
							if (Oy>Ny) Oy=Ny;

	if ((F=(double **)malloc(sizeof(double)*Ny))==NULL ||
	    (*F=(double *)malloc(sizeof(double)*Ny*Nx))==NULL ||
	    (G=(double **)malloc(sizeof(double *)*Oy))==NULL ||
	    (*G=(double *)malloc(sizeof(double)*Oy*Ox))==NULL)
	    Error("memory allocation error.");

	for (y=0; y<Ny; y++) {
	    if (y>0) F[y]=F[y-1]+Nx;

	    f=F[y]; for (x=0; x<Nx; x++) *(f++)=(double)C[y][x]-H;
	}
	R2=R*R;
	for (y=0; y<Oy; y++) {
	    if (y>0) G[y]=G[y-1]+Oy;

	    y2=(double)y*(double)y; g=G[y];
	    for (x=0; x<Ox; x++)
		*(g++)=exp(-((double)x*(double)x+y2)/R2)/(M_PI*R2);
	}
	for (y=0; y<Ny; y++) {
	    if ((y0=y-Oy+1)<0) y0=0;
	    if ((y1=y+Oy)>Ny) y1=Ny;

	    for (x=0; x<Nx; x++) {
		if ((x0=x-Ox+1)<0) x0=0;
		if ((x1=x+Ox)>Nx) x1=Nx; 

		h=H;
		for (v=y0; v<y; v++) {
		    f=F[v]+x0; g=G[y-v]+x-x0;
		    for (u=x0; u<x; u++) h+=(*(f++)*(*(g--)));
		    for (; u<x1; u++) h+=(*(f++)*(*(g++)));
		}
		for (; v<y1; v++) {
		    f=F[v]+x0; g=G[v-y]+x-x0;
		    for (u=x0; u<x; u++) h+=(*(f++)*(*(g--)));
		    for (; u<x1; u++) h+=(*(f++)*(*(g++)));
		}
		C[y][x]=(Cell)(h+0.5);
	    }
	}
	StoreImageFile(argv[argc-1],Nx,Ny,BPS,C,desc);

	free(*G); free(G); free(*F); free(F);
	for (y=0; y<Ny; y++) free(C[y]); free(C);
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
