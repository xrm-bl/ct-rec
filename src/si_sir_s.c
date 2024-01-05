
#include <stdio.h>
#include <stdlib.h>
#include "cell.h"
#include "csi.h"
//#include "sif.h"
#include "tiffio.h"

extern void	Error(),ReadSliceImage();

#define LEN	1024

void Store16TiffFile(char *wname, int wX, int wY, int wBPS, unsigned short *data16, char *wdesc)
{
	TIFF *image;
	long i;

	image = TIFFOpen(wname, "w");

	TIFFSetField(image, TIFFTAG_IMAGEWIDTH, wX);
	TIFFSetField(image, TIFFTAG_IMAGELENGTH, wY);
	TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, wBPS);
	TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, 1);
	TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, wdesc);
	TIFFSetField(image, TIFFTAG_ARTIST, "tif_f2i");

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data16 + i*wX, wX * sizeof(unsigned short));
	}

	TIFFClose(image);
}


int	main(argc,argv)
int	argc;
char	**argv;
{
	char	**paths,form[LEN],path[LEN],*bbs="bad block size.";
	int	Nx,Ny,Nz,BPS,Bx,By,Bz,y,H,Nh,V,Nv,v,h,d,z,x;
	Cell	**cell;
	double	**sum,xy,hy,xv,hv,pels;
	
	unsigned short	*data16;
	int		dBPS;
	char		wdesc[300];
	long		lll;

	int		i;

	if (argc!=5 && argc!=7)
	    Error("usage : si_sir orgDir nameFile Bxyz newDir\n"
		  "        si_sir orgDir nameFile Bx By Bz newDir");

	paths=CheckSliceImages(argv[1],argv[2],&Nx,&Ny,&Nz,&BPS);

	if (argc==5) {
	    if ((Bx=By=Bz=atoi(argv[3]))<=0) Error(bbs);
	}
	else
	    if ((Bx=atoi(argv[3]))<=0 ||
		(By=atoi(argv[4]))<=0 ||
		(Bz=atoi(argv[5]))<=0) Error(bbs);

#define ALLOC(type,noe)	(type *)malloc(sizeof(type)*(noe))

	if ((cell=ALLOC(Cell *,(size_t)Ny))==NULL ||
	    (*cell=ALLOC(Cell,(size_t)Ny*(size_t)Nx))==NULL)
	    Error("no allocatable memory for cell.");

	for (y=1; y<Ny; y++) cell[y]=cell[y-1]+Nx;

	Nh=(H=(Nx-1)/Bx)+1;
	Nv=(V=(Ny-1)/By)+1;

	if ((sum=ALLOC(double *,(size_t)Nv))==NULL ||
	    (*sum=ALLOC(double,(size_t)Nv*(size_t)Nh))==NULL)
	    Error("no allocatable memory for sum.");

	for (v=1; v<Nv; v++) sum[v]=sum[v-1]+Nh;

	sprintf(wdesc, "sum\t%d\t%d", Bx, By);
	data16 = (unsigned short*)malloc(sizeof(unsigned short)*Nh*Nv);
	dBPS=16;

	xy=(double)Bx*(double)By;
	h=(Nx-1)%Bx+1; hy=(double)h*(double)By;
	v=(Ny-1)%By+1; xv=(double)Bx*(double)v; hv=(double)h*(double)v;
{
	int	l=1;

	for (d=(Nz-1)/Bz; d>=10; d/=10) ++l;
	(void)sprintf(form,"%s/a%%0%dd.tif",argv[argc-1],l);
}
	for (z=0; z<Nz; z++) {
	    ReadSliceImage(paths[z],Nx,Ny,cell);

	    if (z%Bz==0)
		for (v=0; v<Nv; v++)
		for (h=0; h<Nh; h++) sum[v][h]=0.0;

	    for (y=0; y<Ny; y++) {
		v=y/By; for (x=0; x<Nx; x++) sum[v][x/Bx]+=(double)cell[y][x];
	    }
	    if ((z+1)%Bz==0 || z+1==Nz) {
		pels=xy*(double)(d=z%Bz+1);
		for (v=0; v<V; v++)
//		for (h=0; h<H; h++) cell[v][h]=(Cell)(sum[v][h]/pels+0.5);
		for (h=0; h<H; h++) cell[v][h]=(Cell)(sum[v][h]);

		pels=hy*(double)d;
//		for (v=0; v<V; v++) cell[v][H]=(Cell)(sum[v][H]/pels+0.5);
		for (v=0; v<V; v++) cell[v][H]=(Cell)(sum[v][H]);

		pels=xv*(double)d;
//		for (h=0; h<H; h++) cell[V][h]=(Cell)(sum[V][h]/pels+0.5);
		for (h=0; h<H; h++) cell[V][h]=(Cell)(sum[V][h]);

//		cell[V][H]=(Cell)(sum[V][H]/hv/(double)d+0.5);
		cell[V][H]=(Cell)(sum[V][H]);

		lll=0;
		for (v=0; v<Nv; v++){
			for (h=0; h<Nh; h++) {
				*(data16+lll)=sum[v][h];
				lll=lll+1;
			}
		}
		(void)sprintf(path,form,z/Bz+1);
//		StoreImageFile(path,Nh,Nv,BPS,cell,NULL);
		Store16TiffFile(path, Nh, Nv, dBPS, data16, wdesc);
	    }
	}

	free(data16);
	
// append to log file
	FILE		*f;
	if((f = fopen("cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(f,"%s ",argv[i]);
	fprintf(f,"\n");
	fclose(f);

	return 0;
}
