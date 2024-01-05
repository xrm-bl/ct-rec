
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fft.h"
#include "tiffio.h"
#include "rif_f.h"

#define MA(cnt,ptr)	malloc((cnt)*sizeof(*(ptr)))

static float *data32;
static int Nx, Ny, Nz, BPS;
static char *desc;


#ifndef	COMPLEX
typedef	struct {
		float	r,i;
	} COMPLEX;
#endif

static void Error(char *msg)
{
	fputs(msg, stderr);
	fputc('\n', stderr);
	exit(1);
}

int existFile(const char* path)
{
	FILE* fp = fopen(path, "r");
	if (fp == NULL) {
		return 0;
	}

	fclose(fp);
	return 1;
}

void Read32TiffFile(char* rname, int iHead)
{
	TIFF* image;
	long i, j;
	unsigned short spp, com, pm, pc, rps;
	float* rline;

	image = TIFFOpen(rname, "r");

	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &Nx);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, &Ny);
	TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &BPS);
//	TIFFGetField(image, TIFFTAG_COMPRESSION, &com);
//	TIFFGetField(image, TIFFTAG_PHOTOMETRIC, &pm);
//	TIFFGetField(image, TIFFTAG_SAMPLESPERPIXEL, &spp);
//	TIFFGetField(image, TIFFTAG_ROWSPERSTRIP, &rps);
//	TIFFGetField(image, TIFFTAG_PLANARCONFIG, &pc);
	TIFFGetField(image, TIFFTAG_IMAGEDESCRIPTION, &desc);

	if(iHead==1){
		TIFFClose(image);
		return;
	}
	if ((data32 = (float*)malloc(sizeof(float) * Nx * Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		exit(1);
	}
	if ((rline = (float*)_TIFFmalloc(TIFFScanlineSize(image))) == NULL) {
		printf("cannot allocate memory for line scan\n");
		exit(1);
	}

	for (i = 0; i < Ny; i++) {
		if (TIFFReadScanline(image, rline, i, 0) < 0) {
			printf("cannot get tif line -> %ld\n", i);
			exit(1);
		}
		for (j = 0; j < Nx; j++) {
			*(data32 + i * Nx + j) = *(rline + j);
		}
	}

	_TIFFfree(rline);

	TIFFClose(image);
}

void Store32TiffFile(char *wname, int wX, int wY, int wBPS, float *data32, char *wdesc)
{
	TIFF *image;
	long i;

	image = TIFFOpen(wname, "w");

	TIFFSetField(image, TIFFTAG_IMAGEWIDTH, wX);
	TIFFSetField(image, TIFFTAG_IMAGELENGTH, wY);
	TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(image, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP );
	TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, 1);
	TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, wdesc);
	TIFFSetField(image, TIFFTAG_ARTIST, "rec_gf");
//	TIFFSetField(image, TIFFTAG_MINSAMPLEVALUE, mmmin );
//	TIFFSetField(image, TIFFTAG_MAXSAMPLEVALUE, mmmax );

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data32 + i*wX, wX * sizeof(float));
	}

	TIFFClose(image);
}

/*----------------------------------------------------------------------*/

//#define LEN	1024

int	main(int argc,char **argv)
{
	int		Lx,Ly,Lz,z,nx,ny,y,x,s;
	double	R,d,PI2R2,y2,x2_y2,g,bias=0.0;
	Complex	*Ex,*Ey,*Ez,*F;
	COMPLEX	***C;
	FOM	**D;

	char	fh[25];
	long	i, j;
	long	l_sta, l_dst;
	double	d_min, d_max, ccc, dmd1, dmd2, dmd3, dmd4, dmmin, dmmax;
	int		dml;
	

	if (argc!=4)
	    Error("usage : rec_gf orgDir radius newDir");

	if ((R=atof(argv[2]))<0.0) Error("bad Gaussian radius.");

	l_sta = -1;
	l_dst = -1;
	for (i = 0; i<99999; i++) {
#ifdef WINDOWS
		sprintf(fh, "%s\\rec%05ld.tif", argv[1], i);
#else
		sprintf(fh, "%s/rec%05ld.tif", argv[1], i);
#endif
//			fprintf(stderr,"%s\t%d\r",fh,existFile(fh));
		if (existFile(fh)==1) {
			if (l_sta == -1) {
				l_sta = i;
			}
			else {
				l_dst = i;
			}
//			Read32TiffFile(fh,1);
			if (ReadImageFile_Float(fh,&Nx,&Ny,NULL,&desc))
	   			(void)fprintf(stderr, "%s : containing non-float pixel values (warning).\n", fh);
//			fprintf(stderr, "%s\r", fh);
			fprintf(stderr, "%s\t%s\r", fh, desc);
		}
	}
	if ((data32 = (float*)malloc(sizeof(float) * Nx * Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		exit(1);
	}
	fprintf(stderr, "\n");
	
	nx=Nx;
	ny=Ny;
	
	Nz=l_dst-l_sta+1;

	if ((Lx=SetUpFFT(Nx,&Ex))==0 ||
	    (Ly=SetUpFFT(Ny,&Ey))==0 ||
	    (Lz=SetUpFFT(Nz,&Ez))==0 ||
	(F=(Complex *)malloc(sizeof(Complex)*((Lx>Ly)?(Lx>Lz)?Lx:Lz:
						      (Ly>Lz)?Ly:Lz)))==NULL ||
	    (C=(COMPLEX ***)malloc(sizeof(COMPLEX **)*Nz))==NULL ||
	    (*C=(COMPLEX **)malloc(sizeof(COMPLEX *)*Nz*Ly))==NULL ||
	    (**C=(COMPLEX *)malloc(sizeof(COMPLEX)*Nz*Ly*Lx))==NULL)
	    Error("memory allocation error.");

	fprintf(stderr,"%d\t%d\t%d\t%d\t%d\t%d\n",Nx,Ny,Nz,Lx,Ly,Lz);
	d=1.0/((double)Lx*(double)Ly*(double)Lz);
	for (z=0; z<Nz; z++) {
	    if (z>0) C[z]=C[z-1]+Ly;

#ifdef WINDOWS
		sprintf(fh, "%s\\rec%05ld.tif", argv[1], z+l_sta);
#else
		sprintf(fh, "%s/rec%05ld.tif", argv[1], z+l_sta);
#endif
		Read32TiffFile(fh,0);
//		if (ReadImageFile_Float(fh,&Nx,&Ny,&D,&desc))
//	  			(void)fprintf(stderr, "%s : containing non-float pixel values (warning).\n", fh);
		fprintf(stderr, "%s\r", fh);
//		fprintf(stderr, "\n");

	    for (y=0; y<ny; y++) {
//		fprintf(stderr, "%d\r", y);
			if (z>0 || y>0) C[z][y]=C[z][y-1]+Lx;

			for (x=0; x<nx; x++) {
//			    F[x].r=((double)D[y][x])*d; F[x].i=0.0;
			    F[x].r=((double)*(data32 + y*nx + x))*d; F[x].i=0.0;
			}
			for (; x<Lx; x++) F[x].r=F[x].i=0.0;

			FFT(-1,Lx,Ex,F);
			for (x=0; x<Lx; x++) {
			    C[z][y][x].r=F[x].r; C[z][y][x].i=F[x].i;
			}
		}
//		free(data32);

//		fprintf(stderr, "\n");
	    for (; y<Ly; y++) C[z][y]=C[z][y-1]+Lx;

	    for (x=0; x<Lx; x++) {
//			fprintf(stderr, "%d\r", x);
			for (y=0; y<ny; y++) {
			    F[y].r=C[z][y][x].r; F[y].i=C[z][y][x].i;
			}
			for (; y<Ly; y++) F[y].r=F[y].i=0.0;

			FFT(-1,Ly,Ey,F);
			for (y=0; y<Ly; y++) {
			    C[z][y][x].r=F[y].r; C[z][y][x].i=F[y].i;
			}
	    }
//		fprintf(stderr, "\n");
	}
//	fprintf(stderr,"\n");
	PI2R2=M_PI*M_PI*R*R;
	for (y=0; y<Ly; y++) {
		fprintf(stderr,"y: %05d / %05d\r", y, Ly-1);
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
//	fprintf(stderr,"\n");

//	if ((data32 = (float*)malloc(sizeof(float)*Nx*Ny)) == NULL) {
//		printf("cannot allocate memory for input 32bit TIFF image\n");
//		exit(1);
//	}
	for (z=0; z<Nz; z++) {
		fprintf(stderr,"z: %05d / %05d\r", z, Nz-1);
		d_min=0.0;
		d_max=0.0;
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
	    	for (x=0; x<Nx; x++){
	    		*(data32 + y*Nx + x)=(float)F[x].r;
	    		ccc=*(data32 + y*Nx + x);
				if (d_min>ccc) d_min = ccc;
				if (d_max<ccc) d_max = ccc;
	    	}
		}
#ifdef WINDOWS
		(void)sprintf(fh, "%s\\rec%05ld.tif", argv[1], z+l_sta);
#else
		(void)sprintf(fh, "%s/rec%05ld.tif", argv[1], z+l_sta);
#endif
//		Read32TiffFile(fh,1);
			if (ReadImageFile_Float(fh,&Nx,&Ny,NULL,&desc))
	   			(void)fprintf(stderr, "%s : containing non-float pixel values (warning).\n", fh);
		sscanf(desc, "%lf\t%lf\t%d\t%lf\t%lf\t%lf", &dmd1, &dmd2, &dml, &dmd4, &dmmin, &dmmax);
		sprintf(desc, "%f\t%f\t%d\t%f\t%lf\t%lf", (float)dmd1, (float)dmd2, dml, dmd4, d_min, d_max);
//		sprintf(desc,"%s\t%lf\t%lf", desc, d_min, d_max);
//		(void)printf("%s\r",desc);

#ifdef WINDOWS
		(void)sprintf(fh, "%s\\rec%05ld.tif", argv[3], z+l_sta);
#else
		(void)sprintf(fh, "%s/rec%05ld.tif", argv[3], z+l_sta);
#endif
		Store32TiffFile(fh, Nx, Ny, 32, data32, desc);
//		fprintf(stderr,"%s\r",fh);

	}
	fprintf(stderr,"\n");

	// append to log file
	FILE		*f;
	if((f = fopen("cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(f,"%s ",argv[i]);
	fprintf(f,"\n");
	fclose(f);

	free(**C); free(*C); free(C); free(F); free(Ez); free(Ey); free(Ex);
    free(data32); free(desc);

	return 0;
}
