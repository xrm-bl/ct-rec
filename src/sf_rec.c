#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "tiffio.h"
#include "cbp.h"

#define MA(cnt,ptr)	malloc((cnt)*sizeof(*(ptr)))

static int		Nx, Ny, Nt, M, BPS;
static char		*desc;
static float	*data32;

static void Error(char *msg)
{
	fputs(msg, stderr);
	fputc('\n', stderr);
	exit(1);
}

/*----------------------------------------------------------------------*/

int existFile(const char* path)
{
	FILE* fp = fopen(path, "r");
	if (fp == NULL) {
		return 0;
	}

	fclose(fp);
	return 1;
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
	TIFFSetField(image, TIFFTAG_ARTIST, "p_rec");
//	TIFFSetField(image, TIFFTAG_MINSAMPLEVALUE, mmmin );
//	TIFFSetField(image, TIFFTAG_MAXSAMPLEVALUE, mmmax );

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data32 + i*wX, wX * sizeof(float));
	}

	TIFFClose(image);
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
	if ((rline = (float*)_TIFFmalloc(TIFFScanlineSize(image))) == NULL) {
		printf("cannot allocate memory for line scan\n");
		exit(1);
	}
//		fprintf(stderr, "1\r");

	for (i = 0; i < Ny; i++) {
		if (TIFFReadScanline(image, rline, i, 0) < 0) {
			printf("cannot get tif line -> %d\n", i);
			exit(1);
		}
		for (j = 0; j < Nx; j++) {
			*(data32 + i * Nx + j) = *(rline + j);
		}
	}

	_TIFFfree(rline);
	TIFFClose(image);
	return;
}


float calc_center(){
int		i;
long	j, dx, dx0, N1;
double	mm1, mx1, mm2, mx2;
double	*p000, *p180, rsum, rd, rmsd0, rmsd;
float	center;


// calculate center with RMSD
		p000=(double *)malloc(Nx*sizeof(double));
		p180=(double *)malloc(Nx*sizeof(double));
		for(j=0;j<Nx;++j){
			*(p000+j)=*(data32+0*Nx+j);
			*(p180+j)=*(data32+(Nt-1)*Nx+j);
		}
		
		N1=Nx-1;
		rsum=0.0;
		for(j=0;j<Nx;++j){
			rd=*(p000+j)-*(p180+N1-j);
			rsum=rsum+rd*rd;
		}
		rmsd0=sqrt(rsum/(double)Nx);
		dx0=0;

		for(dx=1;dx<=Nx/2;++dx){
			rsum=0.0;
			for(j=0;j<N1-dx;++j){
				rd=*(p000+dx+j)-*(p180+N1-j);
				rsum=rsum+rd*rd;
			}
			rmsd=sqrt(rsum/(double)(Nx-dx));
			if(rmsd0>rmsd){
				rmsd0=rmsd;
				dx0=dx;
//				printf("%d\t%lf\t%lf\n", dx,rmsd0,rmsd);
			}

			rsum=0.0;
			for(j=0;j<N1-dx;++j){
				rd=*(p000+j)-*(p180-dx+N1-j);
				rsum=rsum+rd*rd;
			}
			rmsd=sqrt(rsum/(double)(Nx-dx));
			if(rmsd0>rmsd){
				rmsd0=rmsd;
				dx0=-dx;
//				printf("%d\t%lf\t%lf\n", -dx,rmsd0,rmsd);
			}
		}
		center=(N1+dx0)/2.0;
//		printf("%d\t%lf\t%lf\n", dx0,rmsd0,center);

	free(p000); free(p180);
	return(center);
}


/*----------------------------------------------------------------------*/

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC	1000000
extern long clock();
#endif

#define CLOCK()		((double)clock()/(double)CLOCKS_PER_SEC)

int	main(int argc, char *argv[])
{
	int			m, n;
	long		i, j, p_sta, p_dst;
	Float		**P, **F;			// full size of reconstructed image
	char		wdesc[300];
	char		fh[25], fo[25];
	int			vv, hh, jx;
	char		*comm = NULL;
	double		data_max, data_min;
	double		Dr,RC,RA0,Ct;
	double		Clock, t1,t2,t3;					// timer setting
//	FILE		*fs[8192];
	
//	if (_setmaxstdio(8192)!=8192)
//		Error("set max stdio error");
	
	//	printf("%d\n", argc);
	Clock=CLOCK();
	if (argc == 3 || argc == 6 ) {
			sprintf(fh, "%s", argv[1]);
			if (existFile(fh)) {
				(void)Read32TiffFile(fh,1);
				fprintf(stderr, "%s\r", fh);
			}
			else {
				fprintf(stderr, "\nNo %s.\n", fh);
				return 1;
				
			}
		Nt=Ny;
		if (argc == 6) {
			Dr = atof(argv[3]);
			RC = atof(argv[4]);
			RA0 = atof(argv[5]);
		}
	}
	else {
			fprintf(stderr, "usage : sf_rec s.tif rec.tif {Dr RC RA0} \n");
		return 1;
	}

	printf("%d\t%d\t%d\t", Nx, Ny, Nt);

	if ((data32 = (float *)malloc(sizeof(float) * Nx * Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		return 1;
	}
	(void)Read32TiffFile(fh,0);

	t1=CLOCK()-Clock;
//	fprintf(stderr, "\t%lf\n",t1);

//	return 0;
	
	// initilaize for CBP
	if ((P=InitCBP(Nx,Nt))==NULL){
		Error("memory allocation error.");
	}

	t1=CLOCK();

		for (i=0;i<Nt;i++){
			for (n=0; n<Nx; n++){
//				if(abs(*(po+Nx*Ny*i+m*Nx+n))<100.){
					P[i][n]=*(data32+i*Nx+n);
//				}
			}
		}
// CBP

	if (argc == 3) {
		Dr = 1.0;
		RC = calc_center();
		RA0 = 0.0;
	}

		fprintf(stderr, "%f\t%f\t%f\n",Dr,RC,RA0);

	Clock=CLOCK();
	F=CBP(1.0,-RC,RA0);
	t2=CLOCK()-Clock;
//	fprintf(stderr, "\t%lf\n",t2);

	free(data32);
	data32=(float *)malloc(Nx*Nx*sizeof(float));

// Store CT images
			data_max =-32000.;
			data_min = 32000.;

			Clock=CLOCK();
			for(vv=0; vv<Nx; vv++){
				for (hh=0; hh<Nx; hh++){
					*(data32+Nx*vv+hh) = F[vv][hh]*10000./Dr;	/* unit change  um -> cm */;
					if(data_max<*(data32+Nx*vv+hh)) data_max=*(data32+Nx*vv+hh);
					if(data_min>*(data32+Nx*vv+hh)) data_min=*(data32+Nx*vv+hh);
//	printf("\n%d\t%d\t%f\t%f",hh,vv,F[vv][hh],*(data32+Nx*vv+hh)*Dr/10000.);
				}
			}
//	printf("\n");
//return 0;
			if ((comm=(char *)malloc(150))==NULL){
				Error("comment memory allocation error.");
			}
			sprintf(comm,"%f\t%f\t%d\t%f\t%f\t%f",Dr, RC, Nt, RA0, (float)data_min, (float)data_max);
			sprintf(fo, "%s", argv[2]);
			(void)Store32TiffFile(fo, Nx, Nx, 32, data32, comm);
			 free(comm);
			t3=CLOCK()-Clock;
			fprintf(stderr, "\rstore:\t%s\t%lf\t%lf", argv[2], t2, t3);


	printf("\t%f\n",CLOCK()-t1);
	printf("finish. \n");
	free(data32);

//	// append to log file
//	FILE		*ff;
//	if ((ff = fopen("cmd-hst.log", "a")) == NULL) {
//		return(-1);
//	}
//	for (i = 0; i<argc; ++i) fprintf(ff, "%s ", argv[i]);
//	fprintf(ff, "\n");
//	fclose(ff);

	return 0;
}