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
static Float	*data32;

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
	char		dirin[25], dirout[25], wdesc[300];
	char		fh[25], fo[25];
	int			z1, z2;
	Float		*po, f_temp;
	int			vv, hh, jx;
	char		*comm = NULL;
	double		data_max, data_min;
	double		Dr,RC,RA0,Ct;
	double		Clock, t1,t2,t3;					// timer setting
	FILE		*fs[8192];
	
	if (_setmaxstdio(8192)!=8192)
		Error("set max stdio error");
	
	//	printf("%d\n", argc);
	Clock=CLOCK();
	if (argc == 6 || argc == 9 ) {
		p_sta = -1;
		p_dst = -1;
		for (i = 1; i<100000; i++) {
			sprintf(fh, "%s/p%05d.tif", argv[1], i);
			if (existFile(fh)) {
				if (p_sta == -1) {
					p_sta = i;
					Read32TiffFile(fh,1);
					fprintf(stderr, "%s\r", fh);
				}
				else {
					p_dst = i;
				}
			}
		}
		Nt=p_dst-p_sta+1;
		z1=0;
		z2=Ny;
		Ct=0.0;
		if (argc == 6) {
			Dr = atof(argv[3]);
			RC = atof(argv[4]);
			RA0 = atof(argv[5]);
		}
		if (argc == 9) {
			Dr = atof(argv[3]);
			z1 = atoi(argv[4]);
			RC = atof(argv[5]);
			z2 = atoi(argv[6]);
			Ct = (atof(argv[7])-RC)/(double)(z2-z1);
			RA0 = atof(argv[8]);
		}
	}
	else {
		fprintf(stderr, "usage : p_rec p/ rec/ Dr RC RA0 \nusage : p_rec p/ rec/ Dr L1 C1 L2 C2 RA0");
		return 1;
	}
	printf("%d\t%d\t%d\t%d\t%d\n", Nx, Ny, Nt,z1,z2);

//open output files
	for(m=z1;m<z2;++m){
		sprintf(fh,"s%05d.sin", (int)(m));
		if((fs[m] = fopen(fh,"wb")) == NULL){ 
			printf("can not open %s for output\n", fh);
			return(2);
		}
	}

	if ((data32 = (float *)malloc(sizeof(float) * Nx * Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		return 1;
	}
	
	// store p-data from float tiff files
	for(i=0;i<Nt;++i) {
		sprintf(fh, "%s/p%05d.tif", argv[1], i+1);
		fprintf(stderr, "\rread:\t%s\t", fh);
		(void)Read32TiffFile(fh,0);

		for(m=z1;m<z2;m++) {
			if(j=fwrite(data32+Nx*m,sizeof(float),Nx,fs[m])>Nx*sizeof(float)){
				fprintf(stderr, " Error writing data to %s at %d\n", fh, j);
				fclose(fs[m]);
				return(j+1);
			}
		}
	}
	for(m=z1;m<z2;++m){
		fclose(fs[m]);
	}
	t1=CLOCK()-Clock;
	fprintf(stderr, "\t%lf\n",t1);

	free(data32);
	
//	return 0;
	
	// initilaize for CBP
	if ((P=InitCBP(Nx,Nt))==NULL){
		Error("memory allocation error.");
	}

	po=(float *)malloc(Nx*Nt*sizeof(float));
	data32=(float *)malloc(Nx*Nx*sizeof(float));

	// loop z1 to z2
	
	t1=CLOCK();
		for(m=z1;m<z2;++m){
			sprintf(fh,"s%05d.sin", (int)(m));
			if((fs[m] = fopen(fh,"rb")) == NULL){ 
				printf("can not open %s for output\n", fh);
				return(2);
			}
			if ((j = fread(po,sizeof(float),Nx*Nt,fs[m]))!=Nx*Nt){
				printf("error in reading sinogram (%d)\n", (int)j);
				return(-2);
			}
			fclose(fs[m]);
			if(remove(fh)==-1) printf( "%s cannot be removed.",fh);

			for (i=0;i<Nt;i++){
				for (n=0; n<Nx; n++){
//					if(abs(*(po+Nx*Ny*i+m*Nx+n))<100.){
						P[i][n]=*(po+i*Nx+n);
//					}
				}
			}
// CBP
			Clock=CLOCK();
			F=CBP(1.0,-RC,RA0);
			t2=CLOCK()-Clock;

// Store CT images
			data_max =-32000.;
			data_min = 32000.;

			Clock=CLOCK();
			for(vv=0; vv<Nx; vv++){
				for (hh=0; hh<Nx; hh++){
					*(data32+Nx*vv+hh) = F[vv][hh]*10000./Dr;	/* unit change  um -> cm */;
					if(data_max<*(data32+Nx*vv+hh)) data_max=*(data32+Nx*vv+hh);
					if(data_min>*(data32+Nx*vv+hh)) data_min=*(data32+Nx*vv+hh);
	printf("\n%d\t%d\t%f\t%f",hh,vv,F[vv][hh],*(data32+Nx*vv+hh)*Dr/10000.);
				}
			}
	printf("\n");
return 0;
			if ((comm=(char *)malloc(150))==NULL){
				Error("comment memory allocation error.");
			}
			sprintf(comm,"%f\t%f\t%d\t%f\t%f\t%f",Dr, RC, Nt, RA0, (float)data_min, (float)data_max);
			sprintf(fo, "%s/rec%05d.tif", argv[2], m);
			(void)Store32TiffFile(fo, Nx, Nx, 32, data32, comm);
			 free(comm);
			t3=CLOCK()-Clock;
			fprintf(stderr, "\rstore:\t%s/rec%05d.tif\t%lf\t%lf", argv[2], m, t2, t3);

			RC=RC+Ct;
		}

	printf("\t%f\n",CLOCK()-t1);
	printf("finish. \n");
	free(po);free(data32);

	// append to log file
	FILE		*ff;
	if ((ff = fopen("cmd-hst.log", "a")) == NULL) {
		return(-1);
	}
	for (i = 0; i<argc; ++i) fprintf(ff, "%s ", argv[i]);
	fprintf(ff, "\n");
	fclose(ff);

	return 0;
}