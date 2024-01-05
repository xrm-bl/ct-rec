#include <stdio.h>
#include <stdlib.h>
#include "tiffio.h"
#include "rif_f.h"

#define MA(cnt,ptr)	malloc((cnt)*sizeof(*(ptr)))

static float *data32;
static int Nx, Ny, BPS;
static char *desc;

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
	TIFFSetField(image, TIFFTAG_ARTIST, "rec_gf");
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

int	main(int argc, char *argv[])
{
	int		dBPS, x1, y1, x2, y2, x, y, x0, y0, cNx, cNy;
	long	i, j, dml;
	long	l_sta, l_dst;
	double	d_min, d_max, ccc, dmd, dmmin, dmmax;
	long	val;
	double	base, div;
	char	fh[25], fo[25], wdesc[300];
	unsigned short	swid;
	float	*dataout;
	int		m, n;
	FOM	**F;

	
	//	printf("%d\n", argc);
	if (argc == 3 || argc == 7) {
		l_sta = -1;
		l_dst = -1;
		for (i = 0; i<99999; i++) {
			sprintf(fh, "%s/rec%05ld.tif", argv[1], i);
//						fprintf(stderr,"%s\r",fh);
			if (existFile(fh)) {
				if (l_sta == -1) {
					l_sta = i;
				}
				else {
					l_dst = i;
				}
//				Read32TiffFile(fh,1);
				if (ReadImageFile_Float(fh,&Nx,&Ny,NULL,&desc))
	    			(void)fprintf(stderr, "%s : containing non-float pixel values (warning).\n", fh);
//				sscanf(desc, "%f\t%f\t%d\t%f\t%lf\t%lf", &dmd, &dmd, &dml, &dmd, &dmmin, &dmmax);
				fprintf(stderr, "%s\r", fh);
			}
		}
		fprintf(stderr, "\n");
		//	fprintf(stderr,"\n%s\n",desc);
		//	fprintf(stderr,"%d\t%d\t%lf\t%lf\n", l_sta, l_dst, d_min, d_max);
		x1 = 0;
		y1 = 0;
		x2 = Nx-1;
		y2 = Ny-1;
		if (argc == 7) {
			x1 = atoi(argv[3]);
			y1 = atoi(argv[4]);
			x2 = atoi(argv[5]);
			y2 = atoi(argv[6]);
			if (x1>x2) Error("x range");
			if (y1>y2) Error("y range");
			if (x2>Nx-1) Error("x2 > Nx");
			if (y2>Ny-1) Error("y2 > Ny");
		}
	}
	else {
		fputs("usage : rec_crop in/ out/ x1 y1 x2 y2\n", stderr);
		exit(1);
	}

	printf("%ld\t%ld\n", l_sta, l_dst);
	//	printf("%lf\t%lf\n", div, base);

	cNx = x2 - x1 + 1; if(cNx%2==1) {cNx=cNx-1; x2=x2-1;}
	cNy = y2 - y1 + 1; if(cNy%2==1) {cNy=cNy-1; y2=y2-1;}

	dataout = (float*)malloc(sizeof(float)*cNx*cNy);

	for (i = l_sta; i < l_dst + 1; ++i) {
		sprintf(fh, "%s/rec%05ld.tif", argv[1], i);
		fprintf(stderr, "%s -> ", fh);
		Read32TiffFile(fh,0);
//		if (ReadImageFile_Float(fh,&Nx,&Ny,&F,&desc))
//	  			(void)fprintf(stderr, "%s : containing non-float pixel values (warning).\n", fh);

		j=0;
		for (m = y1; m <= y2; m++) {
			for (n = x1; n <= x2; n++) {
				*(dataout + j) = *(data32 + m*Nx + n);
				j=j+1;
			}
		}

		sprintf(fo, "%s/rec%05ld.tif", argv[2], i);
		Store32TiffFile(fo, cNx, cNy, 32, dataout, desc);

		fprintf(stderr, "%s\r", fo);

		free(data32);

	}
	printf("\nfinish.\n");
	free(dataout);

	// append to log file
	FILE		*f;
	if ((f = fopen("cmd-hst.log", "a")) == NULL) {
		return(-1);
	}
	for (i = 0; i<argc; ++i) fprintf(f, "%s ", argv[i]);
	fprintf(f, "\n");
	fclose(f);

	return 0;
}

