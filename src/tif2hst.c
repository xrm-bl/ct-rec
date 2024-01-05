#include <stdio.h>
#include <stdlib.h>
#include "tiffio.h"

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
	double	d_min, d_max, ccc, dmd, dmmin, dmmax, wd;
	long	val, iwd;
	double	base, div;
	char	fh[25], fo[25], wdesc[300];
	unsigned short	swid;
	unsigned long	*data;
	int		m, n;

	//	printf("%d\n", argc);
	if (argc == 2 ) {
		l_sta = -1;
		l_dst = -1;
		d_min = 10000.;
		d_max = -10000.;
		for (i = 0; i<99999; i++) {
			sprintf(fh, "%s/rec%05ld.tif", argv[1], i);
			//			fprintf(stderr,"%s\t",fh);
			if (existFile(fh)) {
				if (l_sta == -1) {
					l_sta = i;
				}
				else {
					l_dst = i;
				}
				Read32TiffFile(fh,1);
				sscanf(desc, "%lf\t%lf\t%ld\t%lf\t%lf\t%lf", &dmd, &dmd, &dml, &dmd, &dmmin, &dmmax);
				if (d_min>dmmin) d_min = dmmin;
				if (d_max<dmmax) d_max = dmmax;
				fprintf(stderr, "%s\r", fh);
			}
		}
		fprintf(stderr, "\n");
		x1 = 0;
		y1 = 0;
		x2 = Nx-1;
		y2 = Ny-1;
//		if (argc == 6) {
//			d_min = atof(argv[4]);
//			d_max = atof(argv[5]);
//			if (d_min>d_max) Error("d range");
//		}
		if (argc == 6) {
//			d_min = atof(argv[4]);
//			d_max = atof(argv[5]);
			x1 = atoi(argv[2]);
			y1 = atoi(argv[3]);
			x2 = atoi(argv[4]);
			y2 = atoi(argv[5]);
			if (x1>x2) Error("x range");
			if (y1>y2) Error("y range");
		}
	}
	else {
		fputs("usage : tif2hst rec/ (x1 y1 x2 y2)\n", stderr);
		return 1;
	}

	printf("%ld\t%ld\t%lf\t%lf\n", l_sta, l_dst, d_min, d_max);

	wd=d_max-d_min;
	iwd=1+(int)(100*wd);

	cNx = x2 - x1 + 1;
	cNy = y2 - y1 + 1;

	data  = (unsigned long*)malloc(sizeof(unsigned long)*iwd);

	for (m = 0; m < iwd; m++) {
		*(data + m) = 0;
	}

	for (i = l_sta; i < l_dst + 1; ++i) {
		sprintf(fh, "%s/rec%05ld.tif", argv[1], i);
//		fprintf(stderr, "%s(32bit float) -> ", fh);
		Read32TiffFile(fh,0);

		j=0;
			for (m = y1; m <= y2; m++) {
				for (n = x1; n <= x2; n++) {
					val = (int)(100*(*(data32 + m*Nx + n)-d_min));
					*(data + val) = *(data + val) + 1;
					j=j+1;
				}
			}

		//		printf("%s\n",desc);
//		fprintf(stderr, "%s\t%lf\t%lf\n", desc, d_min, d_max);
		//		printf("%s\t%d\t%d\t%d\t",fo,cNx,cNy,dBPS);
		//		printf("%s\n",desc,div,base);


		free(data32);

	}
//	fprintf(stderr,"\nfinish.\n");

	for (m = 0; m < iwd; m++) {
		printf("%8.2f, %ld\n", (m+1)/100.+d_min,*(data + m));
	}

	free(data);

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