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

void Store8TiffFile(char *wname, int wX, int wY, int wBPS, unsigned char *data8, char *wdesc)
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
		TIFFWriteRawStrip(image, i, data8 + i*wX, wX * sizeof(unsigned char));
	}

	TIFFClose(image);
}

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
	unsigned short	*data16;
	unsigned char	*data8;
	int		m, n;
	FOM	**F;

	
	//	printf("%d\n", argc);
	if (argc == 4 || argc == 6 || argc == 10) {
		l_sta = -1;
		l_dst = -1;
		d_min = 10000.;
		d_max = -10000.;
		for (i = 0; i<99999; i++) {
			sprintf(fh, "%s/rec%05ld.tif", argv[2], i);
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
				sscanf(desc, "%f\t%f\t%d\t%f\t%lf\t%lf", &dmd, &dmd, &dml, &dmd, &dmmin, &dmmax);
				if (d_min>dmmin) d_min = dmmin;
				if (d_max<dmmax) d_max = dmmax;
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
		if (argc == 6) {
			d_min = atof(argv[4]);
			d_max = atof(argv[5]);
			if (d_min>d_max) Error("d range");
		}
		if (argc == 10) {
			d_min = atof(argv[4]);
			d_max = atof(argv[5]);
			x1 = atoi(argv[6]);
			y1 = atoi(argv[7]);
			x2 = atoi(argv[8]);
			y2 = atoi(argv[9]);
			if (x1>x2) Error("x range");
			if (y1>y2) Error("y range");
		}
	}
	else {
		fputs("usage : tif_f2i bit rec/ out/ (LACmin LACmax) (x1 y1 x2 y2)\n", stderr);
		exit(1);
	}

	dBPS = atoi(argv[1]);

	printf("%ld\t%ld\t%d\t%lf\t%lf\n", l_sta, l_dst, dBPS, d_min, d_max);
	//	printf("%lf\t%lf\n", div, base);

	if (dBPS == 0)  return 0;
	if (dBPS == 8)  swid = 255;
	if (dBPS == 16) swid = 65535;

	div = ((double)swid) / (d_max - d_min);
	base = d_min;

	cNx = x2 - x1 + 1; if(cNx%2==1) {cNx=cNx-1; x2=x2-1;}
	cNy = y2 - y1 + 1; if(cNy%2==1) {cNy=cNy-1; y2=y2-1;}

	data16 = (unsigned short*)malloc(sizeof(unsigned short)*cNx*cNy);
	data8  = (unsigned char*)malloc(sizeof(unsigned char)*cNx*cNy);

	for (m = 0; m < cNy; m++) {
		for (n = 0; n < cNx; n++) {
			*(data16 + m*cNx + n) = 0;
			*(data8 + m*cNx + n) = 0;
		}
	}

	for (i = l_sta; i < l_dst + 1; ++i) {
		sprintf(fh, "%s/rec%05ld.tif", argv[2], i);
		fprintf(stderr, "%s(32bit float) -> ", fh);
		Read32TiffFile(fh,0);
//		if (ReadImageFile_Float(fh,&Nx,&Ny,&F,&desc))
//	  			(void)fprintf(stderr, "%s : containing non-float pixel values (warning).\n", fh);

		j=0;
		if(dBPS==8){
			for (m = y1; m <= y2; m++) {
				for (n = x1; n <= x2; n++) {
					val = (long)(div*((*(data32 + m*Nx + n)) - base));
//					val = (long)(div*((F[m][n]) - base));
					if (val<1)    val = 0;
					if (val>255) val = 255;
					*(data8 + j) = (unsigned char)(val);
					j=j+1;
				}
			}
		}else{
			for (m = y1; m <= y2; m++) {
				for (n = x1; n <= x2; n++) {
					val = (long)(div*((*(data32 + m*Nx + n)) - base));
//					val = (long)(div*((F[m][n]) - base));
					if (val<1)    val = 0;
					if (val>65535) val = 65535;
					*(data16 + j) = (unsigned short)(val);
					j=j+1;
				}
			}
		}

		//		printf("%s\n",desc);
		sprintf(wdesc, "%s\t%lf\t%lf", desc, d_min, d_max);
//		fprintf(stderr, "%s\t%lf\t%lf\n", desc, d_min, d_max);
		//		printf("%s\t%d\t%d\t%d\t",fo,cNx,cNy,dBPS);
		//		printf("%s\n",desc,div,base);

		if (dBPS == 8){
			sprintf(fo, "%s/ro%05ld.tif", argv[3], i);
			Store8TiffFile(fo, cNx, cNy, dBPS, data8, wdesc);
		}
		else {
			sprintf(fo, "%s/rh%05ld.tif", argv[3], i);
			Store16TiffFile(fo, cNx, cNy, dBPS, data16, wdesc);
		}

//		StoreTiffFile(fo, cNx, cNy, dBPS, data16, data8, wdesc);
		fprintf(stderr, "%s(%d bit)\r", fo, dBPS);

		free(data32);
//	    for (j=Ny-1; j>=0; j--) free(F[j]); free(F);

	}
	printf("\nfinish.\n");


	free(data16); free(data8);

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