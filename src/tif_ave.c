
#include <stdio.h>
#include <stdlib.h>
#include "tiffio.h"
#include "cell.h"
#include "rif.h"

#define INTEL

static unsigned short *data;
static int Nx, Ny, BPS;
static char *desc;

static void Error(msg)
char        *msg;
	{
		fputs(msg,stderr);
		fputc('\n',stderr);
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
	TIFFSetField(image, TIFFTAG_ARTIST, "tif_f2i_libtiff");

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data16 + i*wX, wX * sizeof(unsigned short));
	}

	TIFFClose(image);
}

void Read16TiffFile(char* rname, int iHead)
{
	TIFF* image;
	long i, j;
	unsigned short spp, com, pm, pc, rps;
	unsigned short* rline;

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
	if ((data = (unsigned short*)malloc(sizeof(unsigned short) * Nx * Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		exit(1);
	}
	if ((rline = (unsigned short*)_TIFFmalloc(TIFFScanlineSize(image))) == NULL) {
		printf("cannot allocate memory for line scan\n");
		exit(1);
	}

	for (i = 0; i < Ny; i++) {
		if (TIFFReadScanline(image, rline, i, 0) < 0) {
			printf("cannot get tif line -> %ld\n", i);
			exit(1);
		}
		for (j = 0; j < Nx; j++) {
			*(data + i * Nx + j) = *(rline + j);
		}
	}

	_TIFFfree(rline);

	TIFFClose(image);
}

/*----------------------------------------------------------------------*/

int	main(argc,argv)
int	argc;
char	**argv;
{
	long	vv,hh,jj,i;
	long 	*sumdata;
	unsigned short *outimg;
	Cell	**cell;

	if (argc<4) {
	    fputs("usage : tif_ave input1.tif input2.tif ... output.tif\n",stderr);
	    return (1);
	}

//	printf("%d\n", argc-2);
// read input images
	for(i=1;i<argc-1;++i){
//		Read16TiffFile(argv[i],0);
	    ReadImageFile(argv[i],&Nx,&Ny,&BPS,&cell,&desc);
		if (i==1){
			sumdata = (long *) malloc(Nx*Ny*sizeof(long));
			for(jj=0;jj<Nx*Ny;++jj) *(sumdata+jj)=0;
		}
		printf("read %s\t%d\t%d\r",argv[i],Nx,Ny);
		jj=0;
		for(vv=0;vv<Ny;++vv){
			for(hh=0;hh<Nx;++hh){
				*(sumdata+jj)=*(sumdata+jj)+cell[vv][hh];
//			printf("%d\t%d\t%d\t%d\t%d\n", jj, hh, vv, *(sumdata+jj), cell[vv][hh]);
				jj=jj+1;
			}
		}
		//free(data);
		for (vv=0; vv<Ny ; vv++) free(cell[vv]);
		free(cell);
	}


// initialize outimg
	outimg = (unsigned short *)malloc(Nx*Ny*sizeof(short));
	jj=0;
	for (vv=0; vv<Ny; vv++){			// vertical loop (in reconstructed image)
		for (hh=0; hh<Nx; hh++){		// horizontal loop
//			*(outimg+jj) = (unsigned short)(*(sumdata+jj));
			*(outimg+jj) = (unsigned short)(*(sumdata+jj)/(argc-2));
//			printf("%d\t%d\t%d\t%d\t%d\n", jj, hh, vv, *(sumdata+jj), cell[vv][hh]);
			jj=jj+1;
		}
	}

	BPS=16;
	Store16TiffFile(argv[argc-1], Nx, Ny, BPS, outimg, desc);
	printf("\nwrite %s \n",argv[argc-1]);

	return 0;
}
