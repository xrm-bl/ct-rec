// program ct_sino
// 
// Required files are q???.img, dark.img.
// output file "s????.sin"

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <string.h>
#include <time.h>
//#include "sif_f.h"
#include "tiffio.h"

/*----------------------------------------------------------------------*/
#ifndef M_PI
#define M_PI				3.1415926535897932385
#endif
#define INTEL
#define HiPic_Header_Size	64
#define MAX_SHOT			90010

/*----------------------------------------------------------------------*/
struct HiPic_Header{
	char			head[2];
	short			comment_length;
	short			width;
	short			height;
	short			x_offset;
	short			y_offset;
	short			type;
	char			reserved[63 - 13];
	char			*comment;
};
typedef struct HiPic_Header Header;

/*----------------------------------------------------------------------*/

// main data for read transmitted images (data[y][x]) 772
unsigned short	*data, *dark, *II01, *II02, *I;
double	*I0;

// image profile from 'output.log'
float	shottime[MAX_SHOT], shotangle[MAX_SHOT];
unsigned short	Ishot[MAX_SHOT], NST, NI0, II0[MAX_SHOT];

// flag for q001.img or q0001.img
int		iFlag;

// projection
unsigned short	N, n_total;

short Nx,Ny;


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
	TIFFSetField(image, TIFFTAG_ARTIST, "ct_prj");
//	TIFFSetField(image, TIFFTAG_MINSAMPLEVALUE, mmmin );
//	TIFFSetField(image, TIFFTAG_MAXSAMPLEVALUE, mmmax );

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data32 + i*wX, wX * sizeof(float));
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


/*----------------------------------------------------------------------*/

static void Error(msg)
char        *msg;
{
	fputs(msg,stderr);
	fputc('\n',stderr);
	exit(1);
}


/*----------------------------------------------------------------------*/
int read_hipic(fname, h)
char	*fname;
Header	*h;
{
		short		x_size, y_size;
		int		j;
		FILE	*fi;

//open input files
		if((fi = fopen(fname,"rb")) == NULL)
		{
			printf("can not open %s for input\n", fname);
			free(h->comment); 
			return(-10);
		}
// read comment and image data from file1
		if (fread(h, sizeof(char), HiPic_Header_Size, fi) != HiPic_Header_Size)
		{
			fprintf(stderr, "EOF in %s (header)\n", fname);
			free(h->comment); 
			return(-1);
		}
		if (strncmp(h->head, "IM", 2))
		{
			fclose(fi); 
			free(h->comment); 
			fprintf(stderr, "This is not HiPIc image file.\n"); 
			fprintf(stderr, "File name = %s\n", fname); 
			return(-11);
		}

		h->comment = (char *) malloc(h->comment_length + 1);
		if (h->comment == NULL)
		{
			fclose(fi); 
			free(h->comment); 
			fprintf(stderr, "Connot allocate memory (HiPic_Header.comment)\n"); 
			fprintf(stderr, "File name = %s\n", fname); 
			return(-12);
		}
		if (fread(h->comment, sizeof(char), h->comment_length, fi) != h->comment_length)
		{ 
			fclose(fi); 
			free(h->comment); 
			fprintf(stderr, "EOF found (comment)\n"); 
			fprintf(stderr, "File name = %s\n", fname); 
			return(-13);
		}

		h->comment[h->comment_length] = '\0';

		Nx = h->width;
		Ny = h->height;

		data = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
		if ((j = fread(data, sizeof(unsigned short), Nx*Ny, fi)) != Nx*Ny){
			fclose(fi);
			free(h->comment); 
			fprintf(stderr, " Error reading %s at %d\n", fname, j);
			fprintf(stderr, "File name = %s\n", fname); 
			return(1);
		}

// check for correct reading of HiPic file
//printf("-------------------------------\n");
//printf("%s\n", fname);
//printf("%s\n", h->head);
//printf("%d\n", h->comment_length);
//printf("%d\n", h->width);
//printf("%d\n", h->height);
//printf("%d\n", h->x_offset);
//printf("%d\n", h->y_offset);
//printf("%d\n", h->type);
//printf("%s\n", h->reserved);
//printf("%s\n", h->comment);
//printf("-------------------------------\n");

		fclose(fi);
		return(0);
	}

/*----------------------------------------------------------------------*/

int read_log(char *dirin)
{
	int			i, j;
	unsigned short		nnn;
	FILE		*f;
	char		lne[100], *ss;
	int			flg_I0;
	char		outlog[100];


// open parameter file
	fprintf(stderr,"read %s/output.log",dirin);
	(void)sprintf(outlog,"%s/output.log",dirin);
	if((f = fopen(outlog,"r")) == NULL){
		printf("cannot open output.log \n");
		return(-10);
	}

// read parameters
	NI0 = 0;
	NST = 0;
	nnn = 0;
	j   = 0;
	while(j>=0){
		j=j+1;
		lne[0]='\0';
		ss = fgets(lne, 100, f);
		if (i = strlen(lne) > 2){
			sscanf(lne, "%hd %f %f %d", &nnn, &shottime[j], &shotangle[j], &flg_I0);
			if (flg_I0 == 0){
				II0[NI0] = nnn;
				NI0=NI0+1;				// number of I_0
			}
			else{
				Ishot[NST] = nnn;
				NST=NST+1;				// number of I
			}
		}
		else{
			break;
		}
	}
	n_total = (unsigned short)(nnn + 2);
	fprintf(stderr, " nshot = %d, NI0 = %d, total = %d \n", NST, NI0, n_total);
	fclose(f);
	if(n_total<1000) iFlag = 0;
	if(n_total>=1000) iFlag = 1;
	if(n_total>=10000) iFlag = 2;
	if(n_total>MAX_SHOT){
		fprintf(stderr, "Too many projections!");
		return(1);
	}
	return(0);
}

/*----------------------------------------------------------------------*/

int StoreProjection(char *dirin, char *dirout)
{
	int			i, j, k, jx, jy, nshot, *ilp, iplc, x, y;
	double		t1, t2;
	double		*a, *b;
	double		DI1, DI2, DI3, DI4, DI5; 
	Header		h;
	char		path[2048];
	char		fname[20];
	double		I01, I02;

	double		*po;
	double		p_sum,p_ave;
//	FILE		*fi;
//	FOM			**fom;
//	float		*data32;
	unsigned short	*data16;
	double		mmmin, mmmax, XXX;
	char		*comm = NULL;

//		printf("1\n");
	srand((unsigned int)time(NULL));

// p initialization
	po = (double *)malloc(Nx*Ny*sizeof(double));
	a  = (double *)malloc(Nx*Ny*sizeof(double));
	b  = (double *)malloc(Nx*Ny*sizeof(double));

	if ((data16 = (unsigned short *)malloc(sizeof(unsigned short)*Nx*Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		exit(1);
	}

	ilp=(int *)malloc(NST*sizeof(int));
	iplc=0;

// counting for number of projections (NST = nshot)
	nshot = 0;

// loop between I0_1st and I0_2nd
	for ( j = 0; j < NI0-1; ++j){
//		printf("%d\n",j);

//IIO[j] and IIO[j+1] are opened
		if(iFlag==0) sprintf(fname, "%s/q%03d.img", dirin, II0[j]);
		if(iFlag==1) sprintf(fname, "%s/q%04d.img", dirin, II0[j]);
		if(iFlag==2) sprintf(fname, "%s/q%05d.img", dirin, II0[j]);
		if ((i = read_hipic(fname, &h)) != 0){
			printf("something wrong -- return value is %d(II01)", i);
			return(-1);
		}
//		for(i=0;i<Nx*Ny;++i){
//			*(II01+i)=*(data+i);
//		}
		for(jy=0;jy<Ny;++jy){
			for(jx=0;jx<Nx;++jx){
				if((jx>1)&&(jx<Nx-2)){
					DI1=*(data+jy*Nx+jx-2)*0.10;
					DI2=*(data+jy*Nx+jx-1)*0.20;
					DI3=*(data+jy*Nx+jx  )*0.40;
					DI4=*(data+jy*Nx+jx+1)*0.20;
					DI5=*(data+jy*Nx+jx+2)*0.10;
					*(II01+jy*Nx+jx)=DI1+DI2+DI3+DI4+DI5;
				}else{
					*(II01+jy*Nx+jx)=*(data+jy*Nx+jx);
				}
			}
		}
		free(data);
		
		if(iFlag==0) sprintf(fname, "%s/q%03d.img", dirin, II0[j+1]);
		if(iFlag==1) sprintf(fname, "%s/q%04d.img", dirin, II0[j+1]);
		if(iFlag==2) sprintf(fname, "%s/q%05d.img", dirin, II0[j+1]);
		if ((i = read_hipic(fname, &h)) != 0){
			printf("something wrong -- return value is %d(II02)", i);
			return(-1);
		}
//		for(i=0;i<Nx*Ny;++i){
//			*(II02+i)=*(data+i);
//		}
		for(jy=0;jy<Ny;++jy){
			for(jx=0;jx<Nx;++jx){
				if((jx>1)&&(jx<Nx-2)){
					DI1=*(data+jy*Nx+jx-2)*0.10;
					DI2=*(data+jy*Nx+jx-1)*0.20;
					DI3=*(data+jy*Nx+jx  )*0.40;
					DI4=*(data+jy*Nx+jx+1)*0.20;
					DI5=*(data+jy*Nx+jx+2)*0.10;
					*(II02+jy*Nx+jx)=DI1+DI2+DI3+DI4+DI5;
				}else{
					*(II02+jy*Nx+jx)=*(data+jy*Nx+jx);
				}
			}
		}
		free(data);

// I0EV
		t1 = shottime[II0[j]];
		t2 = shottime[II0[j+1]];

// 1 layer(ln)
		for (jx=0;jx<Nx*Ny;++jx){
			I01     = (double)(*(II01+jx)-*(dark+jx));
			I02     = (double)(*(II02+jx)-*(dark+jx));
			*(a+jx) = (double)(((double)(I02   - I01))    / (t2 - t1));
			*(b+jx) = (double)(((double)I01*t2 - (double)I02*t1) / (t2 - t1));
		}
		for ( k = II0[j] + 1; k < II0[j+1]; ++k){
			// obtain p(x) from a[jx], b[jx] using shottime[k]
			if(iFlag==0) sprintf(fname, "%s/q%03d.img", dirin, k);
			if(iFlag==1) sprintf(fname, "%s/q%04d.img", dirin, k);
			if(iFlag==2) sprintf(fname, "%s/q%05d.img", dirin, k);
			if ((i = read_hipic(fname, &h)) != 0){
				printf("something wrong -- return value is %d(II01)", i);
				return(-1);
			}
			for(i=0;i<Nx*Ny;++i){
				*(I+i)=*(data+i);
			}
			free(data);
			
//			printf("%d\n",k);
			*(ilp+nshot)=0;
			for(jx=0;jx<Nx*Ny;++jx){
				*(I0+jx)=(*(a+jx) * shottime[k] + *(b+jx));
				*(po+jx)=*(I+jx)-*(dark+jx)-*(I0+jx);
			}

			nshot = nshot + 1;
//			if((nshot%100)==0) printf("%d\r",nshot);
			
			mmmin=100.0;
			mmmax=-100.0;
			jx=0;
			for(y=0;y<Ny;y++){
				for(x=0;x<Nx;x++){
					XXX=*(po+jx)+5000;
					if(XXX<0) XXX=0;
					if(XXX>65535) XXX=65535;
					*(data16+jx)=(unsigned short)XXX;
					if (mmmin>XXX) mmmin = XXX;
					if (mmmax<XXX) mmmax = XXX;
					jx=jx+1;
				}
			}
			
			if ((comm=(char *)malloc(150))==NULL)
				Error("comment memory allocation error.");

			sprintf(comm,"%lf\t%lf",mmmin, mmmax);
#ifdef WINDOWS
			(void)sprintf(path, "%s\\p%05d.tif", dirout, nshot);
#else
			(void)sprintf(path, "%s/p%05d.tif", dirout, nshot);
#endif
			printf("%s\r",path);
			Store16TiffFile(path, Nx, Ny, 16, data16, comm);

			free(comm);
			
		} // end of k loop
	} // end of j loop
	printf("\n");
	free(a);free(b);free(po);free(data16);
	return (0);
}


/*----------------------------------------------------------------------*/

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC	1000000
extern long clock();
#endif

#define CLOCK()		((double)clock()/(double)CLOCKS_PER_SEC)
#define DESC_LEN	256

int main(argc,argv)
int		argc;
char	**argv;
{
//	double		Clock;					// timer setting
	int			i;
	Header		h;
//	FILE		*fo;
	char		darkfile[100];

// parameter setting
	if (argc!=3){
//		fprintf(stderr, "parameter was wrong!!!\n");
		fprintf(stderr, "usage : %s HiPic/ prj/\n", argv[0]);
//		fprintf(stderr, "default head=q, dark=dark.img\n");
		return(1);
	}

// read shot log file
	if ((i = read_log(argv[1])) != 0){
		fprintf(stderr, "canot read 'output.log' ! (%d)\n", i);
		return(1);
	}
	printf("%s/dark.img\n",argv[1]);

// read dark image
	sprintf(darkfile,"%s/dark.img",argv[1]);
	if ((i = read_hipic(darkfile, &h)) != 0){
		fprintf(stderr, "something wrong in dark file (%d)\n", i);
		return(1);
	}
	dark = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
	for(i=0;i<Nx*Ny;++i){
		*(dark+i)=*(data+i);
//		printf("%d\t%d\t%d\r",i,*(data+i),*(dark+i));
	}
	free(data);
//	fprintf(stderr, "Nx, Ny= %d %d \n", Nx,Ny);

	II01 = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
	II02 = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
	I    = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
	I0   = (double *) malloc(Nx*Ny*sizeof(double));
	fprintf(stderr, "Nx, Ny= %d %d \n", Nx,Ny);

	if((i=StoreProjection(argv[1], argv[2])) !=0){
		printf("something wrong in StoreProjection (%d)\n",i);
		return(1);
	}
	fprintf(stderr, "fin\n");

	free(II01); free(II02);free(I);free(I0);

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
