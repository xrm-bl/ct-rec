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
#define MAX_SHOT			50010

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
unsigned short	*data, *dark, *I;
double	*I0, *I0sum, *I0sum1;

// image profile from 'output.log'
float	shottime[MAX_SHOT], shotangle[MAX_SHOT];
int		Ishot[MAX_SHOT], NST, NI0, II0[MAX_SHOT];

// flag for q001.img or q0001.img
int		iFlag;

// projection
short	N, M, n_total;
short	Nx,Ny;


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
	TIFFSetField(image, TIFFTAG_ARTIST, "ict_prj");
//	TIFFSetField(image, TIFFTAG_MINSAMPLEVALUE, mmmin );
//	TIFFSetField(image, TIFFTAG_MAXSAMPLEVALUE, mmmax );

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data32 + i*wX, wX * sizeof(float));
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
	short		nnn;
	FILE		*f;
	char		lne[100], *ss;
	int			flg_I0;
	char		outlog[100];


// open parameter file
//	fprintf(stderr,"read %s/output.log",dirin);
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
	n_total = (short)(nnn + 2);
//	fprintf(stderr, " nshot = %d, NI0 = %d, total = %d \n", NST, NI0, n_total);
	M=NST;
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
	Header		h;
	char		path[2048];
	char		fname[20];
	double		I01, I02;

	double		*po, *p_sum, *p_ave;
	float		*data32;
	double		mmmin, mmmax, XXX;
	char		*comm = NULL;

// p initialization
	po = (double *)malloc(Nx*Ny*sizeof(double));
	p_sum = (double *)malloc(Ny*sizeof(double));
	p_ave = (double *)malloc(Ny*sizeof(double));
	
	if ((data32 = (float*)malloc(sizeof(float)*Nx*Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		return(1);
	}

	ilp=(int *)malloc(NST*sizeof(int));
	iplc=0;

// counting for number of projections (NST = nshot)
	nshot = 0;

// loop between I0_1st and I0_2nd
	for ( j = 0; j < NI0-1; ++j){
//		printf("%d\n",j);

// make I0 from averaging of Is
		for ( k = II0[j] + 1; k < II0[j+1]; ++k){
			if(iFlag==0) sprintf(fname, "%s/q%03d.img", dirin, k);
			if(iFlag==1) sprintf(fname, "%s/q%04d.img", dirin, k);
			if(iFlag==2) sprintf(fname, "%s/q%05d.img", dirin, k);
			if ((i = read_hipic(fname, &h)) != 0){
				printf("something wrong -- return value is %d(II01)", i);
				return(-1);
			}
			printf("%s\r",fname);
			for(i=0;i<Nx*Ny;++i){
//				*(I0sum1+i)=*(I0sum1+i)+*(data+i)-*(dark+i);
				*(I0sum+i)=*(I0sum+i)+*(data+i)-*(dark+i);
			}
//			for(jy=2;jy<Ny-2;++jy){
//				for (jx=0;jx<Nx;++jx){
//					*(I0sum+Nx*jy+jx)=(0.05**(I0sum1+Nx*(jy-2)+jx)+0.15**(I0sum1+Nx*(jy-1)+jx)+0.6**(I0sum1+Nx*(jy-0)+jx)+0.15**(I0sum1+Nx*(jy+1)+jx)+0.05**(I0sum1+Nx*(jy+1)+jx));
//				}
//			}
			free(data);
		}
		for ( k = II0[j] + 1; k < II0[j+1]; ++k){
			if(iFlag==0) sprintf(fname, "%s/q%03d.img", dirin, k);
			if(iFlag==1) sprintf(fname, "%s/q%04d.img", dirin, k);
			if(iFlag==2) sprintf(fname, "%s/q%05d.img", dirin, k);
			if ((i = read_hipic(fname, &h)) != 0){
				printf("something wrong -- return value is %d(II01)", i);
				return(-1);
			}
			for(i=0;i<Nx*Ny;++i){
				*(I+i)=*(data+i)-*(dark+i);
			}
			free(data);
			
//			printf("%d\n",k);
			*(ilp+nshot)=0;
			for(jx=0;jx<Nx*Ny;++jx){
				*(I0+jx)=*(I0sum+jx)/(double)(II0[j+1] - (II0[j] + 1));
				if ((*(I+jx)) > 1){
					*(po+jx)=log((double)*(I0+jx)/(double)(*(I+jx)));
				}else{
					*(ilp+nshot)=1;
				}
			}

// roupe for correction
			if(*(ilp+nshot)==0){
				for(jy=0;jy<Ny;++jy){
					*(p_sum+jy) = 0.0;
					*(p_ave+jy) = 0.0;
				}
				for(jy=0;jy<Ny;++jy){
					for (jx=0; jx<10; ++jx){
						*(p_sum+jy) = *(p_sum+jy) + *(po+Nx*jy+jx) + *(po+Nx*jy+Nx-(jx+1));
					}
				}
				for(jy=0;jy<Ny;++jy){
					*(p_ave+jy) = *(p_sum+jy)/20.;
					for (jx=0; jx<Nx; ++jx){
						*(po+Nx*jy+jx) = *(po+Nx*jy+jx) - *(p_ave+jy); /* comment out here */
					}
				}
			}
			if(*(ilp+nshot)==1){
				for(jx=0;jx<Nx*Ny;++jx){
					*(po+jx)=0.0;
				}
			}
			nshot = nshot + 1;
//			if((nshot%100)==0) printf("%d\r",nshot);
			
			mmmin=100.0;
			mmmax=-100.0;
			jx=0;
			for(y=0;y<Ny;y++){
				for(x=0;x<Nx;x++){
					XXX=*(po+jx);
					*(data32+jx)=(float)XXX;
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
			Store32TiffFile(path, Nx, Ny, 32, data32, comm);

			free(comm);
			
		} // end of k loop
	} // end of j loop
	printf("\n");

// Store I0
	for(jx=0;jx<Nx*Ny;++jx){
		*(data32+jx)=(float)*(I0+jx);
	}
	(void)sprintf(path, "ict-I0.tif");
	(void)sprintf(comm,"butter I0");
	Store32TiffFile(path, Nx, Ny, 32, data32, comm);

	free(po);free(p_sum);free(p_ave);free(data32);
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
//	printf("%s/dark.img\n",argv[1]);

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

	I     = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
	I0    = (double *) malloc(Nx*Ny*sizeof(double));
	I0sum = (double *) malloc(Nx*Ny*sizeof(double));
	I0sum1= (double *) malloc(Nx*Ny*sizeof(double));
	fprintf(stderr, "Nx, Ny= %d %d \n", Nx,Ny);

	if((i=StoreProjection(argv[1], argv[2])) !=0){
		printf("something wrong in StoreProjection (%d)\n",i);
		return(1);
	}
	fprintf(stderr, "done. \n");

	free(I);free(I0);free(I0sum);free(I0sum1);

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
