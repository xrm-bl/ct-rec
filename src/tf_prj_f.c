// program tf_prj_f
// 
// Required files are q???.tif, dark.tif.
// output file "p?????.tif"

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <string.h>
#include <time.h>
//#include "sif_f.h"
#include "tiffio.h"
#include "cell.h"
#include "rif.h"

/*----------------------------------------------------------------------*/
#ifndef M_PI
#define M_PI				3.1415926535897932385
#endif
#define INTEL
#define MAX_SHOT			60010


/*----------------------------------------------------------------------*/

// main data for read transmitted images (data[y][x]) 772
unsigned short	*data, *dark, *II01, *II02, *I;
double	*I0, I0cor;
int		cFlag;

// image profile from 'output.log'
float	shottime[MAX_SHOT], shotangle[MAX_SHOT];
int		Ishot[MAX_SHOT], NST, NI0, II0[MAX_SHOT];

// flag for q001.img or q0001.img
int		iFlag;

// projection
unsigned short	N, n_total;

static int Nx, Ny, BPS;
static char *desc;


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
	TIFFSetField(image, TIFFTAG_ARTIST, "tf_prj_f");
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
	n_total = (short)(nnn + 2);
	fprintf(stderr, " nshot = %d, NI0 = %d, total = %d \n", NST, NI0, n_total);
	fclose(f);
	if(n_total<1000) iFlag = 0;
	if(n_total>=1000) iFlag = 1;
	if(n_total>MAX_SHOT){
		fprintf(stderr, "Too many projections!");
		return(1);
	}
	return(0);
}

/*----------------------------------------------------------------------*/

int StoreProjection(char *dirin, char *dirout)
{
	int		i, j, k, jx, jy, nshot, *ilp, iplc, x, y;
	long		li;
	double		t1, t2;
	double		*a, *b;
	char		path[2048];
	char		fname[20];
	double		I01, I02;

	double		*po;
	double		p_sum,p_ave;
//	FILE		*fi;
//	FOM			**fom;
	float		*data32;
	double		mmmin, mmmax, XXX;
	char		*comm = NULL;
	double		II01s, II02s, I0s, Is;

	Cell	**cell;

// p initialization
	po = (double *)malloc(Nx*Ny*sizeof(double));
	a  = (double *)malloc(Nx*Ny*sizeof(double));
	b  = (double *)malloc(Nx*Ny*sizeof(double));

	if ((data32 = (float*)malloc(sizeof(float)*Nx*Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		exit(1);
	}

	ilp=(int *)malloc(NST*sizeof(int));
	iplc=0;

// counting for number of projections (NST = nshot)
	nshot = 0;

// loop between I0_1st and I0_2nd
	for ( j = 0; j < NI0-1; ++j){
//		printf("j %d\n",j);

//IIO[j] and IIO[j+1] are opened
		if(iFlag==0) sprintf(fname, "%s/q%03d.tif", dirin, II0[j]);
		if(iFlag==1) sprintf(fname, "%s/q%04d.tif", dirin, II0[j]);
		if(iFlag==2) sprintf(fname, "%s/q%05d.tif", dirin, II0[j]);
		ReadImageFile(fname,&Nx,&Ny,&BPS,&cell,&desc);
		for(jy=0;jy<Ny;++jy){
			for(jx=0;jx<Nx;++jx){
				*(II01+Nx*jy+jx)=cell[jy][jx]-*(dark+Nx*jy+jx);
			}
		}
		for (jy=0; jy<Ny ; jy++) free(cell[jy]);
		free(cell);
		
		if(iFlag==0) sprintf(fname, "%s/q%03d.tif", dirin, II0[j+1]);
		if(iFlag==1) sprintf(fname, "%s/q%04d.tif", dirin, II0[j+1]);
		if(iFlag==2) sprintf(fname, "%s/q%05d.tif", dirin, II0[j+1]);
		ReadImageFile(fname,&Nx,&Ny,&BPS,&cell,&desc);
		for(jy=0;jy<Ny;++jy){
			for(jx=0;jx<Nx;++jx){
				*(II02+Nx*jy+jx)=cell[jy][jx]-*(dark+Nx*jy+jx);
			}
		}
		for (jy=0; jy<Ny ; jy++) free(cell[jy]);
		free(cell);
		

// I0EV
		t1 = shottime[II0[j]];
		t2 = shottime[II0[j+1]];

// 1 layer(ln)
		for (jx=0;jx<Nx*Ny;++jx){
			I01     = (double)(*(II01+jx));
			I02     = (double)(*(II02+jx));
			*(a+jx) = (double)(((double)(I02   - I01))    / (t2 - t1));
			*(b+jx) = (double)(((double)I01*t2 - (double)I02*t1) / (t2 - t1));
		}
		for ( k = II0[j] + 1; k < II0[j+1]; ++k){
//			fprintf(stderr, "%d\r",k);
			// obtain p(x) from a[jx], b[jx] using shottime[k]
			if(iFlag==0) sprintf(fname, "%s/q%03d.tif", dirin, k);
			if(iFlag==1) sprintf(fname, "%s/q%04d.tif", dirin, k);
			if(iFlag==2) sprintf(fname, "%s/q%05d.tif", dirin, k);
			ReadImageFile(fname,&Nx,&Ny,&BPS,&cell,&desc);
//			printf("k %d %s\n",k,fname);
			Is =0.0;
			I0s=0.0;
			for(jy=0;jy<Ny;++jy){
				for(jx=0;jx<Nx;++jx){
					*(I+Nx*jy+jx)=cell[jy][jx]-*(dark+Nx*jy+jx);
					Is=Is+*(I+Nx*jy+jx);
					*(I0+Nx*jy+jx)=*(a+Nx*jy+jx) * shottime[k] + *(b+Nx*jy+jx);
					I0s=I0s+*(I0+Nx*jy+jx);
				}
			}
//			printf("ff\n");
			for (jy=0; jy<Ny ; jy++) free(cell[jy]);
			free(cell);
			if(cFlag==2) I0cor=Is/I0s;   // correct but auto
//			if(cFlag==1) I0cor=Is/I0s;   //
			if(cFlag==0) I0cor=1.0;      // for no correction
			
//			printf("kk %d\n",k);
			*(ilp+nshot)=0;
			for(jx=0;jx<Nx*Ny;++jx){
				if (*(I+jx) > 10){
//					if(*(ilp+nshot)==0){
//						printf("Warning \t");
//						printf("  jx = %d, I0 = %f, I = %d, dark = %d, ln =%d \n", jx, I0[jx], I[jx], dark[jx], ln);
//						printf("  II01 = %d, II02 = %d\n", II01[jx], II02[jx]);
//						printf("  t1   = %f, t2   = %f \n", t1, t2);
//						printf("  a = %f,  b = %f\n", a[jx], b[jx]);
//						printf("  %d\t black\n", k);
//						*(ilp+nshot)=1;
//					}
					*(po+jx)=log(I0cor*(double)*(I0+jx)/(double)*(I+jx));
				}else{
					*(po+jx)=0.0;
				}
			}

//			if(*(ilp+nshot)==1){
//				for(jx=0;jx<Nx*Ny;++jx){
//					*(po+jx)=0.0;
//				}
//			}
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
			printf("%s\t%lf\r",path, I0cor);
			Store32TiffFile(path, Nx, Ny, 32, data32, comm);

			free(comm);
			
		} // end of k loop
	} // end of j loop
	printf("\n");
	free(a);free(b);free(po);free(data32);
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
	int			i, jx, jy;
//	FILE		*fo;
	char		darkfile[100];

		Cell	**cell;

// parameter setting
	cFlag=0; // no correction
	if (argc<3){
//		fprintf(stderr, "parameter was wrong!!!\n");
		fprintf(stderr, "usage : %s raw/ prj/ (I0cor, 0 for auto)\n", argv[0]);
//		fprintf(stderr, "default head=q, dark=dark.img\n");
		return(1);
	}
	
	if (argc==4){
		if(atof(argv[3])!=0.0){
			I0cor=atof(argv[3]);
			cFlag=1; // correct
		} else{
			cFlag=2; // correct but auto
		}
	}

// read shot log file
	if ((i = read_log(argv[1])) != 0){
		fprintf(stderr, "canot read 'output.log' ! (%d)\n", i);
		return(1);
	}
//	printf("%s/dark.tif\n",argv[1]);

// read dark image
	sprintf(darkfile,"%s/dark.tif",argv[1]);
	ReadImageFile(darkfile,&Nx,&Ny,&BPS,&cell,&desc);
	fprintf(stderr, "Nx, Ny= %d %d \n", Nx,Ny);
	dark = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
	for(jy=0;jy<Ny;++jy){
		for(jx=0;jx<Nx;++jx){
			*(dark+Nx*jy+jx)=cell[jy][jx];
		}
	}
	for (jy=0; jy<Ny ; jy++) free(cell[jy]);
	free(cell);

	II01 = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
	II02 = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
	I    = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
	I0   = (double *) malloc(Nx*Ny*sizeof(double));
//	fprintf(stderr, "Nx, Ny= %d %d \n", Nx,Ny);

	if((i=StoreProjection(argv[1], argv[2])) !=0){
		printf("something wrong in StoreProjection (%d)\n",i);
		return(1);
	}
	fprintf(stderr, "fin \n");

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
