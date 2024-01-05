
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rhp.h"
#include "cbp.h"
//#include "sif_f.h"
#include "tiffio.h"


extern void	Error(char *msg);

#ifdef	_WIN32
#include <process.h>
#include <windows.h>

#define THREAD_T	HANDLE
#define FUNCTION_T	unsigned __stdcall
#define RETURN_VALUE	0

#define INIT_THREAD(T,F,A)	\
	(T=(HANDLE)_beginthreadex(NULL,0,F,(void *)(A),0,NULL))==0
#define TERM_THREAD(T)	\
	WaitForSingleObject(T,INFINITE)==WAIT_FAILED || CloseHandle(T)==0
#else
#include <pthread.h>

#define THREAD_T	pthread_t
#define FUNCTION_T	void *
#define RETURN_VALUE	NULL

#define INIT_THREAD(T,F,A)	pthread_create(&(T),NULL,F,(void *)(A))
#define TERM_THREAD(T)		pthread_join(T,NULL)
#endif

#define INIT_MT(T,F,A)	\
if (INIT_THREAD(T,F,A)) Error("multi-threading initialization error.")
#define TERM_MT(T)	\
if (TERM_THREAD(T)) Error("multi-threading termination error.")

#define LEN	2048

typedef struct {
		char	form[LEN];
		int	z,N;
		FOM	**F;
	} Struct;

/* insert start */
	double		Dr,RC,RA0,Ct;
	int			hpNtM;
/* insert end */

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
	TIFFSetField(image, TIFFTAG_ARTIST, "hp_tg");
//	TIFFSetField(image, TIFFTAG_MINSAMPLEVALUE, mmmin );
//	TIFFSetField(image, TIFFTAG_MAXSAMPLEVALUE, mmmax );

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data32 + i*wX, wX * sizeof(float));
	}

	TIFFClose(image);
}

static FUNCTION_T	Store(void *a)
{
	Struct	*s=(Struct *)a;
	char	path[LEN];

/* insert start */
	char	*comm = NULL;
	double	mmmin, mmmax, XXX;
	int		x,y;
	long	ll;
/* insert end */

	float *data32;
	int Nx, Ny;

	(void)sprintf(path,s->form,s->z);
	Nx=s->N;
	Ny=s->N;

	if ((data32 = (float*)malloc(sizeof(float)*Nx*Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		exit(1);
	}

/* insert start */

	mmmin=100.0;
	mmmax=-100.0;
	ll=0;
	for(y=0;y<Ny;++y){
		for(x=0;x<Nx;++x){
			*(data32+ll)=s->F[y][x];
			XXX=*(data32+ll);
			if (mmmin>XXX) mmmin = XXX;
			if (mmmax<XXX) mmmax = XXX;
			ll=ll+1;
//			printf("%d\t%d\t%lf\n",x,y,XXX);
		}
	}
	if ((comm=(char *)malloc(150))==NULL)
		Error("comment memory allocation error.");

	sprintf(comm,"%f\t%f\t%d\t%f\t%lf\t%lf",Dr*10000.0, RC-Ct, hpNtM, RA0, mmmin, mmmax);

//	StoreImageFile_Float(path,s->N,s->N,s->F,comm);
	Store32TiffFile(path,Nx, Ny, 32, data32, comm);
	(void)printf("%d\t%s\n",s->z,comm);
/* insert end */

	free(data32); free(comm);
	return RETURN_VALUE;
}

#define EPS	1e-9
#define DEG	3.14159265358979323846/180.0

int	main(int argc,char **argv)
{
	HiPic		hp;
	double		y0,dy,Dz=0.0, cc;
	int		Nz,z1,z2,y1,y2,y,z,t,x,len,i;
	Float		**P,***W,*w,*P_F,**F;
	Struct		S;
	FOM		*hp_T,fom,*S_F;
	THREAD_T	T;

	if (argc!=6 && argc!=9){
		fprintf(stderr,"usage : hp_tg HiPic/ Dr RC RA0 rec/\nusage : hp_tg HiPic/ Dr L1 C1 L2 C2 RA0 rec/");
	    Error(" ");  
	}

	InitReadHiPic(argv[1],&hp);

	if (hp.Nt<2) Error("too few views.");

	if ((Dr=0.0001*atof(argv[2]))<EPS)   // um -> cm
	    Error("bad horizontal interval of detectors on HiPic image.");

	if (argc==6){
	    Nz=hp.Ny; z1=0; z2=Nz-1; RC=atof(argv[3]); Ct=0.0;
	}
	else{
		z1=atoi(argv[3]);
		z2=atoi(argv[5]);
		if (z1<0) z1=0;
		if (z2>hp.Ny) z2=hp.Ny-1;
		RC=atof(argv[4]);
		Ct = (atof(argv[6])-atof(argv[4]))/(double)(z2-z1);
	}

//	RC=RC-1.0;
	RA0=atof(argv[argc-2])*DEG;

	P=InitCBP(hp.Nx,hp.Nt);

	hpNtM=hp.Nt;
	
	if (Dz<0.0) {
	    y0=0.5+Dz; if ((y1=z1+(int)floor(y0))<0) y1=0;
	}
	else {
	    y0=0.5; if ((y1=z1+(int)floor(y0-Dz))<0) y1=0;
	}
	y=(y2=(z2<hp.Ny)?z2:hp.Ny-1)-y1+1;
	if ((W=(Float ***)malloc(sizeof(Float **)*y))==NULL ||
	    (*W=(Float **)malloc(sizeof(Float *)*y*hp.Nt))==NULL ||
	    (**W=(Float *)malloc(sizeof(Float)*y*hp.Nt*hp.Nx))==NULL ||
	    (S.F=(FOM **)malloc(sizeof(FOM *)*hp.Nx))==NULL ||
	    (*(S.F)=(FOM *)malloc(sizeof(FOM)*hp.Nx*hp.Nx))==NULL)
	    Error("memory allocation error.");

	for (z=0; z<y; z++) {
	    if (z!=0) W[z]=W[z-1]+hp.Nt;

	    for (t=z==0; t<hp.Nt; t++) W[z][t]=W[z][t-1]+hp.Nx;
	}
	for (y=1; y<hp.Nx; y++) S.F[y]=S.F[y-1]+hp.Nx;

	for (t=0; t<hp.Nt; t++) {
	    ReadHiPic(&hp,t); fprintf(stderr, "%d / %d\r",t+1,hp.Nt);

	    for (y=y1; y<=y2; y++) {
		w=W[y-y1][t]; hp_T=hp.T[y];
		for (x=0; x<hp.Nx; x++)
		    *(w++)=((fom=(*(hp_T++)))>0.0)?-log(fom):0.0;
	    }
	}
	TermReadHiPic(&hp); printf("\n");

	S.N=hp.Nx;

//	len=1; for (z=10; z<Nz; z*=10) ++len;
//	(void)sprintf(S.form,"%s/rec%05d.tif",argv[argc-1],len);

#ifdef WINDOWS
	(void)sprintf(S.form,"%s\\rec%%0%dd.tif",argv[argc-1],5);
#else
	(void)sprintf(S.form,"%s/rec%%0%dd.tif",argv[argc-1],5);
#endif

	dy=Dz/(double)hp.Nt; 
	for (z=z1; z<=z2; z++) {
	    for (t=0; t<hp.Nt; t++) {
		P_F=P[t];
		if ((y=z+(int)floor(y0-dy*(double)t))<y1 || y>y2)
		    for (x=0; x<hp.Nx; x++) *(P_F++)=0.0;
		else {
		    w=W[y-y1][t]; for (x=0; x<hp.Nx; x++) *(P_F++)=(*(w++));
		}
	    }
	    F=CBP(Dr,-RC,RA0);
		RC=RC+Ct;

	    if (z!=z1) TERM_MT(T);

	    S.z=z;

	    for (y=0; y<hp.Nx; y++) {
		S_F=S.F[y]; P_F=F[y];
		for (x=0; x<hp.Nx; x++) *(S_F++)=(*(P_F++));
	    }
	    INIT_MT(T,Store,&S);
	}
	TERM_MT(T);

	free(*(S.F)); free(S.F); free(**W); free(*W); free(W); TermCBP();

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
