
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rhp.h"
#include "cbp.h"
//#include "sif_f.h"
//#include "cell.h"
//#include "sif.h"
#include "tiffio.h"
#include "sort_filter_omp.h"

extern void	Error(char *msg),
		RangeList(char *rl,size_t limit,char *target);

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
if (INIT_THREAD(T,F,A)) Error(#F " : multi-threading initialization error.")

#define TERM_MT(T,F)	\
	if (TERM_THREAD(T)) Error(#F " : multi-threading termination error.")

#define LEN	2048

static int	N,BPS,Z;
static double	B,S,F1,F2;
//static Cell	C;
static FOM	**fom;
static char	path[LEN];

/* insert start */
	double		Dr,DO,RA;
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
	TIFFSetField(image, TIFFTAG_ARTIST, "ofct_srec");
//	TIFFSetField(image, TIFFTAG_MINSAMPLEVALUE, mmmin );
//	TIFFSetField(image, TIFFTAG_MAXSAMPLEVALUE, mmmax );

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data32 + i*wX, wX * sizeof(float));
	}

	TIFFClose(image);
}

static FUNCTION_T	Store(void *argc)
{
	double	f1,f2;
	int	y,x;
//s	Cell	**c;
	char	desc[LEN];

/* insert start */
	char	*comm = NULL;
	double	mmmin, mmmax, XXX;
	long	ll;
/* insert end */

	float *data32;
	int Nx, Ny;

	Nx=N;
	Ny=N;
	
	if ((data32 = (float*)malloc(sizeof(float)*Nx*Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		exit(1);
	}
	mmmin=100.0;
	mmmax=-100.0;
	ll=0;
	for (y=0; y<Ny; y++){
		for (x=0; x<Nx; x++){
			*(data32+ll)=fom[y][x];
			XXX=*(data32+ll);
			if (mmmin>XXX) mmmin = XXX;
			if (mmmax<XXX) mmmax = XXX;
			ll=ll+1;
		}
	}

	if ((comm=(char *)malloc(150))==NULL)
		Error("comment memory allocation error.");
	sprintf(comm,"%f\t%f\t%d\t%f\t%lf\t%lf",Dr*10000.0, -DO, hpNtM, RA, mmmin, mmmax);

//		StoreImageFile_Float(path,N,N,fom,comm);
		Store32TiffFile(path,Nx, Ny, 32, data32, comm);
	    (void)printf("%d\t%s\n",Z,comm);

	free(data32); free(comm);
	return RETURN_VALUE;
}

static double	Log(double d)
{
	return (d>0.0)?log(d):0.0;
}

int	main(int argc,char **argv)
{
	HiPic	hp;
	int	Ox,Oy,r0,r1,r2,r3,r4,r5,L,z0,z1,M,i;

    int kernel_size = 5; // Default kernel size
    int num_threads = 40; // Default number of threads

	if (argc!=4 && argc!=8)
	    Error("usage : ofct_srec HiPic/ Rc Oy {rangeList Dr RA0 rec/}"
	    );

	InitReadHiPic(argv[1],&hp);

	if (hp.Nt%2) (void)fputs("bad number of views (warning).\n",stderr);

	if ((Ox=2*atoi(argv[2])-hp.Nx)+hp.Nx<=0 || hp.Nx<=Ox ||
	    (Oy=atoi(argv[3]))+hp.Ny<=0 || hp.Ny<=Oy) Error("bad offset.");

	if (Ox<0) {
	    N=hp.Nx-Ox; r0=(-Ox); r1=hp.Nx-1; r2=0; r3=r4=(-Ox); r5=hp.Nx;
	}
	else {
	    N=hp.Nx+Ox; r0=0; r1=N-1; r2=hp.Nx; r3=N; r4=Ox; r5=r2;
	}
	if (Oy<0) {
	    L=hp.Ny+Oy; z0=0; z1=(-Oy);
	}
	else {
	    L=hp.Ny-Oy; z0=Oy; z1=0;
	}
	M=hp.Nt/2;
	hpNtM=M;
	(void)fprintf(stderr,"%d\t%d\t%d\t",N,L,M);
	if (argc==4) (void)fprintf(stderr,"%f GB\n",4*N*L*M/1000000000.);
	(void)fprintf(stderr,"\n");

	if (argc!=4)
{
	char		*target;
	int		l,z,m,y,r,x, i,j;
//	double		Dr,DO,RA;
	Float		**P,***SG,*sg,**F;
	FOM		*T;
	THREAD_T	t;
	size_t		ac=argc;

	if ((target=(char *)malloc(sizeof(char)*L))==NULL)
	    Error("memory allocation error for range list.");

	RangeList(argv[4],(size_t)L-1,target);

	l=0; for (z=0; z<L; z++) if (target[z]) ++l;

	if (l==0) Error("no slice.");

	Dr=atof(argv[5])/10000.;
	DO=(double)(1-N)/2.0;
	RA=atof(argv[6]);

	P=InitCBP(N,M);

	if ((SG=(Float ***)malloc(sizeof(Float **)*l))==NULL ||
	    (*SG=(Float **)malloc(sizeof(Float *)*l*M))==NULL ||
	    (**SG=(Float *)malloc(sizeof(Float)*l*M*N))==NULL ||
	    (fom=(FOM **)malloc(sizeof(FOM *)*N))==NULL ||
	    (*fom=(FOM *)malloc(sizeof(FOM)*N*N))==NULL)
	    Error("memory allocation error for sinograms or tomograms.");

	for (m=1; m<M; m++) SG[0][m]=SG[0][m-1]+N;

	for (z=1; z<l; z++) {
	    SG[z]=SG[z-1]+M;

	    for (m=0; m<M; m++) SG[z][m]=SG[z][m-1]+N;
	}
	for (y=1; y<N; y++) fom[y]=fom[y-1]+N;

	for (m=0; m<M; m++) {
	    ReadHiPic(&hp,m);
	    for (l=z=0; z<L; z++)
		if (target[z]) {
		    sg=SG[l++][m]+r0; T=hp.T[z+z0];
		    for (r=0; r<hp.Nx; r++) sg[r]=(-Log(T[r]));
		}

	    ReadHiPic(&hp,m+M);
	    for (l=z=0; z<L; z++)
		if (target[z]) {
		    sg=SG[l++][m]; T=hp.T[z+z1]+r1;
		    for (r=r2; r<r3; r++) sg[r]=(-Log(T[-r]));
#ifdef	OCT_SBS	/* side by side */
	r=r4+(r5-r4-1)/2; if ((r5-r4)%2) sg[r]=(sg[r]-Log(T[-r]))/2.0;

	while (++r<r5) sg[r]=(-Log(T[-r]));
#else
	for (r=r4; r<r5; r++) sg[r]=
#ifdef	OCT_SA	/* simple average */
	(sg[r]-Log(T[-r]))/2.0;
#else		/* linear mixing */
	((double)(r5-r)*sg[r]+(double)(r-r4+1)*(-Log(T[-r])))/(double)(r5-r4+1);
#endif
#endif
		}
	}
	for (l=z=0; z<L; z++)
	    if (target[z]) {
			sg=SG[l][0];
			for (m=0; m<M; m++)
			for (r=0; r<N; r++) P[m][r]=(*(sg++));

/* ----------------  ring removal start ---------------- */
/*                                                       */
    float		*image_data = NULL, *result_data = NULL;
    // Get kernel size from environment variable
    kernel_size = get_kernel_size_from_env();
	// Get number of threads from environment variable
    num_threads = get_num_threads_from_env();
    // Allocate memory
	image_data = (float *)malloc(N * M * sizeof(float));
	result_data = (float *)malloc(N * M * sizeof(float));

	for (j=0; j<M; j++){
		for (i=0; i<N; i++){
			*(image_data+N*j+i)=P[j][i];
		}
	}
	// Execute OpenMP image processing
	if (sort_filter_restore_omp(image_data, result_data, N, M, kernel_size, num_threads) != 0) {
		fprintf(stderr, "OpenMP image processing failed\n");
		return 5;
	}
	for (j=0; j<M; j++){
		for (i=0; i<N; i++){
				P[j][i]=*(result_data+N*j+i);
		}
	}
    if (image_data) free(image_data);
    if (result_data) free(result_data);
/* ----------------  ring removal finish --------------- */
/*                                                       */
			F=CBP(Dr,DO,RA);

			if (++l!=1) TERM_MT(t,Store);

			for (y=0; y<N; y++)
			for (x=0; x<N; x++) fom[y][x]=F[y][x];
	    	
#ifdef WINDOWS
			(void)sprintf(path,"%s\\rec%05d.tif",argv[argc-1],Z=z);
#else
	    	(void)sprintf(path,"%s/rec%05d.tif",argv[argc-1],Z=z);
#endif
	    	INIT_MT(t,Store,ac);
		}

	
	// append to log file
	FILE		*f;
	if((f = fopen("cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(f,"%s ",argv[i]);
	fprintf(f,"   %% kernel_size %d",kernel_size);
	fprintf(f,"\n");
	fclose(f);


	TERM_MT(t,Store);

	free(*fom); free(fom); free(**SG); free(*SG); free(SG); TermCBP();

	free(target);
}
	TermReadHiPic(&hp);

	return 0;
}
