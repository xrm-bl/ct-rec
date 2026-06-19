
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rhp.h"
#include "cbp.h"
//#include "sif_f.h"
//#include "cell.h"
//#include "sif.h"
#include "tiffio.h"
#ifdef USE_GPU
  #include "sort_filter_g.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_gpu
#else
  #include "sort_filter_omp.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_omp
#endif

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
#include <unistd.h>

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


/* ==========================================================================
 *  メモリ量取得(プラットフォーム別) ― hp_tg と同じチャンク判定に使用
 * ========================================================================== */
static unsigned long long get_available_memory_bytes(void)
{
#ifdef _WIN32
	MEMORYSTATUSEX st; st.dwLength=sizeof(st);
	if (GlobalMemoryStatusEx(&st)) return (unsigned long long)st.ullAvailPhys;
	return 0ull;
#else
	long pages = sysconf(_SC_AVPHYS_PAGES);
	long psize = sysconf(_SC_PAGESIZE);
	if (pages>0 && psize>0) return (unsigned long long)pages*(unsigned long long)psize;
	return 0ull;
#endif
}
static unsigned long long get_total_memory_bytes(void)
{
#ifdef _WIN32
	MEMORYSTATUSEX st; st.dwLength=sizeof(st);
	if (GlobalMemoryStatusEx(&st)) return (unsigned long long)st.ullTotalPhys;
	return 0ull;
#else
	long pages = sysconf(_SC_PHYS_PAGES);
	long psize = sysconf(_SC_PAGESIZE);
	if (pages>0 && psize>0) return (unsigned long long)pages*(unsigned long long)psize;
	return 0ull;
#endif
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

	/* ---- chunk 関連 ---- */
	int		*sel, nsel, slice_done;
	int		rows_per_chunk, maxrows, nbands;

	if ((target=(char *)malloc(sizeof(char)*L))==NULL)
	    Error("memory allocation error for range list.");

	RangeList(argv[4],(size_t)L-1,target);

	/* 選択スライスの実 z を sel[] に集約(コンパクト index = SG の行) */
	if ((sel=(int *)malloc(sizeof(int)*L))==NULL)
	    Error("memory allocation error for slice list.");
	nsel=0;
	for (z=0; z<L; z++) if (target[z]) sel[nsel++]=z;

	if (nsel==0) Error("no slice.");

	Dr=atof(argv[5])/10000.;
	DO=(double)(1-N)/2.0;
	RA=atof(argv[6]);

	P=InitCBP(N,M);

	if ((fom=(FOM **)malloc(sizeof(FOM *)*N))==NULL ||
	    (*fom=(FOM *)malloc(sizeof(FOM)*(size_t)N*N))==NULL)
	    Error("memory allocation error for tomogram.");
	for (y=1; y<N; y++) fom[y]=fom[y-1]+N;

	/* ============================================================
	 *  メモリ上限チェック -> 1 バンドあたりのスライス数を決定
	 *  巨大配列 SG = (スライス数) x M x N x sizeof(Float)
	 * ============================================================ */
	{
		unsigned long long mem_avail = get_available_memory_bytes();
		unsigned long long mem_total = get_total_memory_bytes();
		unsigned long long mem_ref   = mem_avail ? mem_avail : mem_total;

		double frac=0.8;
		{ const char *e=getenv("HPTG_MEM_FRACTION"); if (e){ double v=atof(e); if (v>0.05 && v<=0.95) frac=v; } }

		unsigned long long bytes_per_slice = (unsigned long long)M*(unsigned long long)N*sizeof(Float);

		/* SG 以外の固定オーバヘッド見積り */
		unsigned long long overhead =
		      (unsigned long long)M*N*sizeof(Float)          /* P            */
		    + (unsigned long long)N*N*sizeof(FOM)            /* fom          */
		    + 2ull*(unsigned long long)N*M*sizeof(float);    /* ring 作業    */

		unsigned long long budget = (unsigned long long)((double)mem_ref*frac);
		if (budget>overhead) budget -= overhead;
		else                 budget = bytes_per_slice;	/* 最低 1 スライス */

		{ const char *e=getenv("HPTG_MEM_LIMIT_MB");
		  if (e){ double mb=atof(e); if (mb>0) budget=(unsigned long long)(mb*1024.0*1024.0); } }

		{ const char *e=getenv("HPTG_CHUNK_ROWS");
		  if (e && atoi(e)>0) rows_per_chunk=atoi(e);
		  else {
		      unsigned long long rr = bytes_per_slice ? (budget/bytes_per_slice) : (unsigned long long)nsel;
		      if (rr<1) rr=1;
		      rows_per_chunk = (rr>(unsigned long long)nsel)? nsel : (int)rr;
		  }
		}
		if (rows_per_chunk<1) rows_per_chunk=1;
		if (rows_per_chunk>nsel) rows_per_chunk=nsel;

		maxrows = rows_per_chunk;
		nbands  = (nsel + rows_per_chunk - 1) / rows_per_chunk;

		fprintf(stderr,
		    "memory: avail=%.2f GB  total=%.2f GB  SG/slice=%.2f MB  slices/band=%d  bands=%d  (%s)\n",
		    mem_avail/1073741824.0, mem_total/1073741824.0,
		    bytes_per_slice/1048576.0, rows_per_chunk, nbands,
		    (nbands>1)?"multi-pass":"single-pass");
	}

	/* ---- SG(最大バンド分)を確保 ---- */
	if ((SG=(Float ***)malloc(sizeof(Float **)*maxrows))==NULL ||
	    (*SG=(Float **)malloc(sizeof(Float *)*(size_t)maxrows*M))==NULL ||
	    (**SG=(Float *)malloc(sizeof(Float)*(size_t)maxrows*M*N))==NULL)
	    Error("memory allocation error for sinograms.");

	for (m=1; m<M; m++) SG[0][m]=SG[0][m-1]+N;
	for (z=1; z<maxrows; z++) {
	    SG[z]=SG[z-1]+M;
	    for (m=0; m<M; m++) SG[z][m]=SG[z][m-1]+N;
	}

	slice_done=0;

	/* ============================================================
	 *  バンドループ(必要に応じて複数パス)
	 * ============================================================ */
	for (int cs=0; cs<nsel; cs+=rows_per_chunk) {
	    int ce = cs+rows_per_chunk; if (ce>nsel) ce=nsel;
	    int k;

	    /* ---- 投影読み込み: バンド内スライスの SG を構築 ---- */
	    for (m=0; m<M; m++) {
		ReadHiPic(&hp,m);
		for (k=cs; k<ce; k++) {
		    z=sel[k]; sg=SG[k-cs][m]+r0; T=hp.T[z+z0];
		    for (r=0; r<hp.Nx; r++) sg[r]=(-Log(T[r]));
		}

		ReadHiPic(&hp,m+M);
		for (k=cs; k<ce; k++) {
		    z=sel[k]; sg=SG[k-cs][m]; T=hp.T[z+z1]+r1;
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
		fprintf(stderr,"band[%d-%d] read %d / %d\r",cs,ce-1,m+1,M);
	    }
	    fprintf(stderr,"\n");

	    /* ---- バンド内スライスを再構成 ---- */
	    for (k=cs; k<ce; k++) {
		z=sel[k];
		sg=SG[k-cs][0];
		for (m=0; m<M; m++)
		for (r=0; r<N; r++) P[m][r]=(*(sg++));

/* ----------------  black projection correction start ---------------- */
/*                                                                       */
{
			int		blk_m, blk_r, blk_good;
			double	blk_sum;
			int		*blk_flag;
			double	*blk_avg;

			blk_flag = (int *)malloc(M * sizeof(int));
			blk_avg  = (double *)malloc(N * sizeof(double));

			/* initialize average profile */
			for (blk_r = 0; blk_r < N; blk_r++) blk_avg[blk_r] = 0.0;
			blk_good = 0;

			/* detect black projections: all pixels == 0 */
			for (blk_m = 0; blk_m < M; blk_m++){
				blk_sum = 0.0;
				for (blk_r = 0; blk_r < N; blk_r++){
					blk_sum += P[blk_m][blk_r];
				}
				if (blk_sum == 0.0){
					blk_flag[blk_m] = 1;
					(void)fprintf(stderr, "Warning\t black\t m=%d\n", blk_m);
				} else {
					blk_flag[blk_m] = 0;
					blk_good++;
					for (blk_r = 0; blk_r < N; blk_r++){
						blk_avg[blk_r] += P[blk_m][blk_r];
					}
				}
			}

			/* replace black projections with average profile */
			if (blk_good > 0){
				for (blk_r = 0; blk_r < N; blk_r++){
					blk_avg[blk_r] /= (double)blk_good;
				}
				for (blk_m = 0; blk_m < M; blk_m++){
					if (blk_flag[blk_m] == 1){
						for (blk_r = 0; blk_r < N; blk_r++){
							P[blk_m][blk_r] = blk_avg[blk_r];
						}
					}
				}
			}

			free(blk_flag);
			free(blk_avg);
}
/* ----------------  black projection correction finish --------------- */
/*                                                                       */

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
	if (SORT_FILTER_RESTORE(image_data, result_data, N, M, kernel_size, num_threads) != 0) {
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

			if (slice_done++ != 0) TERM_MT(t,Store);

			for (y=0; y<N; y++)
			for (x=0; x<N; x++) fom[y][x]=F[y][x];

#ifdef WINDOWS
			(void)sprintf(path,"%s\\rec%05d.tif",argv[argc-1],Z=z);
#else
	    	(void)sprintf(path,"%s/rec%05d.tif",argv[argc-1],Z=z);
#endif
	    	INIT_MT(t,Store,ac);
	    }
	}

	
	// append to log file
	FILE		*f;
	if((f = fopen("cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(f,"%s ",argv[i]);
	fprintf(f,"   %% kernel_size %d",kernel_size);
	fprintf(f,"   %% bands %d",nbands);
	fprintf(f,"\n");
	fclose(f);


	TERM_MT(t,Store);

	free(*fom); free(fom); free(**SG); free(*SG); free(SG); TermCBP();

	free(sel);
	free(target);
}
	TermReadHiPic(&hp);

	return 0;
}
