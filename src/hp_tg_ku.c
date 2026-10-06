/* ============================================================================
 * hp_tg_ku.c  (HiPic 入力 + メインメモリ上限チェック / チャンク分割版)
 * ----------------------------------------------------------------------------
 * 入力は rhp.h のリーダ(InitReadHiPic / ReadHiPic / TermReadHiPic / HiPic)
 * 経由で読む。出力は従来どおり rec*.tif (32bit float TIFF)。
 *
 * メモリ上限チェックとチャンク(複数パス)実行:
 *   巨大配列 W (= 行数 x Nt x Nx x sizeof(Float)) が空き物理メモリに
 *   収まらない場合、z(スライス=行)方向に分割して複数回に分けて処理する。
 *   収まる場合は従来どおり 1 パス(投影は 1 回だけ読む)。
 *   分割時は各パスで投影を読み直す(I/O は分割数に比例)。
 *   調整用環境変数:
 *       HPTG_MEM_FRACTION : 空きメモリに対する使用率 (既定 0.9)
 *       HPTG_MEM_LIMIT_MB : W に割り当てる上限を MB で直接指定(優先)
 *       HPTG_CHUNK_ROWS   : 1 チャンクあたりの行数を直接指定(最優先)
 *
 * 注意: 本プログラムは Dz==0 固定(スライス z = 行 y)。
 *
 * スライスごとのオーバーヘッド削減(2026-10):
 *   旧版はリング除去用の image_data / result_data と Store() の書き出し
 *   バッファをスライスごとに malloc/free しており、1 スライスあたり約
 *   1.4 GB 分のページフォルトが sys 時間の大半を占めていた。現在は
 *   - リング除去を P(cbp の連続バッファ、GPU 版は pinned)上でインプレース
 *     に実行し、作業バッファと往復コピーを廃止、
 *   - Store() は S.F が連続領域なのでコピーせず直接書き出し、
 *   - W -> P の投影コピー、黒投影検出、F -> S.F の結果コピーは OpenMP 並列。
 *   各要素の演算順序は従来と同一なので結果はビット単位で不変。
 * ==========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rhp.h"
#include "cbp.h"
//#include "sif_f.h"
#include "tiffio.h"
#include "tifwrite.h"
#ifdef USE_GPU
  #include "sort_filter_g.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_gpu
#else
  #include "sort_filter_omp.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_omp
#endif


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
#include <unistd.h>

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


/* ==========================================================================
 *  メモリ量取得(プラットフォーム別)
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
	TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, wdesc);
	TIFFSetField(image, TIFFTAG_ARTIST, "hp_tg");
//	TIFFSetField(image, TIFFTAG_MINSAMPLEVALUE, mmmin );
//	TIFFSetField(image, TIFFTAG_MAXSAMPLEVALUE, mmmax );

		ct_write_raw_strips(image, data32, (uint32_t)wX, (uint32_t)wY, sizeof(float));

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

	/* S.F は N*N の連続領域(S.F[y]=S.F[0]+y*N)なので、FOM==float なら
	   コピーせずそのまま書き出す。FOM!=float のときだけ変換バッファを
	   使う(初回に確保し、以後は使い回す)。 */
	if (sizeof(FOM)==sizeof(float)) {
	    data32=(float *)s->F[0];
	}
	else {
	    static float	*conv=NULL;
	    static size_t	conv_n=0;
	    size_t		need=(size_t)Nx*(size_t)Ny;

	    if (conv_n<need) {
		free(conv);
		if ((conv=(float *)malloc(sizeof(float)*need))==NULL) {
		    printf("cannot allocate memory for output 32bit TIFF image\n");
		    exit(1);
		}
		conv_n=need;
	    }
	    ll=0;
	    for(y=0;y<Ny;++y)
		for(x=0;x<Nx;++x) conv[ll++]=(float)s->F[y][x];
	    data32=conv;
	}

/* insert start */

	mmmin=100.0;
	mmmax=-100.0;
	for(ll=0;ll<(long)Nx*(long)Ny;++ll){
		XXX=*(data32+ll);
		if (mmmin>XXX) mmmin = XXX;
		if (mmmax<XXX) mmmax = XXX;
	}
	if ((comm=(char *)malloc(150))==NULL)
		Error("comment memory allocation error.");

	sprintf(comm,"%f\t%f\t%d\t%f\t%lf\t%lf",Dr*10000.0, RC-Ct, hpNtM, RA0, mmmin, mmmax);

//	StoreImageFile_Float(path,s->N,s->N,s->F,comm);
	Store32TiffFile(path,Nx, Ny, 32, data32, comm);
	(void)printf("%d\t%s\n",s->z,comm);
/* insert end */

	free(comm);
	return RETURN_VALUE;
}

#define EPS	1e-9
#define DEG	3.14159265358979323846/180.0

int	main(int argc,char **argv)
{
	HiPic		hp;
	double		y0,dy,Dz=0.0;
	int		Nz,z1,z2,y,z,t,x,i,j;
	int		Nx,Ny,Nt;
	Float		**P,***W,*w,*P_F,**F;
	Struct		S;
	FOM		*hp_T,fom,*S_F;
	THREAD_T	T;

    int kernel_size = 5; // Default kernel size
    int num_threads = 40; // Default number of threads

	if (argc!=6 && argc!=9){
		fprintf(stderr,"usage : hp_tg HiPic/ Dr RC RA0 rec/\nusage : hp_tg HiPic/ Dr L1 C1 L2 C2 RA0 rec/");
	    Error(" ");  
	}

	InitReadHiPic(argv[1],&hp);
	Nx=hp.Nx; Ny=hp.Ny; Nt=hp.Nt;

	if (Nt<2) Error("too few views.");

	if ((Dr=0.0001*atof(argv[2]))<EPS)   // um -> cm
	    Error("bad horizontal interval of detectors on HiPic image.");

	if (argc==6){
	    Nz=Ny; z1=0; z2=Nz-1; RC=atof(argv[3]); Ct=0.0;
	}
	else{
		z1=atoi(argv[3]);
		z2=atoi(argv[5]);
		if (z1<0) z1=0;
		if (z2>Ny) z2=Ny-1;
		RC=atof(argv[4]);
		Ct = (atof(argv[6])-atof(argv[4]))/(double)(z2-z1);
	}

//	RC=RC-1.0;
	RA0=atof(argv[argc-2])*DEG;

	P=InitCBP(Nx,Nt);

	hpNtM=Nt;

	/* ---- 行範囲 y0 / dy(Dz==0 固定) ---- */
	if (Dz<0.0) { y0=0.5+Dz; }
	else        { y0=0.5; }
	dy=Dz/(double)Nt;

	/* ============================================================
	 *  メモリ上限チェック -> 1 チャンクあたりの行数を決定
	 * ============================================================ */
	unsigned long long mem_avail = get_available_memory_bytes();
	unsigned long long mem_total = get_total_memory_bytes();
	unsigned long long mem_ref   = mem_avail ? mem_avail : mem_total;

	double frac=0.95;
	{ const char *e=getenv("HPTG_MEM_FRACTION"); if (e){ double v=atof(e); if (v>0.05 && v<=0.95) frac=v; } }

	unsigned long long bytes_per_row = (unsigned long long)Nt*(unsigned long long)Nx*sizeof(Float);

	/* W 以外の固定オーバヘッド見積り */
	unsigned long long overhead =
	      (unsigned long long)Nt*Nx*sizeof(Float)            /* P            */
	    + 2ull*(unsigned long long)Nx*Nx*sizeof(Float)       /* S.F + CBP(F) */
	    + 2ull*(unsigned long long)Nx*Nt*sizeof(float);      /* ring 作業    */

	unsigned long long budget = (unsigned long long)((double)mem_ref*frac);
	if (budget>overhead) budget -= overhead;
	else                 budget = bytes_per_row;	/* 最低 1 行 */

	{ const char *e=getenv("HPTG_MEM_LIMIT_MB");
	  if (e){ double mb=atof(e); if (mb>0) budget=(unsigned long long)(mb*1024.0*1024.0); } }

	int total_rows = z2-z1+1;
	int rows_per_chunk;
	{ const char *e=getenv("HPTG_CHUNK_ROWS");
	  if (e && atoi(e)>0) rows_per_chunk=atoi(e);
	  else {
	      unsigned long long r = bytes_per_row ? (budget/bytes_per_row) : (unsigned long long)total_rows;
	      if (r<1) r=1;
	      rows_per_chunk = (r>(unsigned long long)total_rows)? total_rows : (int)r;
	  }
	}
	if (rows_per_chunk<1) rows_per_chunk=1;
	if (rows_per_chunk>total_rows) rows_per_chunk=total_rows;

	int maxrows = rows_per_chunk;
	int nchunks = (total_rows + rows_per_chunk - 1) / rows_per_chunk;

	fprintf(stderr,
	    "memory: avail=%.2f GB  total=%.2f GB  W/row=%.2f MB  rows/chunk=%d  chunks=%d  (%s)\n",
	    mem_avail/1073741824.0, mem_total/1073741824.0,
	    bytes_per_row/1048576.0, rows_per_chunk, nchunks,
	    (nchunks>1)?"multi-pass":"single-pass");

	/* ---- W(最大チャンク分)と S.F を確保 ---- */
	if ((W=(Float ***)malloc(sizeof(Float **)*maxrows))==NULL ||
	    (*W=(Float **)malloc(sizeof(Float *)*(size_t)maxrows*Nt))==NULL ||
	    (**W=(Float *)malloc(sizeof(Float)*(size_t)maxrows*Nt*Nx))==NULL ||
	    (S.F=(FOM **)malloc(sizeof(FOM *)*Nx))==NULL ||
	    (*(S.F)=(FOM *)malloc(sizeof(FOM)*(size_t)Nx*Nx))==NULL)
	    Error("memory allocation error.");

	for (z=0; z<maxrows; z++) {
	    if (z!=0) W[z]=W[z-1]+Nt;
	    for (t=(z==0); t<Nt; t++) W[z][t]=W[z][t-1]+Nx;
	}
	for (y=1; y<Nx; y++) S.F[y]=S.F[y-1]+Nx;

	S.N=Nx;

	/* ---- リング除去は P 上でインプレース実行 ----
	   前提: Float==float、かつ P[t]=P[0]+t*Nx の連続確保
	   (cbp.cu / cbp_thread_*.c とも満たす)。 */
	if (sizeof(Float)!=sizeof(float))
	    Error("in-place ring removal requires Float==float.");
	for (t=1; t<Nt; t++)
	    if (P[t]!=P[t-1]+Nx) Error("projection buffer is not contiguous.");

#ifdef WINDOWS
	(void)sprintf(S.form,"%s\\rec%%0%dd.tif",argv[argc-1],5);
#else
	(void)sprintf(S.form,"%s/rec%%0%dd.tif",argv[argc-1],5);
#endif

	/* ============================================================
	 *  チャンクループ(必要に応じて複数パス)
	 * ============================================================ */
	for (int cs=z1; cs<=z2; cs+=rows_per_chunk) {
	    int ce = cs+rows_per_chunk-1; if (ce>z2) ce=z2;

	    /* このチャンクで必要となる行範囲(Dz==0 では cy1=cs, cy2=ce) */
	    int cy1 = cs + (int)floor(y0 - dy*(double)(Nt-1));
	    int cy2 = ce + (int)floor(y0);
	    if (cy1<0)     cy1=0;
	    if (cy2>Ny-1)  cy2=Ny-1;

	    /* ---- 投影読み込み: 行 [cy1,cy2] を W へ ----
	       Read() はファイル毎に独立、hp の D/I/OL は読み取り専用なので、
	       投影 t を複数スレッドで並列に読む (ReadHiPicBand)。InfiniBand
	       等のディスクアレイは並列アクセスで帯域が出るため、逐次読みで
	       生じていた読み込みフェーズの遊休を短縮する。
	       既定 16 スレッド、環境変数 HPTG_READ_THREADS で変更可。 */
	    {
		int	rthreads=16,rdone=0;
		{ char *e=getenv("HPTG_READ_THREADS");
		  if (e && atoi(e)>0) rthreads=atoi(e); }
		if (rthreads>Nt) rthreads=Nt;

		#pragma omp parallel num_threads(rthreads)
		{
		    unsigned short	*raw;
		    FOM			**dst;
		    Float		*wp,fv;
		    int			tt,yy,xx,rows=cy2-cy1+1;

		    if ((raw=(unsigned short *)malloc(sizeof(unsigned short)
					*(size_t)Ny*(size_t)Nx))==NULL ||
			(dst=(FOM **)malloc(sizeof(FOM *)*(size_t)rows))==NULL)
			Error("memory allocation error for parallel read.");

		    #pragma omp for schedule(dynamic)
		    for (tt=0; tt<Nt; tt++) {
			/* 透過率を W[.][tt] の行へ直接書かせ、同スレッドで
			   -log 変換する (W の tt 列はスレッド毎に別領域) */
			for (yy=0; yy<rows; yy++) dst[yy]=(FOM *)W[yy][tt];

			ReadHiPicBand(&hp,tt,cy1,cy2,dst,raw);

			for (yy=0; yy<rows; yy++) {
			    wp=W[yy][tt];
			    for (xx=0; xx<Nx; xx++) {
				fv=wp[xx];
				wp[xx]=(fv>0.0)?-log(fv):0.0;
			    }
			}
			#pragma omp critical
			fprintf(stderr,"chunk[%d-%d] read %d / %d\r",
				cs,ce,++rdone,Nt);
		    }
		    free(dst); free(raw);
		}
	    }
	    fprintf(stderr,"\n");

	    /* ---- スライス cs..ce を再構成 ---- */
	    for (z=cs; z<=ce; z++) {
		/* 投影 t ごとに独立 -> OpenMP 並列(行単位のコピー) */
		#pragma omp parallel for
		for (t=0; t<Nt; t++) {
		    Float	*pf=P[t],*wf;
		    int		yy=z+(int)floor(y0-dy*(double)t),xx;

		    if (yy<cy1 || yy>cy2)
			for (xx=0; xx<Nx; xx++) pf[xx]=0.0;
		    else {
			wf=W[yy-cy1][t]; for (xx=0; xx<Nx; xx++) pf[xx]=wf[xx];
		    }
		}

/* ----------------  black projection correction start ---------------- */
/*                                                                       */
{
		int		blk_t, blk_r, blk_good=0, blk_black=0;
		int		*blk_flag;
		double	*blk_avg;

		blk_flag = (int *)malloc((size_t)Nt * sizeof(int));
		blk_avg  = (double *)malloc((size_t)Nx * sizeof(double));

		/* detect black projections: all pixels == 0
		   (投影ごとの和は従来と同じ r 昇順の double 加算 -> 同じ判定) */
		#pragma omp parallel for
		for (blk_t = 0; blk_t < Nt; blk_t++){
			double	s=0.0;
			int	r;

			for (r = 0; r < Nx; r++) s += P[blk_t][r];
			blk_flag[blk_t] = (s == 0.0);
		}
		for (blk_t = 0; blk_t < Nt; blk_t++){
			if (blk_flag[blk_t]) {
				(void)fprintf(stderr, "Warning\t black\t z=%d t=%d\n", z, blk_t);
				blk_black++;
			} else blk_good++;
		}

		/* replace black projections with average profile
		   (平均は黒投影があるときだけ計算。列ごとに t 昇順で加算するので
		    従来の加算順と同一 -> 同じ値) */
		if (blk_black > 0 && blk_good > 0){
			#pragma omp parallel for
			for (blk_r = 0; blk_r < Nx; blk_r++){
				double	s=0.0;
				int	tt;

				for (tt = 0; tt < Nt; tt++)
					if (!blk_flag[tt]) s += P[tt][blk_r];
				blk_avg[blk_r] = s / (double)blk_good;
			}
			for (blk_t = 0; blk_t < Nt; blk_t++){
				if (blk_flag[blk_t] == 1){
					for (blk_r = 0; blk_r < Nx; blk_r++){
						P[blk_t][blk_r] = blk_avg[blk_r];
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
	    // Get kernel size from environment variable
	    kernel_size = get_kernel_size_from_env();
		// Get number of threads from environment variable
	    num_threads = get_num_threads_from_env();
	    /* P 上でインプレース実行(入力=出力)。GPU 版は H2D 後に D2H、
	       CPU 版は列ごとに全行を読んでから書くので、どちらも安全。 */
		if (SORT_FILTER_RESTORE(P[0], P[0], Nx, Nt, kernel_size, num_threads) != 0) {
			fprintf(stderr, "ring removal failed\n");
			return 5;
		}
/* ----------------  ring removal finish --------------- */
/*                                                       */

		F=CBP(Dr,-RC,RA0);
		RC=RC+Ct;

		if (z!=z1) TERM_MT(T);

		S.z=z;

		#pragma omp parallel for
		for (y=0; y<Nx; y++) {
		    FOM		*sf=S.F[y];
		    Float	*ff=F[y];
		    int		xx;

		    for (xx=0; xx<Nx; xx++) sf[xx]=ff[xx];
		}
		INIT_MT(T,Store,&S);
	    }
	}

	TermReadHiPic(&hp);

	// append to log file
	FILE		*f;
	if((f = fopen("cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(f,"%s ",argv[i]);
	fprintf(f,"   %% kernel_size %d",kernel_size);
	fprintf(f,"   %% chunks %d",nchunks);
	fprintf(f,"\n");
	fclose(f);


	TERM_MT(T);

#ifdef USE_GPU
	sort_filter_gpu_release();
#endif
	free(*(S.F)); free(S.F); free(**W); free(*W); free(W); TermCBP();

	return 0;
}
