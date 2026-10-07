/* ============================================================================
 * ofact_srec.c  自動CT装置用のオフセットCT連続再構成 (ofct_srec.c から派生、2026-10)
 * ----------------------------------------------------------------------------
 *   ofact_srec HiPic Rc Oy                                   (サイズ見積もりのみ)
 *   ofact_srec HiPic Rc Oy rangeList Dr RA0 rec16 LACmin16 LACmax16 rec8 LACmin8 LACmax8
 *
 * 出力は 32bit float TIFF ではなく、規格化範囲を引数で与えて
 *   16bit  rec16/rh%05d.tif  (LACmin16 .. LACmax16)
 *    8bit  rec8/ro%05d.tif   (LACmin8  .. LACmax8)
 * を同時に直接書き出す。再構成像(double)は一度 float に丸めてから量子化する
 * ので、量子化規則・8 欄の ImageDescription とも tif_f2i と同一で、
 * 「ofct_srec (32bit) -> tif_f2i」の結果と画素値・記述ともに一致する。
 * 再構成処理(読み込み、合成、リング除去、CBP、パイプライン)は ofct_srec と同じ。
 * ==========================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rhp.h"
#include "cbp.h"
#include "cbp_pipe.h"
//#include "sif_f.h"
//#include "cell.h"
//#include "sif.h"
#include "tiffio.h"
#include "tifwrite.h"
#ifdef USE_GPU
  #include "sort_filter_g.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_gpu
#else
  #include "sort_filter_omp.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_omp
#endif

extern void	Error(char *msg),
		RangeList(char *rl,size_t limit,char *target);

/* スライスごとのオーバーヘッド削減(2026-10、hp_tg_ku.c と同じ変更):
   旧版はリング除去用の image_data / result_data と Store() の float 変換
   バッファをスライスごとに malloc/free しており、1 スライスあたり数 GB 分の
   ページフォルトが発生していた。現在は、リング除去を P(cbp の連続バッファ、
   GPU 版は pinned)上でインプレースに実行して作業バッファと往復コピーを廃止、
   Store() の変換バッファは初回のみ確保して使い回し、SG -> P の投影コピー、
   黒投影検出、F -> fom の結果コピーは OpenMP 並列。各要素の演算順序は従来と
   同一なので結果はビット単位で不変。
   さらに投影バッファを 2 面持ち、スライス z の BeginCBP の直後に z+1 の
   準備(SG -> P、黒投影検出、リング除去)を行う 1 スライス先行パイプライン
   (cbp_pipe.h、hp_tg_ku.c と同じ構造)。 */

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
static char	path16[LEN],path8[LEN];	/* 16bit(rh) / 8bit(ro) の出力パス */
static double	LAC16lo,LAC16hi,LAC8lo,LAC8hi;	/* 規格化範囲 */
static unsigned short	*out16=NULL;		/* 量子化バッファ(1 回確保) */
static unsigned char	*out8=NULL;
static int	cN;				/* 偶数に切り詰めた出力サイズ(tif_f2i と同じ) */

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

static void Store8TiffFile(char *wname, int wX, int wY, unsigned char *data8, char *wdesc)
{
	TIFF *image;

	if ((image = TIFFOpen(wname, "w"))==NULL) Error("cannot open 8bit output file.");

	TIFFSetField(image, TIFFTAG_IMAGEWIDTH, wX);
	TIFFSetField(image, TIFFTAG_IMAGELENGTH, wY);
	TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, wdesc);
	TIFFSetField(image, TIFFTAG_ARTIST, "ofact_srec");

	ct_write_raw_strips(image, data8, (uint32_t)wX, (uint32_t)wY, sizeof(unsigned char));

	TIFFClose(image);
}

static void Store16TiffFile(char *wname, int wX, int wY, unsigned short *data16, char *wdesc)
{
	TIFF *image;

	if ((image = TIFFOpen(wname, "w"))==NULL) Error("cannot open 16bit output file.");

	TIFFSetField(image, TIFFTAG_IMAGEWIDTH, wX);
	TIFFSetField(image, TIFFTAG_IMAGELENGTH, wY);
	TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 16);
	TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, wdesc);
	TIFFSetField(image, TIFFTAG_ARTIST, "ofact_srec");

	ct_write_raw_strips(image, data16, (uint32_t)wX, (uint32_t)wY, sizeof(unsigned short));

	TIFFClose(image);
}

static FUNCTION_T	Store(void *argc)
{
	char	comm[256];
	double	mmmin=100.0,mmmax=-100.0,XXX;
	int	y,x;
	double	div16=65535.0/(LAC16hi-LAC16lo),
		div8 =255.0  /(LAC8hi -LAC8lo );

	/* fom(double) を float に丸めた値で最小・最大と量子化を行う
	   (32bit TIFF を経由して tif_f2i にかけた場合と同じ値になる) */
	for (y=0; y<N; y++)
	    for (x=0; x<N; x++) {
		XXX=(double)(float)fom[y][x];
		if (mmmin>XXX) mmmin=XXX;
		if (mmmax<XXX) mmmax=XXX;
	    }

	#pragma omp parallel for private(x)
	for (y=0; y<cN; y++) {
	    for (x=0; x<cN; x++) {
		double	v=(double)(float)fom[y][x];
		long	val;

		val=(long)(div16*(v-LAC16lo));
		if (val<1) val=0;
		if (val>65535) val=65535;
		out16[(size_t)y*(size_t)cN+x]=(unsigned short)val;

		val=(long)(div8*(v-LAC8lo));
		if (val<1) val=0;
		if (val>255) val=255;
		out8[(size_t)y*(size_t)cN+x]=(unsigned char)val;
	    }
	}

	/* 8 欄の ImageDescription: ofct_srec の 6 欄 + tif_f2i が付ける規格化範囲 2 欄 */
	sprintf(comm,"%f\t%f\t%d\t%f\t%lf\t%lf\t%lf\t%lf",
		Dr*10000.0, -DO, hpNtM, RA, mmmin, mmmax, LAC16lo, LAC16hi);
	Store16TiffFile(path16,cN,cN,out16,comm);

	sprintf(comm,"%f\t%f\t%d\t%f\t%lf\t%lf\t%lf\t%lf",
		Dr*10000.0, -DO, hpNtM, RA, mmmin, mmmax, LAC8lo, LAC8hi);
	Store8TiffFile(path8,cN,cN,out8,comm);

	(void)printf("%d\t%f\t%f\t%d\t%f\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
		Z, Dr*10000.0, -DO, hpNtM, RA, mmmin, mmmax,
		LAC16lo, LAC16hi, LAC8lo, LAC8hi);
	return RETURN_VALUE;
}

static double	Log(double d)
{
	return (d>0.0)?log(d):0.0;
}

/* ---- 1 スライス分のホスト側準備: SG[k] -> Pb、黒投影検出、リング除去 ----
   GPU 版では直前のスライスの BeginCBP の後に呼ばれ、GPU 計算と並行して走る。 */
static int	PrepareSlice(Float **Pb,Float **sgk,int N,int M,
			     int *kernel_size,int *num_threads)
{
	int	m;

	/* 投影 m ごとに独立 -> OpenMP 並列(行単位のコピー) */
	#pragma omp parallel for
	for (m=0; m<M; m++) {
	    Float	*pf=Pb[m],*sf=sgk[m];
	    int		rr;

	    for (rr=0; rr<N; rr++) pf[rr]=sf[rr];
	}

/* ----------------  black projection correction start ---------------- */
	{
		int		blk_m, blk_r, blk_good=0, blk_black=0;
		int		*blk_flag;
		double	*blk_avg;

		blk_flag = (int *)malloc(M * sizeof(int));
		blk_avg  = (double *)malloc(N * sizeof(double));
		if (blk_flag==NULL || blk_avg==NULL)
			Error("memory allocation error for black projection check.");

		/* detect black projections: all pixels == 0
		   (投影ごとの和は従来と同じ r 昇順の double 加算 -> 同じ判定) */
		#pragma omp parallel for
		for (blk_m = 0; blk_m < M; blk_m++){
			double	s=0.0;
			int	r;

			for (r = 0; r < N; r++) s += Pb[blk_m][r];
			blk_flag[blk_m] = (s == 0.0);
		}
		for (blk_m = 0; blk_m < M; blk_m++){
			if (blk_flag[blk_m]) {
				(void)fprintf(stderr, "Warning\t black\t m=%d\n", blk_m);
				blk_black++;
			} else blk_good++;
		}

		/* replace black projections with average profile
		   (黒投影があるときだけ計算。列ごとに m 昇順で加算 = 従来と同じ順) */
		if (blk_black > 0 && blk_good > 0){
			#pragma omp parallel for
			for (blk_r = 0; blk_r < N; blk_r++){
				double	s=0.0;
				int	mm;

				for (mm = 0; mm < M; mm++)
					if (!blk_flag[mm]) s += Pb[mm][blk_r];
				blk_avg[blk_r] = s / (double)blk_good;
			}
			for (blk_m = 0; blk_m < M; blk_m++){
				if (blk_flag[blk_m] == 1){
					for (blk_r = 0; blk_r < N; blk_r++){
						Pb[blk_m][blk_r] = blk_avg[blk_r];
					}
				}
			}
		}

		free(blk_flag);
		free(blk_avg);
	}
/* ----------------  black projection correction finish --------------- */

/* ----------------  ring removal (in place on Pb) ---------------- */
	*kernel_size = get_kernel_size_from_env();
	*num_threads = get_num_threads_from_env();
	if (SORT_FILTER_RESTORE(Pb[0], Pb[0], N, M, *kernel_size, *num_threads) != 0) {
		fprintf(stderr, "ring removal failed\n");
		return 5;
	}
	return 0;
}

int	main(int argc,char **argv)
{
	HiPic	hp;
	int	Ox,Oy,r0,r1,r2,r3,r4,r5,L,z0,z1,M,i;

    int kernel_size = 5; // Default kernel size
    int num_threads = 40; // Default number of threads

	if (argc!=4 && argc!=13)
	    Error("usage : ofact_srec HiPic/ Rc Oy {rangeList Dr RA0 rec16/ LACmin16 LACmax16 rec8/ LACmin8 LACmax8}"
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
	int		l,z,m,y,r,x, i,j, cur, pipe_on=1;
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
	DO=-(double)(N+1)/2.0;	/* rotation center = -(N+1)/2, matching ofct_rec's r0=-(dcnt+0.5); was -(N-1)/2 (1 px off) */
	RA=atof(argv[6]);

	/* ---- 8/16bit 直接出力の規格化範囲とバッファ ---- */
	LAC16lo=atof(argv[8]);  LAC16hi=atof(argv[9]);
	LAC8lo =atof(argv[11]); LAC8hi =atof(argv[12]);
	if (!(LAC16hi>LAC16lo)) Error("LACmin16 must be smaller than LACmax16.");
	if (!(LAC8hi >LAC8lo )) Error("LACmin8 must be smaller than LACmax8.");
	cN=(N%2)?N-1:N;			/* tif_f2i と同じ偶数切り詰め */
	if ((out16=(unsigned short *)malloc(sizeof(unsigned short)*(size_t)cN*(size_t)cN))==NULL ||
	    (out8 =(unsigned char  *)malloc(sizeof(unsigned char )*(size_t)cN*(size_t)cN))==NULL)
	    Error("memory allocation error for 8/16bit output.");

	P=InitCBP(N,M);
	CBP_PIPE_INIT(P,N,M);

	if ((fom=(FOM **)malloc(sizeof(FOM *)*N))==NULL ||
	    (*fom=(FOM *)malloc(sizeof(FOM)*(size_t)N*N))==NULL)
	    Error("memory allocation error for tomogram.");
	for (y=1; y<N; y++) fom[y]=fom[y-1]+N;

	/* ---- 投影バッファ 2 面(1 スライス先行パイプライン) ----
	   Pbuf[0] は cbp の P、Pbuf[1] は AllocCBPProjection() の 2 面目
	   (GPU 版は pinned)。リング除去はこの上でインプレース実行するので、
	   Float==float と P[m]=P[0]+m*N の連続確保が前提
	   (cbp.cu / cbp_thread_*.c とも満たす)。 */
	Float	**Pbuf[2];

	Pbuf[0]=P;
	if ((Pbuf[1]=CBP_PIPE_ALLOC())==NULL)
	    Error("memory allocation error for the second projection buffer.");
	/* HPTG_PIPELINE=0 で先行パイプラインを止め、同じ関数を逐次順
	   (準備 -> Begin -> End -> 保存)で呼ぶ。切り分け・比較用。 */
	{ char *e=getenv("HPTG_PIPELINE"); if (e && atoi(e)==0) pipe_on=0; }
	if (!pipe_on) fprintf(stderr,"HPTG_PIPELINE=0: sequential order\n");
	if (sizeof(Float)!=sizeof(float))
	    Error("in-place ring removal requires Float==float.");
	for (i=0; i<2; i++)
	    for (m=1; m<M; m++)
		if (Pbuf[i][m]!=Pbuf[i][m-1]+N)
		    Error("projection buffer is not contiguous.");

	/* ============================================================
	 *  メモリ上限チェック -> 1 バンドあたりのスライス数を決定
	 *  巨大配列 SG = (スライス数) x M x N x sizeof(Float)
	 * ============================================================ */
	{
		unsigned long long mem_avail = get_available_memory_bytes();
		unsigned long long mem_total = get_total_memory_bytes();
		unsigned long long mem_ref   = mem_avail ? mem_avail : mem_total;

		double frac=0.9;
		{ const char *e=getenv("HPTG_MEM_FRACTION"); if (e){ double v=atof(e); if (v>0.05 && v<=0.95) frac=v; } }

		unsigned long long bytes_per_slice = (unsigned long long)M*(unsigned long long)N*sizeof(Float);

		/* SG 以外の固定オーバヘッド見積り */
		unsigned long long overhead =
		      (unsigned long long)M*N*sizeof(Float)          /* P            */
		    + (unsigned long long)N*N*sizeof(FOM);           /* fom          */

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

	    /* ---- 投影読み込み: バンド内スライスの SG を構築 ----
	       hp_tg と同様に投影ペア m / m+M を複数スレッドで並列に読む
	       (ReadHiPicBand)。各スレッドの出力は SG[.][m] 列で m 毎に
	       独立なので競合しない。dark/入射補正・-log・オフセット合成の
	       式は従来と同一で、結果はビット一致する。
	       既定 16 スレッド、環境変数 HPTG_READ_THREADS で変更可。 */
	    {
		int	rthreads=16,rdone=0,
			rows=sel[ce-1]-sel[cs]+1;	/* バンドの被覆行数 */
		{ char *e=getenv("HPTG_READ_THREADS");
		  if (e && atoi(e)>0) rthreads=atoi(e); }
		if (rthreads>M) rthreads=M;

		#pragma omp parallel num_threads(rthreads)
		{
		    unsigned short	*raw;
		    FOM			**dst,*T2;
		    Float		*sg2;
		    int			mm,kk,rr,yy,y1a,y1b;

		    if ((raw=(unsigned short *)malloc(sizeof(unsigned short)
					*(size_t)hp.Ny*(size_t)hp.Nx))==NULL ||
			(dst=(FOM **)malloc(sizeof(FOM *)*(size_t)rows))==NULL ||
			(dst[0]=(FOM *)malloc(sizeof(FOM)
					*(size_t)rows*(size_t)hp.Nx))==NULL)
			Error("memory allocation error for parallel read.");
		    for (yy=1; yy<rows; yy++) dst[yy]=dst[yy-1]+hp.Nx;

		    #pragma omp for schedule(dynamic)
		    for (mm=0; mm<M; mm++) {
			y1a=sel[cs]+z0;
			ReadHiPicBand(&hp,mm,y1a,y1a+rows-1,dst,raw);
			for (kk=cs; kk<ce; kk++) {
			    T2=dst[sel[kk]+z0-y1a]; sg2=SG[kk-cs][mm]+r0;
			    for (rr=0; rr<hp.Nx; rr++) sg2[rr]=(-Log(T2[rr]));
			}

			y1b=sel[cs]+z1;
			ReadHiPicBand(&hp,mm+M,y1b,y1b+rows-1,dst,raw);
			for (kk=cs; kk<ce; kk++) {
			    T2=dst[sel[kk]+z1-y1b]+r1; sg2=SG[kk-cs][mm];
			    for (rr=r2; rr<r3; rr++) sg2[rr]=(-Log(T2[-rr]));
#ifdef	OCT_SBS	/* side by side */
	rr=r4+(r5-r4-1)/2; if ((r5-r4)%2) sg2[rr]=(sg2[rr]-Log(T2[-rr]))/2.0;

	while (++rr<r5) sg2[rr]=(-Log(T2[-rr]));
#else
	for (rr=r4; rr<r5; rr++) sg2[rr]=
#ifdef	OCT_SA	/* simple average */
	(sg2[rr]-Log(T2[-rr]))/2.0;
#else		/* linear mixing */
	((double)(r5-rr)*sg2[rr]+(double)(rr-r4+1)*(-Log(T2[-rr])))/(double)(r5-r4+1);
#endif
#endif
			}
			#pragma omp critical
			fprintf(stderr,"band[%d-%d] read %d / %d\r",
				cs,ce-1,++rdone,M);
		    }
		    free(dst[0]); free(dst); free(raw);
		}
	    }
	    fprintf(stderr,"\n");

	    /* ---- バンド内スライスを再構成 (1 スライス先行のパイプライン) ----
	       BeginCBP(k) を投げた直後に k+1 の準備を行い、EndCBP(k) で結果を
	       受け取ってから BeginCBP(k+1) を投げる。F -> fom のコピーと Store は
	       k+1 の GPU 計算と並行する。バンド境界でパイプラインは空になる。 */
	    cur=0;
	    if (pipe_on) {
		if (PrepareSlice(Pbuf[cur],SG[0],N,M,&kernel_size,&num_threads)) return 5;
		CBP_PIPE_SELECT(Pbuf[cur]); CBP_PIPE_BEGIN(Dr,DO,RA);
	    }

	    for (k=cs; k<ce; k++) {
		z=sel[k];
		if (pipe_on) {
		    if (k+1<ce &&
			PrepareSlice(Pbuf[cur^1],SG[k+1-cs],N,M,&kernel_size,&num_threads))
			return 5;

		    F=CBP_PIPE_END();			/* スライス z の結果 */

		    if (k+1<ce) {
			CBP_PIPE_SELECT(Pbuf[cur^1]); CBP_PIPE_BEGIN(Dr,DO,RA);
		    }
		}
		else {					/* 逐次順(比較用) */
		    if (PrepareSlice(Pbuf[0],SG[k-cs],N,M,&kernel_size,&num_threads)) return 5;
		    CBP_PIPE_SELECT(Pbuf[0]); CBP_PIPE_BEGIN(Dr,DO,RA);
		    F=CBP_PIPE_END();
		}

		if (slice_done++ != 0) TERM_MT(t,Store);

		#pragma omp parallel for
		for (y=0; y<N; y++) {
		    FOM		*ff=fom[y];
		    Float	*sf=F[y];
		    int		xx;

		    for (xx=0; xx<N; xx++) ff[xx]=sf[xx];
		}

		Z=z;
#ifdef WINDOWS
		(void)sprintf(path16,"%s\\rh%05d.tif",argv[7],z);
		(void)sprintf(path8, "%s\\ro%05d.tif",argv[10],z);
#else
		(void)sprintf(path16,"%s/rh%05d.tif",argv[7],z);
		(void)sprintf(path8, "%s/ro%05d.tif",argv[10],z);
#endif
		INIT_MT(t,Store,ac);
		cur^=1;
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

#ifdef USE_GPU
	sort_filter_gpu_release();
#endif
	free(out16); free(out8);
	free(*fom); free(fom); free(**SG); free(*SG); free(SG); TermCBP();

	free(sel);
	free(target);
}
	TermReadHiPic(&hp);

	return 0;
}
