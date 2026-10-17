
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

#ifdef HIS_INPUT
#define PROG_NAME	"his_srec"	/* a.his + conv.bat + output.log via rhis_hp.c */
#define main ofct_srec_main	/* his_main.c: 360-degree offset scan -> here, 180 -> hp_tg_main */
#else
#define PROG_NAME	"ofct_srec"
#endif

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
static char	path[LEN];

/* insert start */
	static double		Dr,DO,RA;
	static int			hpNtM;
	static double		RcZ;	/* axis position (detector pixels) of the slice being stored */
/* insert end */

/* ---- 傾いた回転軸 (2026-10-17) ----
   ofct_srec HiPic L1 C1 L2 C2 Oy rangeList Dr RA0 rec の形では、スライス z の
   軸位置を Rc(z) = C1 + (C2-C1)(z-L1)/(L2-L1) を整数に丸めて使う(貼り合わせは
   画素単位)。幅 N(z) = Nx+|2Rc(z)-Nx| はスライスごとに変わるので、選択スライス
   の最大幅 N を出力サイズとし、狭いスライスのシノグラムは N 幅のバッファの
   中央に置いて左右を行の端の値で埋める(ゼロだと cbp の打ち切り判定 DetectPad
   が両端の列を見て常に「打ち切りなし」とし、幅最大のスライスと補正の有無が
   食い違う。リング除去の列方向フィルタも境界で混ざる)。回転中心は常に画像
   中央で、CBP は N で 1 回だけ初期化する。C1 == C2 なら詰め物がなく従来形式と
   同じ結果になる。ImageDescription の 2 列目にはそのスライスの Rc(z)+0.5 を
   記録する(ofct_rec の dcnt+0.5 と同じ規約。軸が中央より右なら従来の (N+1)/2
   と同値)。

   ---- 小数の軸位置 (2026-10-18) ----
   Rc = ci + f (ci = floor(Rc), 0 <= f < 1)。貼り合わせの枠(幅 N と各領域の境界)は
   整数部 ci で決め、対向像 B だけを 2f 画素右へずらして置く(線形補間、端は
   端値で保持)。像 A は不動なので鏡映の軸は f だけ右へ動き、回転中心を
   (N+1)/2 + f にする。f = 0 では補間の重みが 0 になり従来とビット一致する。
   ofct_rec の f_center = dcnt + 0.5 と同じ規約。 */
/* stitching geometry of ofct_srec for the axis position Rc; 0 = bad offset */
static int	Stitch(int Rc,int Nx,int *n,int *r0,int *r1,int *r2,int *r3,int *r4,int *r5)
{
	int	Ox=2*Rc-Nx;

	if (Ox+Nx<=0 || Nx<=Ox) return 0;
	if (Ox<0) { *n=Nx-Ox; *r0=(-Ox); *r1=Nx-1; *r2=0; *r3=*r4=(-Ox); *r5=Nx; }
	else      { *n=Nx+Ox; *r0=0; *r1=*n-1; *r2=Nx; *r3=*n; *r4=Ox; *r5=Nx; }
	return 1;
}


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

static void Store32TiffFile(char *wname, int wX, int wY, int wBPS, float *data32, char *wdesc)
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
	TIFFSetField(image, TIFFTAG_ARTIST, PROG_NAME);
//	TIFFSetField(image, TIFFTAG_MINSAMPLEVALUE, mmmin );
//	TIFFSetField(image, TIFFTAG_MAXSAMPLEVALUE, mmmax );

		ct_write_raw_strips(image, data32, (uint32_t)wX, (uint32_t)wY, sizeof(float));

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

	static float	*data32=NULL;	/* FOM(double) -> float 変換バッファ(使い回し) */
	static size_t	data32_n=0;
	int Nx, Ny;

	Nx=N;
	Ny=N;

	if (data32_n<(size_t)Nx*(size_t)Ny) {
		free(data32);
		if ((data32 = (float*)malloc(sizeof(float)*(size_t)Nx*(size_t)Ny)) == NULL) {
			printf("cannot allocate memory for output 32bit TIFF image\n");
			exit(1);
		}
		data32_n=(size_t)Nx*(size_t)Ny;
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
	sprintf(comm,"%f\t%f\t%d\t%f\t%lf\t%lf",Dr*10000.0, RcZ, hpNtM, RA, mmmin, mmmax);

//		StoreImageFile_Float(path,N,N,fom,comm);
		Store32TiffFile(path,Nx, Ny, 32, data32, comm);
	    (void)printf("%d\t%s\n",Z,comm);

	free(comm);
	return RETURN_VALUE;
}

static double	Log(double d)
{
	return (d>0.0)?log(d):0.0;
}

/* -log of the opposing view at the (fractional) detector position x,
   linear interpolation, edges held.  w == 0 returns -Log(Tb[i]) exactly,
   so an integer axis position reproduces the old result bit for bit. */
static double	Bval(const FOM *Tb,int Nx,double x)
{
	int	i=(int)floor(x);
	double	w=x-(double)i,a;

	if (i<0) { i=0; w=0.0; } else if (i>=Nx-1) { i=Nx-1; w=0.0; }
	a=-Log(Tb[i]);
	if (w==0.0) return a;
	return (1.0-w)*a+w*(-Log(Tb[i+1]));
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
	int	Oy,L,z0,z1,M,i,Rc0,tilt,L1=0,L2=0,r0,r1,r2,r3,r4,r5;
	double	C1=0.0,C2=0.0,Rcf0;
	char	*a_oy,*a_rl,*a_dr,*a_ra;

    int kernel_size = 5; // Default kernel size
    int num_threads = 40; // Default number of threads

	if (argc!=4 && argc!=8 && argc!=11)
	    Error("usage : ofct_srec HiPic/ Rc Oy {rangeList Dr RA0 rec/}\n"
		  "usage : ofct_srec HiPic/ L1 C1 L2 C2 Oy rangeList Dr RA0 rec/\n"
		  "        (tilted axis: Rc = C1 at layer L1 and C2 at L2, rounded per slice;\n"
		  "         every slice is written at the size of the widest one)"
	    );

	tilt=(argc==11);
	if (tilt) {
	    L1=atoi(argv[2]); C1=atof(argv[3]); L2=atoi(argv[4]); C2=atof(argv[5]);
	    if (L1==L2 && C1!=C2) Error("tilted axis: L1 == L2 but C1 != C2.");
	    Rcf0=C1;
	    a_oy=argv[6]; a_rl=argv[7]; a_dr=argv[8]; a_ra=argv[9];
	}
	else {
	    Rcf0=atof(argv[2]);		/* may be fractional */
	    a_oy=argv[3]; a_rl=argv[4]; a_dr=argv[5]; a_ra=argv[6];
	}
	Rc0=(int)floor(Rcf0);

	InitReadHiPic(argv[1],&hp);

	if (hp.Nt%2) (void)fputs("bad number of views (warning).\n",stderr);

	if (!Stitch(Rc0,hp.Nx,&N,&r0,&r1,&r2,&r3,&r4,&r5) ||
	    (Oy=atoi(a_oy))+hp.Ny<=0 || hp.Ny<=Oy) Error("bad offset.");
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

	RangeList(a_rl,(size_t)L-1,target);

	/* 選択スライスの実 z を sel[] に集約(コンパクト index = SG の行) */
	if ((sel=(int *)malloc(sizeof(int)*L))==NULL)
	    Error("memory allocation error for slice list.");
	nsel=0;
	for (z=0; z<L; z++) if (target[z]) sel[nsel++]=z;

	if (nsel==0) Error("no slice.");

	/* ---- スライスごとの軸位置と貼り合わせ幅(傾き対応) ----
	   最大幅を N(出力サイズ)にし、狭いスライスは pad だけ右へずらして
	   中央に置く。傾きなしでは全スライス同じで pad = 0。 */
	int	*RcS,*NS,*padS,*r0S,*r1S,*r2S,*r3S,*r4S,*r5S,Nmax=0;
	double	*RcF,*fS,*DOS,Rmin=0.0,Rmax=0.0;

	if ((RcS=(int *)malloc(sizeof(int)*9*(size_t)nsel))==NULL ||
	    (RcF=(double *)malloc(sizeof(double)*3*(size_t)nsel))==NULL)
	    Error("memory allocation error for slice geometry.");
	NS=RcS+nsel; padS=NS+nsel; r0S=padS+nsel; r1S=r0S+nsel; r2S=r1S+nsel;
	r3S=r2S+nsel; r4S=r3S+nsel; r5S=r4S+nsel;
	fS=RcF+nsel; DOS=fS+nsel;
	for (i=0; i<nsel; i++) {
	    z=sel[i];
	    RcF[i]=(!tilt)?Rcf0:(L1==L2)?C1:C1+(C2-C1)*(double)(z-L1)/(double)(L2-L1);
	    RcS[i]=(int)floor(RcF[i]);		/* integer part: stitching frame */
	    fS[i]=RcF[i]-(double)RcS[i];	/* fraction: opposing view shifted by 2f, centre by f */
	    if (!Stitch(RcS[i],hp.Nx,NS+i,r0S+i,r1S+i,r2S+i,r3S+i,r4S+i,r5S+i)) {
		(void)fprintf(stderr,"slice %d: axis position %g gives a bad offset.\n",z,RcF[i]);
		Error("bad offset.");
	    }
	    if (NS[i]>Nmax) Nmax=NS[i];
	    if (i==0 || RcF[i]<Rmin) Rmin=RcF[i];
	    if (i==0 || RcF[i]>Rmax) Rmax=RcF[i];
	}
	for (i=0; i<nsel; i++) padS[i]=(Nmax-NS[i])/2;
	N=Nmax;
	if (tilt)
	    (void)fprintf(stderr,"tilted axis: Rc %g .. %g over the selected slices; "
			  "output %d x %d, narrower slices centred and zero-padded\n",Rmin,Rmax,N,N);

	Dr=atof(a_dr)/10000.;
	DO=-(double)(N+1)/2.0;	/* rotation center = -(N+1)/2, matching ofct_rec's r0=-(dcnt+0.5); was -(N-1)/2 (1 px off) */
	for (i=0; i<nsel; i++) DOS[i]=DO-fS[i];	/* fractional part of the axis: centre moves by f */
	RA=atof(a_ra);

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
			    Float	*row=SG[kk-cs][mm];

			    /* narrower slice: centred in the N-wide row; the margins are
			       filled with the row's end values after the second view
			       (below), so that the truncation check of cbp (DetectPad looks
			       at the outermost columns) and the ring filter see the same
			       edges as a slice reconstructed at its own width */
			    T2=dst[sel[kk]+z0-y1a]; sg2=row+padS[kk]+r0S[kk];
			    for (rr=0; rr<hp.Nx; rr++) sg2[rr]=(-Log(T2[rr]));
			}

			y1b=sel[cs]+z1;
			ReadHiPicBand(&hp,mm+M,y1b,y1b+rows-1,dst,raw);
			for (kk=cs; kk<ce; kk++) {
			    int	r1=r1S[kk],r2=r2S[kk],r3=r3S[kk],r4=r4S[kk],r5=r5S[kk];
			    const FOM	*Tb=dst[sel[kk]+z1-y1b];
			    double	f2=2.0*fS[kk];
/* opposing view B at output position rr: detector position r1-rr, shifted by 2f */
#define BV(rr)	Bval(Tb,hp.Nx,(double)(r1-(rr))+f2)

			    sg2=SG[kk-cs][mm]+padS[kk];
			    for (rr=r2; rr<r3; rr++) sg2[rr]=BV(rr);
#ifdef	OCT_SBS	/* side by side */
	rr=r4+(r5-r4-1)/2; if ((r5-r4)%2) sg2[rr]=(sg2[rr]+BV(rr))/2.0;

	while (++rr<r5) sg2[rr]=BV(rr);
#else
	for (rr=r4; rr<r5; rr++) sg2[rr]=
#ifdef	OCT_SA	/* simple average */
	(sg2[rr]+BV(rr))/2.0;
#else		/* linear mixing */
	((double)(r5-rr)*sg2[rr]+(double)(rr-r4+1)*BV(rr))/(double)(r5-r4+1);
#endif
#endif
#undef BV
			    {	/* margins: hold the end values of the stitched row */
				Float	*row=SG[kk-cs][mm],vl=row[padS[kk]],vr=row[padS[kk]+NS[kk]-1];
				for (rr=0; rr<padS[kk]; rr++) row[rr]=vl;
				for (rr=padS[kk]+NS[kk]; rr<N; rr++) row[rr]=vr;
			    }
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
		CBP_PIPE_SELECT(Pbuf[cur]); CBP_PIPE_BEGIN(Dr,DOS[cs],RA);
	    }

	    for (k=cs; k<ce; k++) {
		z=sel[k];
		if (pipe_on) {
		    if (k+1<ce &&
			PrepareSlice(Pbuf[cur^1],SG[k+1-cs],N,M,&kernel_size,&num_threads))
			return 5;

		    F=CBP_PIPE_END();			/* スライス z の結果 */

		    if (k+1<ce) {
			CBP_PIPE_SELECT(Pbuf[cur^1]); CBP_PIPE_BEGIN(Dr,DOS[k+1],RA);
		    }
		}
		else {					/* 逐次順(比較用) */
		    if (PrepareSlice(Pbuf[0],SG[k-cs],N,M,&kernel_size,&num_threads)) return 5;
		    CBP_PIPE_SELECT(Pbuf[0]); CBP_PIPE_BEGIN(Dr,DOS[k],RA);
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

		RcZ=RcF[k]+0.5;		/* after the join: Store reads it. Rc+0.5 = the centre
					   convention of ofct_rec (dcnt+0.5) and of the old
					   ofct_srec ((N+1)/2) for an axis right of the middle */
#ifdef WINDOWS
		(void)sprintf(path,"%s\\rec%05d.tif",argv[argc-1],Z=z);
#else
		(void)sprintf(path,"%s/rec%05d.tif",argv[argc-1],Z=z);
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
#ifdef HIS_INPUT
	{ extern const char *HisScanNote(void); fprintf(f,"%s",HisScanNote()); }
#endif
	fprintf(f,"\n");
	fclose(f);


	TERM_MT(t,Store);

#ifdef USE_GPU
	sort_filter_gpu_release();
#endif
	free(*fom); free(fom); free(**SG); free(*SG); free(SG); TermCBP();

	free(RcF); free(RcS);
	free(sel);
	free(target);
}
	TermReadHiPic(&hp);

	return 0;
}
