/* ============================================================================
 * p_rec.c  -  32bit float 投影像 (p?????.tif) からの連続再構成
 * ----------------------------------------------------------------------------
 *   p_rec p/ rec/ Dr RC RA0
 *   p_rec p/ rec/ Dr L1 C1 L2 C2 RA0
 *
 * 入力: p/ の p%05d.tif (ct_prj_f などが出力する -log 透過率、32bit float、
 *       ストリップ形式)。番号は 0 または 1 始まりのどちらでもよく、見つかった
 *       最小番号から最大番号までを投影 0..Nt-1 として使う。
 * 出力: rec/rec%05d.tif (32bit float、単位変換 ×10000/Dr 済み)。
 *
 * メモリ上限チェックとバンド(複数パス)実行は hp_tg と同じ
 * (HPTG_MEM_FRACTION / HPTG_MEM_LIMIT_MB / HPTG_CHUNK_ROWS)。
 *
 * 2026-10 の見直し:
 *   - 旧版は読み込みループの番号付けにより最初のファイルを読まず、投影 0 が
 *     常にゼロ(空)になっていた。全ファイルを使うよう修正(結果が変わる唯一の点)。
 *   - 旧版の異常値除外は int 版 abs() を float に適用していたため機能して
 *     いなかった(NaN/Inf は通過、|v|>=101 は前スライスの値が残る)。
 *     fabsf で判定し、除外画素は 0 にする。通常データでは結果不変。
 *   - 読み込みを投影ファイル単位で並列化(HPTG_READ_THREADS、既定 16)し、
 *     バンドに必要な行のストリップだけを読む(multi-pass でも全行を読まない)。
 *   - リング除去は CBP の投影バッファ上でインプレース実行、出力バッファは
 *     1 回だけ確保、単位変換と min/max は OpenMP 並列、TIFF 書き出しは
 *     別スレッド、1 スライス先行パイプライン(cbp_pipe.h)。hp_tg と同じ構造。
 *     これらは各要素の演算順序を変えないので結果はビット単位で不変。
 *   - HPTG_PIPELINE=0 で先行パイプラインを止めて逐次順で実行(比較用)。
 * ==========================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "tiffio.h"
#include "tifwrite.h"
#include "cbp.h"
#include "cbp_pipe.h"
#ifdef USE_GPU
  #include "sort_filter_g.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_gpu
#else
  #include "sort_filter_omp.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_omp
#endif

#ifdef _WIN32
#include <windows.h>
#include <process.h>

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

/* 異常値除外のしきい値: |v| がこれ以上の画素は 0 として扱う(旧版の意図を踏襲) */
#define P_REJECT	100.0f

static int	Nx, Ny, Nt;

static void Error(char *msg)
{
	fputs(msg, stderr);
	fputc('\n', stderr);
	exit(1);
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

/*----------------------------------------------------------------------*/

static int existFile(const char* path)
{
	FILE* fp = fopen(path, "r");
	if (fp == NULL) return 0;
	fclose(fp);
	return 1;
}

static void Store32TiffFile(char *wname, int wX, int wY, float *data32, char *wdesc)
{
	TIFF *image;

	if ((image = TIFFOpen(wname, "w"))==NULL) Error("cannot open output TIFF.");

	TIFFSetField(image, TIFFTAG_IMAGEWIDTH, wX);
	TIFFSetField(image, TIFFTAG_IMAGELENGTH, wY);
	TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(image, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP );
	TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, wdesc);
	TIFFSetField(image, TIFFTAG_ARTIST, "p_rec");

	ct_write_raw_strips(image, data32, (uint32_t)wX, (uint32_t)wY, sizeof(float));

	TIFFClose(image);
}

/* ---- 先頭ファイルから Nx, Ny を得る ---- */
static void ReadTifHeader(const char *rname)
{
	TIFF		*image;
	uint32_t	w=0,h=0;
	uint16_t	bps=0,sf=SAMPLEFORMAT_UINT;

	if ((image = TIFFOpen(rname, "r"))==NULL) Error("cannot open the first p-file.");
	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &w);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, &h);
	TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &bps);
	TIFFGetFieldDefaulted(image, TIFFTAG_SAMPLEFORMAT, &sf);
	TIFFClose(image);
	if (bps!=32 || sf!=SAMPLEFORMAT_IEEEFP) Error("p-files must be 32bit float TIFF.");
	Nx=(int)w; Ny=(int)h;
}

/* ---- p ファイルの行 [y0,y1] だけを読む ----
   dst + (y-y0)*stride に行 y を置く。必要なストリップだけを
   TIFFReadEncodedStrip で読むので、バンド分割時も全行は読まない。
   スレッド安全(ファイル毎に独立なハンドルとバッファ)。 */
static void ReadTifRows(const char *rname, int y0, int y1, Float *dst, size_t stride)
{
	TIFF		*image;
	uint32_t	w=0,h=0,rps=0;
	tmsize_t	ssize;
	float		*buf;
	int		y;

	if ((image = TIFFOpen(rname, "r"))==NULL) {
	    fprintf(stderr, "cannot open %s\n", rname); exit(1);
	}
	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &w);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, &h);
	if ((int)w!=Nx || (int)h!=Ny) {
	    fprintf(stderr, "%s: image size mismatch (%ux%u, expected %dx%d)\n",
		    rname, w, h, Nx, Ny); exit(1);
	}
	TIFFGetFieldDefaulted(image, TIFFTAG_ROWSPERSTRIP, &rps);
	if (rps==0 || rps>h) rps=h;
	ssize=TIFFStripSize(image);
	if ((buf=(float *)_TIFFmalloc(ssize))==NULL) Error("no memory for a TIFF strip.");

	for (y=y0; y<=y1; ) {
	    tstrip_t	strip=(tstrip_t)(y/(int)rps);
	    int		first=(int)(strip*rps),
			last=first+(int)rps-1;
	    tmsize_t	n;

	    if (last>Ny-1) last=Ny-1;
	    if ((n=TIFFReadEncodedStrip(image, strip, buf, ssize))<0) {
		fprintf(stderr, "%s: cannot read strip %u\n", rname, (unsigned)strip); exit(1);
	    }
	    if (n<(tmsize_t)(last-first+1)*(tmsize_t)Nx*(tmsize_t)sizeof(float)) {
		fprintf(stderr, "%s: short strip %u\n", rname, (unsigned)strip); exit(1);
	    }
	    for (; y<=last && y<=y1; y++)
		memcpy(dst+(size_t)(y-y0)*stride, buf+(size_t)(y-first)*(size_t)Nx,
		       sizeof(float)*(size_t)Nx);
	}
	_TIFFfree(buf);
	TIFFClose(image);
}

/* ---- 1 スライス分のホスト側準備: po_band のスライス -> Pb、異常値除外、
   リング除去(インプレース)。GPU 版では直前スライスの BeginCBP の後に呼ばれ、
   GPU 計算と並行して走る。 */
static int	PrepareSlice(Float **Pb, const Float *slice, int *kernel_size, int *num_threads)
{
	int	l;

	#pragma omp parallel for
	for (l=0; l<Nt; l++) {
	    const Float	*src=slice+(size_t)l*(size_t)Nx;
	    Float	*pf=Pb[l];
	    int		n;

	    for (n=0; n<Nx; n++) {
		Float v=src[n];
		pf[n]=(fabsf(v)<P_REJECT)?v:0.0f;	/* NaN/Inf/外れ値は 0 */
	    }
	}

	*kernel_size = get_kernel_size_from_env();
	*num_threads = get_num_threads_from_env();
	/* Pb 上でインプレース実行(入力=出力)。GPU 版は H2D 後に D2H、
	   CPU 版は列ごとに全行を読んでから書くので、どちらも安全。 */
	if (SORT_FILTER_RESTORE(Pb[0], Pb[0], Nx, Nt, *kernel_size, *num_threads) != 0) {
	    fprintf(stderr, "ring removal failed\n");
	    return 5;
	}
	return 0;
}

/* ---- 書き出しスレッド ---- */
typedef struct {
	char	path[LEN];
	float	*out;
	int	N;
	char	comm[150];
} StoreArg;

static FUNCTION_T	Store(void *a)
{
	StoreArg	*s=(StoreArg *)a;

	Store32TiffFile(s->path, s->N, s->N, s->out, s->comm);
	return RETURN_VALUE;
}

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC	1000000
extern long clock();
#endif

#define CLOCK()		((double)clock()/(double)CLOCKS_PER_SEC)

int	main(int argc, char *argv[])
{
	long		i, p_sta, p_dst;
	Float		**P, **F, **Pbuf[2], *po_band, *out32;
	char		fh[LEN];
	int		z1, z2, cur, pipe_on=1, slice_done=0;
	double		Dr, RC, RA0, Ct, RCcur;
	double		t1;
	THREAD_T	T;
	StoreArg	S;

	/* ---- chunk 関連 ---- */
	int		rows_per_chunk, maxrows, nbands;
	long long	total_rows;

	int kernel_size = 5;	/* 既定(環境変数 KERNEL_SIZE で上書き) */
	int num_threads = 40;

	if (argc != 6 && argc != 9) {
		fprintf(stderr, "usage : p_rec p/ rec/ Dr RC RA0 \nusage : p_rec p/ rec/ Dr L1 C1 L2 C2 RA0\n");
		return 1;
	}

	/* ---- p ファイルの番号範囲 (0 始まりも 1 始まりも可) ---- */
	p_sta = -1; p_dst = -1;
	for (i = 0; i<100000; i++) {
		sprintf(fh, "%s/p%05ld.tif", argv[1], i);
		if (existFile(fh)) {
			if (p_sta == -1) p_sta = i;
			p_dst = i;
		}
	}
	if (p_sta<0) Error("no p?????.tif in the input directory.");
	sprintf(fh, "%s/p%05ld.tif", argv[1], p_sta);
	ReadTifHeader(fh);
	Nt=(int)(p_dst-p_sta+1);

	z1=0; z2=Ny; Ct=0.0;
	Dr = atof(argv[3]);
	if (argc == 6) {
		RC = atof(argv[4]);
		RA0 = atof(argv[5]);
	}
	else {
		z1 = atoi(argv[4]);
		RC = atof(argv[5]);
		z2 = atoi(argv[6]);
		Ct = (atof(argv[7])-atof(argv[5]))/(double)(z2-z1);
		RA0 = atof(argv[8]);
	}
	if (z1<0) z1=0;
	if (z2>Ny) z2=Ny;
	printf("%d\t%d\t%d\t%d\t%d\n", Nx, Ny, Nt, z1, z2-1);

	/* ============================================================
	 *  メモリ上限チェック -> 1 バンドあたりの行(スライス)数を決定
	 *  巨大配列 po_band = (行数) x Nt x Nx x sizeof(Float)
	 *  -> 行方向に分割し、各バンドで必要な行だけを p ファイル群から読む。
	 * ============================================================ */
	unsigned long long mem_avail = get_available_memory_bytes();
	unsigned long long mem_total = get_total_memory_bytes();
	unsigned long long mem_ref   = mem_avail ? mem_avail : mem_total;

	double frac=0.9;
	{ const char *e=getenv("HPTG_MEM_FRACTION"); if (e){ double v=atof(e); if (v>0.05 && v<=0.95) frac=v; } }

	unsigned long long bytes_per_row = (unsigned long long)Nt*(unsigned long long)Nx*sizeof(Float);

	/* po_band 以外の固定オーバヘッド見積り */
	unsigned long long overhead =
	      2ull*(unsigned long long)Nt*Nx*sizeof(Float)       /* P x 2 面     */
	    + (unsigned long long)Nx*Nx*sizeof(float)            /* 出力バッファ */
	    + 16ull*(unsigned long long)Nx*sizeof(float)*64;     /* 読み込みストリップ(概算) */

	unsigned long long budget = (unsigned long long)((double)mem_ref*frac);
	if (budget>overhead) budget -= overhead;
	else                 budget = bytes_per_row;	/* 最低 1 行 */

	{ const char *e=getenv("HPTG_MEM_LIMIT_MB");
	  if (e){ double mb=atof(e); if (mb>0) budget=(unsigned long long)(mb*1024.0*1024.0); } }

	total_rows = (long long)z2 - (long long)z1;
	if (total_rows < 1) total_rows = 1;

	{ const char *e=getenv("HPTG_CHUNK_ROWS");
	  if (e && atoi(e)>0) rows_per_chunk=atoi(e);
	  else {
	      unsigned long long rr = bytes_per_row ? (budget/bytes_per_row) : (unsigned long long)total_rows;
	      if (rr<1) rr=1;
	      rows_per_chunk = (rr>(unsigned long long)total_rows)? (int)total_rows : (int)rr;
	  }
	}
	if (rows_per_chunk<1) rows_per_chunk=1;
	if ((long long)rows_per_chunk>total_rows) rows_per_chunk=(int)total_rows;

	maxrows = rows_per_chunk;
	nbands  = (int)((total_rows + rows_per_chunk - 1) / rows_per_chunk);

	fprintf(stderr,
	    "memory: avail=%.2f GB  total=%.2f GB  po/row=%.2f MB  rows/band=%d  bands=%d  (%s)\n",
	    mem_avail/1073741824.0, mem_total/1073741824.0,
	    bytes_per_row/1048576.0, rows_per_chunk, nbands,
	    (nbands>1)?"multi-pass":"single-pass");

	/* po_band(バンド分)= maxrows x Nt x Nx、スライス主(スライス m の投影が連続) */
	if ((po_band = (Float *)malloc(sizeof(Float)*(size_t)maxrows*(size_t)Nt*(size_t)Nx)) == NULL)
		Error("cannot allocate memory for the projection band.");

	/* 出力バッファ(全スライスで使い回し。Store の完了を待ってから上書きする) */
	if ((out32 = (float *)malloc(sizeof(float)*(size_t)Nx*(size_t)Nx)) == NULL)
		Error("cannot allocate memory for the output image.");
	S.out=out32; S.N=Nx;

	/* ---- CBP と投影バッファ 2 面(1 スライス先行パイプライン) ---- */
	if ((P=InitCBP(Nx,Nt))==NULL) Error("memory allocation error.");
	CBP_PIPE_INIT(P,Nx,Nt);
	Pbuf[0]=P;
	if ((Pbuf[1]=CBP_PIPE_ALLOC())==NULL)
	    Error("memory allocation error for the second projection buffer.");
	if (sizeof(Float)!=sizeof(float))
	    Error("in-place ring removal requires Float==float.");
	for (i=0; i<2; i++) {
	    int l;
	    for (l=1; l<Nt; l++)
		if (Pbuf[i][l]!=Pbuf[i][l-1]+Nx)
		    Error("projection buffer is not contiguous.");
	}
	{ char *e=getenv("HPTG_PIPELINE"); if (e && atoi(e)==0) pipe_on=0; }
	if (!pipe_on) fprintf(stderr,"HPTG_PIPELINE=0: sequential order\n");

	/* ---- 読み込みスレッド数 (hp_tg と同じ環境変数) ---- */
	int rthreads=16;
	{ char *e=getenv("HPTG_READ_THREADS"); if (e && atoi(e)>0) rthreads=atoi(e); }
	if (rthreads>Nt) rthreads=Nt;

	/* ============================================================
	 *  バンドループ(必要に応じて複数パス)
	 * ============================================================ */
	for (int cs=z1; cs<z2; cs+=rows_per_chunk) {
		int	ce = cs+rows_per_chunk; if (ce>z2) ce=z2;
		int	rdone=0, l, m;

		/* ---- このバンドの行 [cs,ce) を全 p ファイルから並列に読む ---- */
		t1=CLOCK();
		#pragma omp parallel for schedule(dynamic) num_threads(rthreads)
		for (l=0; l<Nt; l++) {
		    char	path[LEN];

		    sprintf(path, "%s/p%05ld.tif", argv[1], p_sta+(long)l);
		    /* 行 m -> po_band[((m-cs)*Nt + l)*Nx] */
		    ReadTifRows(path, cs, ce-1,
				po_band+(size_t)l*(size_t)Nx, (size_t)Nt*(size_t)Nx);
		    #pragma omp critical
		    fprintf(stderr, "\rband[%d-%d] read %d / %d", cs, ce-1, ++rdone, Nt);
		}
		fprintf(stderr, "\t%.1f s\n", CLOCK()-t1);

		/* ---- スライス cs..ce-1 を再構成 (1 スライス先行のパイプライン) ----
		   BeginCBP(m) を投げた直後に m+1 の準備を行い、EndCBP(m) で結果を
		   受け取ってから BeginCBP(m+1) を投げる。単位変換と Store は m+1 の
		   GPU 計算と並行する。バンド境界でパイプラインは空になる。 */
#define SLICE(mm)	(po_band+(size_t)((mm)-cs)*(size_t)Nt*(size_t)Nx)
		cur=0;
		if (pipe_on) {
		    if (PrepareSlice(Pbuf[cur],SLICE(cs),&kernel_size,&num_threads)) return 5;
		    CBP_PIPE_SELECT(Pbuf[cur]); CBP_PIPE_BEGIN(1.0,-RC,RA0);
		}
		for (m=cs; m<ce; m++) {
		    double	Clock=CLOCK(), data_max, data_min;
		    int		vv;

		    RCcur=RC;			/* このスライスの回転中心 (= C1+(m-z1)*Ct) */
		    if (pipe_on) {
			if (m+1<ce &&
			    PrepareSlice(Pbuf[cur^1],SLICE(m+1),&kernel_size,&num_threads)) return 5;

			F=CBP_PIPE_END();	/* スライス m の結果 */

			RC=RC+Ct;		/* 次スライス(次バンド先頭を含む)の回転中心 */
			if (m+1<ce) {
			    CBP_PIPE_SELECT(Pbuf[cur^1]); CBP_PIPE_BEGIN(1.0,-RC,RA0);
			}
		    }
		    else {			/* 逐次順(比較用) */
			if (PrepareSlice(Pbuf[0],SLICE(m),&kernel_size,&num_threads)) return 5;
			CBP_PIPE_SELECT(Pbuf[0]); CBP_PIPE_BEGIN(1.0,-RC,RA0);
			F=CBP_PIPE_END();
			RC=RC+Ct;
		    }

		    /* 前スライスの書き出しが out32 を読み終わるのを待つ */
		    if (slice_done++ != 0) TERM_MT(T);

		    /* 単位変換 (um -> cm) と min/max。行ごとに並列、最後に逐次で畳む
		       (min/max は順序に依らないので結果は不変) */
		    {
			double	*rmin=(double *)malloc(sizeof(double)*(size_t)Nx),
				*rmax=(double *)malloc(sizeof(double)*(size_t)Nx);

			if (rmin==NULL || rmax==NULL) Error("no memory for min/max.");
			#pragma omp parallel for
			for (vv=0; vv<Nx; vv++) {
			    float	*o=out32+(size_t)Nx*(size_t)vv;
			    Float	*f=F[vv];
			    double	mn=32000., mx=-32000.;
			    int		hh;

			    for (hh=0; hh<Nx; hh++) {
				o[hh] = f[hh]*10000./Dr;
				if (mx<o[hh]) mx=o[hh];
				if (mn>o[hh]) mn=o[hh];
			    }
			    rmin[vv]=mn; rmax[vv]=mx;
			}
			data_max=-32000.; data_min=32000.;
			for (vv=0; vv<Nx; vv++) {
			    if (data_max<rmax[vv]) data_max=rmax[vv];
			    if (data_min>rmin[vv]) data_min=rmin[vv];
			}
			free(rmin); free(rmax);
		    }

		    sprintf(S.comm,"%f\t%f\t%d\t%f\t%f\t%f",Dr, RCcur, Nt, RA0, (float)data_min, (float)data_max);
#ifdef WINDOWS
		    sprintf(S.path, "%s\\rec%05d.tif", argv[2], m);
#else
		    sprintf(S.path, "%s/rec%05d.tif", argv[2], m);
#endif
		    INIT_MT(T,Store,&S);
		    fprintf(stderr, "\rstore:\t%s\t%.3f s", S.path, CLOCK()-Clock);
		    cur^=1;
		}
#undef SLICE
	}
	if (slice_done) TERM_MT(T);

	printf("\nfinish.\n");
	free(po_band);
	free(out32);
#ifdef USE_GPU
	sort_filter_gpu_release();
#endif
	TermCBP();

	// append to log file
	FILE		*ff;
	if ((ff = fopen("cmd-hst.log", "a")) == NULL) {
		return(-1);
	}
	for (i = 0; i<argc; ++i) fprintf(ff, "%s ", argv[i]);
	fprintf(ff, "   %% kernel_size %d", kernel_size);
	fprintf(ff, "   %% bands %d", nbands);
	fprintf(ff, "\n");
	fclose(ff);

	return 0;
}
