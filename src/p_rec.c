#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "tiffio.h"
#include "cbp.h"
#ifdef USE_GPU
  #include "sort_filter_g.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_gpu
#else
  #include "sort_filter_omp.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_omp
#endif

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

#define MA(cnt,ptr)	malloc((cnt)*sizeof(*(ptr)))

static long long	Nx, Ny, Nt, M;
static int			BPS;
static char			*desc;
static Float		*data32;

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

int existFile(const char* path)
{
	FILE* fp = fopen(path, "r");
	if (fp == NULL) {
		return 0;
	}

	fclose(fp);
	return 1;
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
	TIFFSetField(image, TIFFTAG_ARTIST, "p_rec");
//	TIFFSetField(image, TIFFTAG_MINSAMPLEVALUE, mmmin );
//	TIFFSetField(image, TIFFTAG_MAXSAMPLEVALUE, mmmax );

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data32 + i*wX, wX * sizeof(float));
	}

	TIFFClose(image);
}

void Read32TiffFile(char* rname, int iHead)
{
	TIFF* image;
	long i, j;
	unsigned short spp, com, pm, pc, rps;
	float* rline;

	image = TIFFOpen(rname, "r");

	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &Nx);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, &Ny);
	TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &BPS);
//	TIFFGetField(image, TIFFTAG_COMPRESSION, &com);
//	TIFFGetField(image, TIFFTAG_PHOTOMETRIC, &pm);
//	TIFFGetField(image, TIFFTAG_SAMPLESPERPIXEL, &spp);
//	TIFFGetField(image, TIFFTAG_ROWSPERSTRIP, &rps);
//	TIFFGetField(image, TIFFTAG_PLANARCONFIG, &pc);
	TIFFGetField(image, TIFFTAG_IMAGEDESCRIPTION, &desc);

	if(iHead==1){
		TIFFClose(image);
		return;
	}
	if ((rline = (float*)_TIFFmalloc(TIFFScanlineSize(image))) == NULL) {
		printf("cannot allocate memory for line scan\n");
		exit(1);
	}
//		fprintf(stderr, "1\r");

	for (i = 0; i < Ny; i++) {
		if (TIFFReadScanline(image, rline, i, 0) < 0) {
			printf("cannot get tif line -> %d\n", i);
			exit(1);
		}
		for (j = 0; j < Nx; j++) {
			*(data32 + i * Nx + j) = *(rline + j);
		}
	}

	_TIFFfree(rline);
	TIFFClose(image);
	return;
}
#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC	1000000
extern long clock();
#endif

#define CLOCK()		((double)clock()/(double)CLOCKS_PER_SEC)

int	main(int argc, char *argv[])
{
	long long	l, m, n;
	long		i, j, p_sta, p_dst;
	Float		**P, **F;			// full size of reconstructed image
	char		fh[256], fo[256];
	int			z1, z2;
	Float		*po_band, *out32;
	int			vv, hh;
	char		*comm = NULL;
	double		data_max, data_min;
	double		Dr,RC,RA0,Ct;
	double		Clock, t1,t2,t3;					// timer setting

	/* ---- chunk 関連 ---- */
	int			rows_per_chunk, maxrows, nbands;
	long long	total_rows;

    int kernel_size = 5; // Default kernel size
    int num_threads = 40; // Default number of threads

//	printf("%d\n", argc);
	Clock=CLOCK();
	if (argc == 6 || argc == 9 ) {
		p_sta = -1;
		p_dst = -1;
		for (i = 1; i<100000; i++) {
			sprintf(fh, "%s/p%05ld.tif", argv[1], i);
			if (existFile(fh)) {
				if (p_sta == -1) {
					p_sta = i;
					Read32TiffFile(fh,1);
					fprintf(stderr, "%s\r", fh);
				}
				else {
					p_dst = i;
				}
			}
		}
		Nt=p_dst-p_sta+1;
		z1=0;
		z2=Ny;
		Ct=0.0;
		if (argc == 6) {
			Dr = atof(argv[3]);
			RC = atof(argv[4]);
			RA0 = atof(argv[5]);
		}
		if (argc == 9) {
			Dr = atof(argv[3]);
			z1 = atoi(argv[4]);
			RC = atof(argv[5]);
			z2 = atoi(argv[6]);
			Ct = (atof(argv[7])-atof(argv[5]))/(double)(z2-z1);
			RA0 = atof(argv[8]);
		}
	}
	else {
		fprintf(stderr, "usage : p_rec p/ rec/ Dr RC RA0 \nusage : p_rec p/ rec/ Dr L1 C1 L2 C2 RA0\n");
		return 1;
	}
	printf("%lld\t%lld\t%lld\t%d\t%d\n", Nx, Ny, Nt,z1,z2-1);

	/* ============================================================
	 *  メモリ上限チェック -> 1 バンドあたりの行(スライス)数を決定
	 *  巨大配列 po(全投影ボリューム) = Nt x (行数) x Nx x sizeof(float)
	 *  -> 行方向に分割し、各バンドで p ファイル群を読み直す。
	 * ============================================================ */
	unsigned long long mem_avail = get_available_memory_bytes();
	unsigned long long mem_total = get_total_memory_bytes();
	unsigned long long mem_ref   = mem_avail ? mem_avail : mem_total;

	double frac=0.8;
	{ const char *e=getenv("HPTG_MEM_FRACTION"); if (e){ double v=atof(e); if (v>0.05 && v<=0.95) frac=v; } }

	unsigned long long bytes_per_row = (unsigned long long)Nt*(unsigned long long)Nx*sizeof(float);

	/* po 以外の固定オーバヘッド見積り */
	unsigned long long overhead =
	      (unsigned long long)Nt*Nx*sizeof(Float)            /* P            */
	    + (unsigned long long)Nx*Ny*sizeof(float)            /* 読み込みバッファ data32 */
	    + (unsigned long long)Nx*Nx*sizeof(float)            /* 出力バッファ out32 */
	    + 2ull*(unsigned long long)Nx*Nt*sizeof(float);      /* ring 作業    */

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

	/* po(バンド分)= Nt x maxrows x Nx */
	if ((po_band = (Float *)malloc(sizeof(Float)*(size_t)Nt*maxrows*Nx)) == NULL) {
		printf("cannot allocate memory for projection band.\n");
		return 1;
	}

	/* 読み込みバッファ(1 投影 = Nx*Ny)。Read32TiffFile が使うグローバル data32 */
	if ((data32 = (float *)malloc(sizeof(float) * Nx * Ny)) == NULL) {
		printf("cannot allocate memory for input 32bit TIFF image\n");
		return 1;
	}

	// initilaize for CBP
	if ((P=InitCBP(Nx,Nt))==NULL){
		Error("memory allocation error.");
	}

	/* ============================================================
	 *  バンドループ(必要に応じて複数パス)
	 * ============================================================ */
	for (long long cs=z1; cs<z2; cs+=rows_per_chunk) {
		long long ce = cs+rows_per_chunk; if (ce>z2) ce=z2;

		/* バンドをゼロ初期化(未読スロット 0..p_sta-1 を 0 に保つ:元コードと同じ) */
		memset(po_band, 0, sizeof(Float)*(size_t)Nt*maxrows*Nx);

		// store p-data (このバンドの行 [cs,ce) のみ) from float tiff files
		t1=CLOCK();
		for(l=p_sta;l<p_dst;++l) {
			sprintf(fh, "%s/p%05lld.tif", argv[1], (long long)l+1);
			fprintf(stderr, "\rband[%lld-%lld] read:\t%s\t", cs, ce-1, fh);
			(void)Read32TiffFile(fh,0);

			for(m=cs;m<ce;m++) {
				for(n=0;n<Nx;n++) {
					*(po_band + ((size_t)l*maxrows + (m-cs))*Nx + n) = *(data32 + m*Nx + n);
				}
			}
		}
		fprintf(stderr, "\t%lf\n",CLOCK()-t1);

		// loop cs to ce
		for(m=cs;m<ce;++m){
			for (l=0;l<Nt;l++){
				for (n=0; n<Nx; n++){
					if(abs(*(po_band + ((size_t)l*maxrows + (m-cs))*Nx + n))<100.){
						P[l][n]=*(po_band + ((size_t)l*maxrows + (m-cs))*Nx + n);
					}
				}
			}

/* ----------------  ring removal start ---------------- */
/*                                                       */
    float		*image_data = NULL, *result_data = NULL;
    // Get kernel size from environment variable
    kernel_size = get_kernel_size_from_env();
	// Get number of threads from environment variable
    num_threads = get_num_threads_from_env();
    // Allocate memory
	image_data = (float *)malloc(Nx * Nt * sizeof(float));
	result_data = (float *)malloc(Nx * Nt * sizeof(float));

	for (j=0; j<Nt; j++){
		for (i=0; i<Nx; i++){
			*(image_data+Nx*j+i)=P[j][i];
		}
	}
	// Execute OpenMP image processing
	if (SORT_FILTER_RESTORE(image_data, result_data, Nx, Nt, kernel_size, num_threads) != 0) {
		fprintf(stderr, "OpenMP image processing failed\n");
		return 5;
	}
	for (j=0; j<Nt; j++){
		for (i=0; i<Nx; i++){
				P[j][i]=*(result_data+Nx*j+i);
		}
	}
    if (image_data) free(image_data);
    if (result_data) free(result_data);
/* ----------------  ring removal finish --------------- */
/*                                                       */

// CBP
			Clock=CLOCK();
			F=CBP(1.0,-RC,RA0);
			t2=CLOCK()-Clock;

// Store CT images
			out32 = (float *)malloc(Nx*Nx*sizeof(float));
			data_max =-32000.;
			data_min = 32000.;

			Clock=CLOCK();
			for(vv=0; vv<Nx; vv++){
				for (hh=0; hh<Nx; hh++){
					*(out32+Nx*vv+hh) = F[vv][hh]*10000./Dr;	/* unit change  um -> cm */;
					if(data_max<*(out32+Nx*vv+hh)) data_max=*(out32+Nx*vv+hh);
					if(data_min>*(out32+Nx*vv+hh)) data_min=*(out32+Nx*vv+hh);
				}
			}
			if ((comm=(char *)malloc(150))==NULL){
				Error("comment memory allocation error.");
			}

			sprintf(comm,"%f\t%f\t%lld\t%f\t%f\t%f",Dr, RC, Nt, RA0, (float)data_min, (float)data_max);
#ifdef WINDOWS
			sprintf(fo, "%s\\rec%05lld.tif", argv[2], m);
#else
			sprintf(fo, "%s/rec%05lld.tif", argv[2], m);
#endif
			(void)Store32TiffFile(fo, Nx, Nx, 32, out32, comm);
			free(out32); free(comm);
			t3=CLOCK()-Clock;
			fprintf(stderr, "\rstore:\t%s/rec%05lld.tif\t%lf\t%lf", argv[2], m, t2, t3);

			RC=RC+Ct;
		}
	}

	printf("\nfinish.\n");
	free(po_band);
	free(data32);
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
