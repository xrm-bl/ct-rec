/*
 * tif2his.c - a%06d.tif の連番を a.his (HiPic, type 2) にパックする。
 *             パック時に median / gaussian フィルタをかけられる (tif_mgf 相当)。
 *
 *   tif2his [median_kernel_size] [gaussian_sigma]
 *     median_kernel_size : 0-25      (省略時 0 = かけない; 0,1 はかけない)
 *     gaussian_sigma     : 0.0-100.0 (省略時 0.0 = かけない)
 *
 *   例) tif2his          … フィルタなしでパック
 *       tif2his 3 0.5    … tif_mgf in out 3 0.5 を全フレームにかけてからパック
 *
 * フィルタ部 (median → gaussian、ミラー境界、double 演算、+0.5 で丸め) は
 * ct-rec の tif_mgf.c のコードをそのまま使っているので、各フレームの結果は
 * tif_mgf で 1 枚ずつ処理した場合とビット単位で一致する。
 * 引数の意味と範囲も tif_mgf と同じ。ただし省略時はフィルタなし
 * (tif_mgf の既定値 median 1 / sigma 0.5 とは異なる)。
 *
 * 入力: カレントディレクトリの a000000.tif または a000001.tif から始まる
 *       連番 (番号が途切れたところで終了)。非圧縮・ストリップ形式・
 *       モノクロ 8/16bit 符号なし整数。II / MM どちらも可。libtiff 不要。
 * 出力: a.his。各フレーム = 64 バイトヘッダ ("IM", コメント長 0, 幅, 高さ,
 *       オフセット 0, type) + 画素 (リトルエンディアン)。先頭フレームの
 *       offset 14 に総フレーム数 (32bit)。rhis.c と同じ解釈。
 *       a.his が既にあれば上書きしない。実行内容は cmd-hst.log に追記。
 *       ※ rhis (his_rec / his_srec 等) が読めるのは 16bit (type 2) のみ。
 *
 * Build (Linux):   gcc -O3 -march=native -fno-math-errno -fopenmp -o tif2his tif2his.c -lm
 * Build (Windows): cl /O2 /openmp tif2his.c
 *   tif_mgf と同様、-ffast-math / -Ofast は付けないこと (gaussian の加算順が
 *   変わり、tif_mgf と結果が一致しなくなる)。
 */
#define _CRT_SECURE_NO_WARNINGS
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

/*======================================================================*/
/* ここから tif_mgf.c (ct-rec) のフィルタ部をそのまま流用                 */
/*======================================================================*/

#ifdef _OPENMP
#include <omp.h>
#endif

#define MEDIAN_KMAX	25		/* CLI limit for the median kernel size */
#define GAUSS_KMAX	51		/* calculate_kernel_size() caps here */

/*----------------------------------------------------------------------*/
/* Compiler helpers.  None of these change any arithmetic; they only tell
   the compiler what it may assume, so that the straight-line min/max code
   of the median networks can be vectorised across neighbouring pixels. */

/* no aliasing between the source and destination images */
#if	defined(__GNUC__) || defined(_MSC_VER)
#define RESTRICT	__restrict
#elif	defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#define RESTRICT	restrict
#else
#define RESTRICT
#endif

/* the window networks must be inlined into the pixel loop, otherwise the
   window array stays in memory and no scalar replacement happens */
#if	defined(__GNUC__)
#define ALWAYS_INLINE	static __inline__ __attribute__((always_inline))
#elif	defined(_MSC_VER)
#define ALWAYS_INLINE	static __forceinline
#else
#define ALWAYS_INLINE	static
#endif

/* iterations of the pixel loops are independent (restrict already says so;
   this is a fallback for the cases where gcc's alias analysis gives up) */
#if	defined(__GNUC__)
#define IVDEP		_Pragma("GCC ivdep")
#else
#define IVDEP
#endif

/*----------------------------------------------------------------------*/
/* index reflection (EDGE_MIRROR), with a final clamp so that a window
   wider than the image cannot read out of bounds */

#define MIRROR(v,n)	do {					\
		if ((v) <  0)  (v) = -(v);			\
		if ((v) >= (n)) (v) = 2*(n)-(v)-2;		\
		if ((v) <  0)  (v) = 0;				\
		if ((v) >= (n)) (v) = (n)-1;			\
	} while (0)

/*----------------------------------------------------------------------*/
/* sorting networks: compare-exchange of two window elements.

   Written as an explicit min/max pair rather than "if (b<a) swap(a,b)".
   The two forms are the same function on integers (they also agree when
   a==b), but the branchless min/max form is what gcc recognises as the
   MIN_EXPR/MAX_EXPR idiom, which it can widen to pminub/pmaxub (8bit) and
   pminuw/pmaxuw (16bit) when the surrounding pixel loop is vectorised.
   The conditional-swap form keeps a data-dependent branch and blocks it. */

#define MSORT(T,a,b)	do {					\
		T lo_=((a)<(b))?(a):(b);			\
		T hi_=((a)<(b))?(b):(a);			\
		(a)=lo_; (b)=hi_;				\
	} while (0)

/* median of 9 (3x3): 19 compare-exchanges.
   ALWAYS_INLINE so that the window array of the caller never has its
   address taken across a call boundary: after inlining, p[] is written and
   read only through constant indices, so gcc's scalar replacement turns it
   into 9 (25) values that the vectoriser can hold in vector lanes. */
#define DEFINE_MED9(NAME,T)						\
ALWAYS_INLINE T NAME(T *p)						\
{									\
	MSORT(T,p[1],p[2]); MSORT(T,p[4],p[5]); MSORT(T,p[7],p[8]);	\
	MSORT(T,p[0],p[1]); MSORT(T,p[3],p[4]); MSORT(T,p[6],p[7]);	\
	MSORT(T,p[1],p[2]); MSORT(T,p[4],p[5]); MSORT(T,p[7],p[8]);	\
	MSORT(T,p[0],p[3]); MSORT(T,p[5],p[8]); MSORT(T,p[4],p[7]);	\
	MSORT(T,p[3],p[6]); MSORT(T,p[1],p[4]); MSORT(T,p[2],p[5]);	\
	MSORT(T,p[4],p[7]); MSORT(T,p[4],p[2]); MSORT(T,p[6],p[4]);	\
	MSORT(T,p[4],p[2]);						\
	return p[4];							\
}

/* median of 25 (5x5): 99 compare-exchanges (ALWAYS_INLINE, see above) */
#define DEFINE_MED25(NAME,T)						\
ALWAYS_INLINE T NAME(T *p)						\
{									\
	MSORT(T,p[0],p[1]);   MSORT(T,p[3],p[4]);   MSORT(T,p[2],p[4]);	\
	MSORT(T,p[2],p[3]);   MSORT(T,p[6],p[7]);   MSORT(T,p[5],p[7]);	\
	MSORT(T,p[5],p[6]);   MSORT(T,p[9],p[10]);  MSORT(T,p[8],p[10]);\
	MSORT(T,p[8],p[9]);   MSORT(T,p[12],p[13]); MSORT(T,p[11],p[13]);\
	MSORT(T,p[11],p[12]); MSORT(T,p[15],p[16]); MSORT(T,p[14],p[16]);\
	MSORT(T,p[14],p[15]); MSORT(T,p[18],p[19]); MSORT(T,p[17],p[19]);\
	MSORT(T,p[17],p[18]); MSORT(T,p[21],p[22]); MSORT(T,p[20],p[22]);\
	MSORT(T,p[20],p[21]); MSORT(T,p[23],p[24]); MSORT(T,p[2],p[5]);	\
	MSORT(T,p[3],p[6]);   MSORT(T,p[0],p[6]);   MSORT(T,p[0],p[3]);	\
	MSORT(T,p[4],p[7]);   MSORT(T,p[1],p[7]);   MSORT(T,p[1],p[4]);	\
	MSORT(T,p[11],p[14]); MSORT(T,p[8],p[14]);  MSORT(T,p[8],p[11]);\
	MSORT(T,p[12],p[15]); MSORT(T,p[9],p[15]);  MSORT(T,p[9],p[12]);\
	MSORT(T,p[13],p[16]); MSORT(T,p[10],p[16]); MSORT(T,p[10],p[13]);\
	MSORT(T,p[20],p[23]); MSORT(T,p[17],p[23]); MSORT(T,p[17],p[20]);\
	MSORT(T,p[21],p[24]); MSORT(T,p[18],p[24]); MSORT(T,p[18],p[21]);\
	MSORT(T,p[19],p[22]); MSORT(T,p[8],p[17]);  MSORT(T,p[9],p[18]);\
	MSORT(T,p[0],p[18]);  MSORT(T,p[0],p[9]);   MSORT(T,p[10],p[19]);\
	MSORT(T,p[1],p[19]);  MSORT(T,p[1],p[10]);  MSORT(T,p[11],p[20]);\
	MSORT(T,p[2],p[20]);  MSORT(T,p[2],p[11]);  MSORT(T,p[12],p[21]);\
	MSORT(T,p[3],p[21]);  MSORT(T,p[3],p[12]);  MSORT(T,p[13],p[22]);\
	MSORT(T,p[4],p[22]);  MSORT(T,p[4],p[13]);  MSORT(T,p[14],p[23]);\
	MSORT(T,p[5],p[23]);  MSORT(T,p[5],p[14]);  MSORT(T,p[15],p[24]);\
	MSORT(T,p[6],p[24]);  MSORT(T,p[6],p[15]);  MSORT(T,p[7],p[16]);\
	MSORT(T,p[7],p[19]);  MSORT(T,p[13],p[21]); MSORT(T,p[15],p[23]);\
	MSORT(T,p[7],p[13]);  MSORT(T,p[7],p[15]);  MSORT(T,p[1],p[9]);	\
	MSORT(T,p[3],p[11]);  MSORT(T,p[5],p[17]);  MSORT(T,p[11],p[17]);\
	MSORT(T,p[9],p[17]);  MSORT(T,p[4],p[10]);  MSORT(T,p[6],p[12]);\
	MSORT(T,p[7],p[14]);  MSORT(T,p[4],p[6]);   MSORT(T,p[4],p[7]);	\
	MSORT(T,p[12],p[14]); MSORT(T,p[10],p[14]); MSORT(T,p[6],p[7]);	\
	MSORT(T,p[10],p[12]); MSORT(T,p[6],p[10]);  MSORT(T,p[6],p[17]);\
	MSORT(T,p[12],p[17]); MSORT(T,p[7],p[17]);  MSORT(T,p[7],p[10]);\
	MSORT(T,p[12],p[18]); MSORT(T,p[7],p[12]);  MSORT(T,p[10],p[18]);\
	MSORT(T,p[12],p[20]); MSORT(T,p[10],p[20]); MSORT(T,p[10],p[12]);\
	return p[12];							\
}

/* insertion sort (k > 5 path) */
#define DEFINE_ISORT(NAME,T)						\
static void NAME(T *arr,int size)					\
{									\
	int	i,j;							\
	T	temp;							\
									\
	for (i = 1; i < size; i++) {					\
	    temp = arr[i];						\
	    j = i-1;							\
	    while (j >= 0 && arr[j] > temp) { arr[j+1] = arr[j]; j--; }	\
	    arr[j+1] = temp;						\
	}								\
}

DEFINE_MED9 (med9_u8 ,unsigned char)
DEFINE_MED9 (med9_u16,unsigned short)
DEFINE_MED25(med25_u8 ,unsigned char)
DEFINE_MED25(med25_u16,unsigned short)
DEFINE_ISORT(isort_u8 ,unsigned char)
DEFINE_ISORT(isort_u16,unsigned short)

/*----------------------------------------------------------------------*/
/* median filter, one output row.  Window is (2*half+1)^2, mirror boundary.
   The OpenMP directive lives in the wrapper below, not in this macro, so
   that no _Pragma() is needed (MSVC /openmp uses the legacy preprocessor). */

#define DEFINE_MEDIAN_ROW(NAME,T,MED9,MED25,ISORT)			\
static void NAME(const T * RESTRICT src,T * RESTRICT dst,		\
		 int width,int height,int half,int y)			\
{									\
	const T	*rows[MEDIAN_KMAX];					\
	T	win[MEDIAN_KMAX*MEDIAN_KMAX];				\
	T	* RESTRICT out=dst+(size_t)y*width;			\
	const int k=2*half+1, n=k*k;					\
	int	xi0=(half<width)?half:width,				\
		xi1=(width-half>xi0)?width-half:xi0,			\
		x,i,j,sy;						\
									\
	for (j = 0; j < k; j++) {					\
	    sy=y+j-half; MIRROR(sy,height);				\
	    rows[j]=src+(size_t)sy*width;				\
	}								\
	/* left border: reflection can trigger */			\
	for (x = 0; x < xi0; x++) {					\
	    int c=0;							\
	    for (j = 0; j < k; j++)					\
	    for (i = 0; i < k; i++) {					\
		int nx=x+i-half; MIRROR(nx,width);			\
		win[c++]=rows[j][nx];					\
	    }								\
	    if	    (k==3) out[x]=MED9(win);				\
	    else if (k==5) out[x]=MED25(win);				\
	    else { ISORT(win,n); out[x]=win[n/2]; }			\
	}								\
	/* Interior: no reflection can trigger here, so the window is read
	   straight through the row base pointers with compile-time constant
	   offsets.  Each x is independent and every access is contiguous in
	   x, which is what lets gcc vectorise these two loops across
	   neighbouring pixels. */					\
	if (k == 3) {							\
	    const T * RESTRICT q0=rows[0];				\
	    const T * RESTRICT q1=rows[1];				\
	    const T * RESTRICT q2=rows[2];				\
	    IVDEP							\
	    for (x = xi0; x < xi1; x++) {				\
		T p[9];							\
		p[0]=q0[x-1]; p[1]=q0[x]; p[2]=q0[x+1];			\
		p[3]=q1[x-1]; p[4]=q1[x]; p[5]=q1[x+1];			\
		p[6]=q2[x-1]; p[7]=q2[x]; p[8]=q2[x+1];			\
		out[x]=MED9(p);						\
	    }								\
	} else if (k == 5) {						\
	    const T * RESTRICT q0=rows[0];				\
	    const T * RESTRICT q1=rows[1];				\
	    const T * RESTRICT q2=rows[2];				\
	    const T * RESTRICT q3=rows[3];				\
	    const T * RESTRICT q4=rows[4];				\
	    IVDEP							\
	    for (x = xi0; x < xi1; x++) {				\
		T p[25];						\
		p[ 0]=q0[x-2]; p[ 1]=q0[x-1]; p[ 2]=q0[x];		\
		p[ 3]=q0[x+1]; p[ 4]=q0[x+2];				\
		p[ 5]=q1[x-2]; p[ 6]=q1[x-1]; p[ 7]=q1[x];		\
		p[ 8]=q1[x+1]; p[ 9]=q1[x+2];				\
		p[10]=q2[x-2]; p[11]=q2[x-1]; p[12]=q2[x];		\
		p[13]=q2[x+1]; p[14]=q2[x+2];				\
		p[15]=q3[x-2]; p[16]=q3[x-1]; p[17]=q3[x];		\
		p[18]=q3[x+1]; p[19]=q3[x+2];				\
		p[20]=q4[x-2]; p[21]=q4[x-1]; p[22]=q4[x];		\
		p[23]=q4[x+1]; p[24]=q4[x+2];				\
		out[x]=MED25(p);					\
	    }								\
	} else {							\
	    for (x = xi0; x < xi1; x++) {				\
		int c=0;						\
		for (j = 0; j < k; j++) {				\
		    const T *r=rows[j]+x-half;				\
		    for (i = 0; i < k; i++) win[c++]=r[i];		\
		}							\
		ISORT(win,n); out[x]=win[n/2];				\
	    }								\
	}								\
	/* right border */						\
	for (x = xi1; x < width; x++) {					\
	    int c=0;							\
	    for (j = 0; j < k; j++)					\
	    for (i = 0; i < k; i++) {					\
		int nx=x+i-half; MIRROR(nx,width);			\
		win[c++]=rows[j][nx];					\
	    }								\
	    if	    (k==3) out[x]=MED9(win);				\
	    else if (k==5) out[x]=MED25(win);				\
	    else { ISORT(win,n); out[x]=win[n/2]; }			\
	}								\
}

DEFINE_MEDIAN_ROW(median_row_u8 ,unsigned char ,med9_u8 ,med25_u8 ,isort_u8)
DEFINE_MEDIAN_ROW(median_row_u16,unsigned short,med9_u16,med25_u16,isort_u16)

static void	median_u8(const unsigned char * RESTRICT src,
			  unsigned char * RESTRICT dst,
			  int width,int height,int half)
{
	int	y;

#pragma omp parallel for schedule(static)
	for (y = 0; y < height; y++)
	    median_row_u8(src,dst,width,height,half,y);
}

static void	median_u16(const unsigned short * RESTRICT src,
			   unsigned short * RESTRICT dst,
			   int width,int height,int half)
{
	int	y;

#pragma omp parallel for schedule(static)
	for (y = 0; y < height; y++)
	    median_row_u16(src,dst,width,height,half,y);
}

/*----------------------------------------------------------------------*/
/* gaussian kernel */

static int	calculate_kernel_size(double sigma)
{
	int	size;

	if (sigma <= 0.0) return 0;
	size=(int)ceil(3.0*sigma)*2+1;
	if (size%2 == 0) size++;
	if (size > GAUSS_KMAX) size=GAUSS_KMAX;
	return size;
}

static void	create_gaussian_kernel_1d(double *kernel,int size,double sigma)
{
	int	center=size/2,i;
	double	sum=0.0,
		two_sigma_sq=2.0*sigma*sigma,
		norm_factor=1.0/sqrt(2.0*3.14159265358979323846*sigma*sigma);

	for (i = 0; i < size; i++) {
	    int distance=i-center;

	    kernel[i]=norm_factor*exp(-(distance*distance)/two_sigma_sq);
	    sum+=kernel[i];
	}
	for (i = 0; i < size; i++) kernel[i]/=sum;
}

/* horizontal pass; each pixel sums its taps in ascending order */
static void	gauss_h(const double * RESTRICT src,double * RESTRICT dst,
			int width,int height,const double * RESTRICT ker,int ks)
{
	const int	half=ks/2;
	int		xi0=(half<width)?half:width,
			xi1=(width-half>xi0)?width-half:xi0,
			y;

#pragma omp parallel for schedule(static)
	for (y = 0; y < height; y++) {
	    const double	* RESTRICT sr=src+(size_t)y*width;
	    double		* RESTRICT dr=dst+(size_t)y*width;
	    int			x,i;
	    double		s;

	    for (x = 0; x < xi0; x++) {
		s=0.0;
		for (i = 0; i < ks; i++) {
		    int sx=x+i-half; MIRROR(sx,width);
		    s+=sr[sx]*ker[i];
		}
		dr[x]=s;
	    }
	    IVDEP
	    for (x = xi0; x < xi1; x++) {
		const double *p=sr+x-half;

		s=0.0;
		for (i = 0; i < ks; i++) s+=p[i]*ker[i];
		dr[x]=s;
	    }
	    for (x = xi1; x < width; x++) {
		s=0.0;
		for (i = 0; i < ks; i++) {
		    int sx=x+i-half; MIRROR(sx,width);
		    s+=sr[sx]*ker[i];
		}
		dr[x]=s;
	    }
	}
}

/* vertical pass; row-blocked instead of column-blocked (same arithmetic) */
static void	gauss_v(const double * RESTRICT src,double * RESTRICT dst,
			int width,int height,const double * RESTRICT ker,int ks)
{
	const int	half=ks/2;
	int		y;

#pragma omp parallel for schedule(static)
	for (y = 0; y < height; y++) {
	    const double	*rows[GAUSS_KMAX];
	    double		* RESTRICT dr=dst+(size_t)y*width;
	    int			x,i,sy;
	    double		s;

	    for (i = 0; i < ks; i++) {
		sy=y+i-half; MIRROR(sy,height);
		rows[i]=src+(size_t)sy*width;
	    }
	    IVDEP
	    for (x = 0; x < width; x++) {
		s=0.0;
		for (i = 0; i < ks; i++) s+=rows[i][x]*ker[i];
		dr[x]=s;
	    }
	}
}

/*----------------------------------------------------------------------*/
/* integer <-> double conversion */

/* The pixel count is passed as a plain int: MSVC's /openmp implements
   OpenMP 2.0, whose parallel-for index must be a signed int.  An 8/16-bit
   image of more than 2G pixels cannot be handled by the rest of the
   program either. */
static void	to_double(const void *data,double * RESTRICT out,int n,int is_16bit)
{
	int	i;

	if (is_16bit) {
	    const unsigned short * RESTRICT src=(const unsigned short *)data;

#pragma omp parallel for schedule(static)
	    for (i = 0; i < n; i++) out[i]=(double)src[i];
	} else {
	    const unsigned char * RESTRICT src=(const unsigned char *)data;

#pragma omp parallel for schedule(static)
	    for (i = 0; i < n; i++) out[i]=(double)src[i];
	}
}

static void	from_double(const double * RESTRICT in,void *data,int n,int is_16bit)
{
	int	i;

	if (is_16bit) {
	    unsigned short * RESTRICT dst=(unsigned short *)data;

#pragma omp parallel for schedule(static)
	    for (i = 0; i < n; i++) {
		double val=in[i]+0.5;

		if (val <     0.0) val=0.0;
		if (val > 65535.0) val=65535.0;
		dst[i]=(unsigned short)val;
	    }
	} else {
	    unsigned char * RESTRICT dst=(unsigned char *)data;

#pragma omp parallel for schedule(static)
	    for (i = 0; i < n; i++) {
		double val=in[i]+0.5;

		if (val <   0.0) val=0.0;
		if (val > 255.0) val=255.0;
		dst[i]=(unsigned char)val;
	    }
	}
}

/*----------------------------------------------------------------------*/

static int	apply_median(void *data,int width,int height,int kernel_size,
			     int is_16bit)
{
	const size_t	npix=(size_t)width*height;
	const int	half=kernel_size/2;
	void		*tmp;

	if (half <= 0) return 0;			/* k=1: identity */

	tmp=malloc(npix*(is_16bit?sizeof(unsigned short):sizeof(unsigned char)));
	if (tmp == NULL) {
	    fprintf(stderr,"Error: Cannot allocate memory for median filter\n");
	    return 1;
	}
	if (is_16bit)
	    median_u16((const unsigned short *)data,(unsigned short *)tmp,
		       width,height,half);
	else
	    median_u8 ((const unsigned char  *)data,(unsigned char  *)tmp,
		       width,height,half);

	memcpy(data,tmp,npix*(is_16bit?sizeof(unsigned short):sizeof(unsigned char)));
	free(tmp);
	return 0;
}

static int	apply_gaussian(void *data,int width,int height,double sigma,
			       int is_16bit)
{
	const size_t	npix=(size_t)width*height;
	int		ks=calculate_kernel_size(sigma);
	double		*f,*t,ker[GAUSS_KMAX];

	if (ks <= 0) return 0;

	f=(double *)malloc(npix*sizeof(double));
	t=(double *)malloc(npix*sizeof(double));
	if (f == NULL || t == NULL) {
	    fprintf(stderr,"Error: Cannot allocate float buffer\n");
	    free(f); free(t);
	    return 1;
	}
	create_gaussian_kernel_1d(ker,ks,sigma);

	to_double(data,f,(int)npix,is_16bit);
	gauss_h(f,t,width,height,ker,ks);
	gauss_v(t,f,width,height,ker,ks);
	from_double(f,data,(int)npix,is_16bit);

	free(t); free(f);
	return 0;
}

/*======================================================================*/
/* ここまで tif_mgf.c 由来                                               */
/*======================================================================*/

#define IN_FMT		"a%06d.tif"
#define OUT_NAME	"a.his"
#define HIS_HDR		64

typedef struct {
	int	w,h,bps;	/* bps: 8 or 16 */
	uint8_t	*pix;		/* w*h*(bps/8) bytes, ネイティブ(リトル)エンディアン */
} Img;

static int	g_be;		/* 読んでいる TIFF がビッグエンディアンか */

static uint16_t	g16(const uint8_t *p)
{
	return g_be ? (uint16_t)((p[0]<<8)|p[1]) : (uint16_t)(p[0]|(p[1]<<8));
}
static uint32_t	g32(const uint8_t *p)
{
	return g_be ? ((uint32_t)p[0]<<24)|((uint32_t)p[1]<<16)|((uint32_t)p[2]<<8)|p[3]
		    : p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);
}
static void	p16(uint8_t *p,uint16_t v) { p[0]=v&0xff; p[1]=v>>8; }
static void	p32(uint8_t *p,uint32_t v)
{
	p[0]=v&0xff; p[1]=(v>>8)&0xff; p[2]=(v>>16)&0xff; p[3]=v>>24;
}

static int	file_exists(const char *fn)
{
	FILE	*fp=fopen(fn,"rb");

	if (fp == NULL) return 0;
	fclose(fp);
	return 1;
}

static uint8_t	*read_all(const char *fn,size_t *size)
{
	FILE	*fp=fopen(fn,"rb");
	uint8_t	*buf;
	long	n;

	if (fp == NULL) return NULL;
	if (fseek(fp,0,SEEK_END) != 0 || (n=ftell(fp)) <= 0) { fclose(fp); return NULL; }
	rewind(fp);
	buf=(uint8_t *)malloc((size_t)n);
	if (buf != NULL && fread(buf,1,(size_t)n,fp) != (size_t)n) { free(buf); buf=NULL; }
	fclose(fp);
	*size=(size_t)n;
	return buf;
}

/* IFD エントリの i 番目の値 (SHORT / LONG のみ) */
static uint32_t	ent_val(const uint8_t *f,size_t fsz,const uint8_t *e,uint32_t i,int *err)
{
	uint16_t	type=g16(e+2);
	uint32_t	cnt=g32(e+4);
	int		sz=(type == 3) ? 2 : (type == 4) ? 4 : 0;
	const uint8_t	*p;

	if (sz == 0 || i >= cnt) { *err=1; return 0; }
	if ((uint64_t)cnt*sz <= 4) p=e+8;
	else {
	    uint32_t off=g32(e+8);

	    if ((uint64_t)off+(uint64_t)cnt*sz > fsz) { *err=1; return 0; }
	    p=f+off;
	}
	return (sz == 2) ? g16(p+2*i) : g32(p+4*i);
}

static int	read_tiff(const char *fn,Img *img)
{
	size_t		fsz,need,pos=0;
	uint8_t		*f=read_all(fn,&fsz);
	const uint8_t	*e_off=NULL,*e_cnt=NULL;
	uint32_t	ifd,nstrip=0,i;
	int		n,k,err=0,comp=1,spp=1,sfmt=1,tiled=0,bpb;

	memset(img,0,sizeof(*img));
	if (f == NULL) { fprintf(stderr,"%s: 読み込めません\n",fn); return -1; }
	if (fsz < 8 || !((f[0]=='I' && f[1]=='I') || (f[0]=='M' && f[1]=='M'))) goto bad;
	g_be=(f[0] == 'M');
	if (g16(f+2) != 42) goto bad;		/* BigTIFF(43) は非対応 */
	ifd=g32(f+4);
	if ((uint64_t)ifd+2 > fsz) goto bad;
	n=g16(f+ifd);
	if ((uint64_t)ifd+2+12ull*n > fsz) goto bad;

	for (k = 0; k < n; k++) {
	    const uint8_t *e=f+ifd+2+12*k;

	    switch (g16(e)) {
	    case 256: img->w  =(int)ent_val(f,fsz,e,0,&err); break;
	    case 257: img->h  =(int)ent_val(f,fsz,e,0,&err); break;
	    case 258: img->bps=(int)ent_val(f,fsz,e,0,&err); break;
	    case 259: comp    =(int)ent_val(f,fsz,e,0,&err); break;
	    case 273: e_off=e; nstrip=g32(e+4); break;
	    case 277: spp     =(int)ent_val(f,fsz,e,0,&err); break;
	    case 279: e_cnt=e; break;
	    case 322: tiled=1; break;
	    case 339: sfmt    =(int)ent_val(f,fsz,e,0,&err); break;
	    }
	}
	if (err || e_off == NULL || e_cnt == NULL || img->w <= 0 || img->h <= 0) goto bad;
	if (comp != 1) { fprintf(stderr,"%s: 圧縮 TIFF (Compression=%d) には非対応です\n",fn,comp); goto fail; }
	if (tiled)     { fprintf(stderr,"%s: タイル形式 TIFF には非対応です\n",fn); goto fail; }
	if (spp != 1)  { fprintf(stderr,"%s: モノクロ 1 サンプル/画素のみ対応です\n",fn); goto fail; }
	if (sfmt != 1 || (img->bps != 8 && img->bps != 16)) {
	    fprintf(stderr,"%s: 8/16bit 符号なし整数のみ対応です (bps=%d, SampleFormat=%d)\n",fn,img->bps,sfmt);
	    goto fail;
	}
	if (img->w > 65535 || img->h > 65535) { fprintf(stderr,"%s: 画像が大きすぎます\n",fn); goto fail; }

	bpb=img->bps/8;
	need=(size_t)img->w*img->h*bpb;
	if ((img->pix=(uint8_t *)malloc(need)) == NULL) { fprintf(stderr,"メモリ不足\n"); goto fail; }

	for (i = 0; i < nstrip && pos < need; i++) {
	    uint32_t off=ent_val(f,fsz,e_off,i,&err),
		     len=ent_val(f,fsz,e_cnt,i,&err);

	    if (err || (uint64_t)off+len > fsz) goto bad;
	    if (len > need-pos) len=(uint32_t)(need-pos);
	    memcpy(img->pix+pos,f+off,len);
	    pos+=len;
	}
	if (pos != need) goto bad;

	if (bpb == 2 && g_be) {			/* フィルタも HIS もリトルエンディアン前提 */
	    size_t j;

	    for (j = 0; j < need; j += 2) {
		uint8_t t=img->pix[j]; img->pix[j]=img->pix[j+1]; img->pix[j+1]=t;
	    }
	}
	free(f);
	return 0;

bad:
	fprintf(stderr,"%s: TIFF の構造が不正か、非対応の形式です\n",fn);
fail:
	free(f);
	free(img->pix);
	img->pix=NULL;
	return -1;
}

int	main(int argc,char *argv[])
{
	char	fn[64];
	int	start=-1,nimg,i,w0=0,h0=0,bps0=0,
		median_kernel_size=0;
	double	gaussian_sigma=0.0;
	FILE	*out,*flog;

	if (argc > 3) {
	    fprintf(stderr,"Usage: %s [median_kernel_size] [gaussian_sigma]\n",argv[0]);
	    fprintf(stderr,"  median_kernel_size: 0-%d (default: 0 = skip; 3 and 5 take the fast path)\n",MEDIAN_KMAX);
	    fprintf(stderr,"  gaussian_sigma: 0.0-100.0 (default: 0.0 = skip)\n");
	    fprintf(stderr,"  packs a%%06d.tif in the current directory into " OUT_NAME "\n");
	    return 1;
	}
	if (argc > 1) {
	    median_kernel_size=atoi(argv[1]);
	    if (median_kernel_size < 0 || median_kernel_size > MEDIAN_KMAX) {
		fprintf(stderr,"Error: Invalid median kernel size (must be 0-%d)\n",MEDIAN_KMAX);
		return 1;
	    }
	}
	if (argc > 2) {
	    gaussian_sigma=atof(argv[2]);
	    if (gaussian_sigma < 0.0 || gaussian_sigma > 100.0) {
		fprintf(stderr,"Error: Invalid gaussian sigma (must be 0.0-100.0)\n");
		return 1;
	    }
	}

	for (i = 0; i <= 1; i++) {
	    snprintf(fn,sizeof fn,IN_FMT,i);
	    if (file_exists(fn)) { start=i; break; }
	}
	if (start < 0) { fprintf(stderr,"a000000.tif / a000001.tif が見つかりません\n"); return 1; }
	for (nimg = 0; nimg < 1000000-start; nimg++) {
	    snprintf(fn,sizeof fn,IN_FMT,start+nimg);
	    if (!file_exists(fn)) break;
	}

	if (file_exists(OUT_NAME)) { fprintf(stderr,"%s が既にあります (上書きしません)\n",OUT_NAME); return 1; }
	if ((out=fopen(OUT_NAME,"wb")) == NULL) { fprintf(stderr,"%s を作成できません\n",OUT_NAME); return 1; }

	fprintf(stderr,"%d frames, median %d, sigma %g\n",nimg,median_kernel_size,gaussian_sigma);

	for (i = 0; i < nimg; i++) {
	    Img		img;
	    uint8_t	hdr[HIS_HDR];
	    size_t	nbytes;
	    int		is_16bit;

	    snprintf(fn,sizeof fn,IN_FMT,start+i);
	    if (read_tiff(fn,&img) != 0) goto abort;

	    if (i == 0) {
		w0=img.w; h0=img.h; bps0=img.bps;
		if (bps0 != 16)
		    fprintf(stderr,"警告: 8bit 画像です。rhis (his_rec/his_srec 等) は 16bit の his しか読めません\n");
	    } else if (img.w != w0 || img.h != h0 || img.bps != bps0) {
		fprintf(stderr,"%s: 画像サイズ/ビット深度が先頭フレームと異なります (%dx%d %dbit)\n",
			fn,img.w,img.h,img.bps);
		free(img.pix);
		goto abort;
	    }

	    /* tif_mgf と同じ順序: median → gaussian */
	    is_16bit=(img.bps == 16);
	    if ((median_kernel_size > 0 &&
		 apply_median(img.pix,img.w,img.h,median_kernel_size,is_16bit) != 0) ||
		(gaussian_sigma > 0.0 &&
		 apply_gaussian(img.pix,img.w,img.h,gaussian_sigma,is_16bit) != 0)) {
		free(img.pix);
		goto abort;
	    }

	    memset(hdr,0,sizeof hdr);
	    hdr[0]='I'; hdr[1]='M';
	    p16(hdr+2,0);				/* comment length */
	    p16(hdr+4,(uint16_t)img.w);
	    p16(hdr+6,(uint16_t)img.h);
	    p16(hdr+8,0);				/* x offset */
	    p16(hdr+10,0);				/* y offset */
	    p16(hdr+12,is_16bit ? 2 : 0);		/* type: 0=8bit, 2=16bit */
	    if (i == 0) p32(hdr+14,(uint32_t)nimg);	/* number of frames */

	    nbytes=(size_t)img.w*img.h*(img.bps/8);
	    if (fwrite(hdr,1,HIS_HDR,out) != HIS_HDR ||
		fwrite(img.pix,1,nbytes,out) != nbytes) {
		fprintf(stderr,"書き込みエラー (ディスク容量を確認してください)\n");
		free(img.pix);
		goto abort;
	    }
	    free(img.pix);

	    if ((i+1)%100 == 0 || i+1 == nimg) fprintf(stderr,"\r%d / %d",i+1,nimg);
	}
	fprintf(stderr,"\n");
	if (fclose(out) != 0) { fprintf(stderr,"書き込みエラー\n"); remove(OUT_NAME); return 1; }

	printf("a%06d.tif - a%06d.tif (%d frames, %dx%d, %dbit) -> %s\n",
	       start,start+nimg-1,nimg,w0,h0,bps0,OUT_NAME);

	/* append to log file (tif_mgf と同じ形式) */
	if ((flog=fopen("cmd-hst.log","a")) != NULL) {
	    for (i = 0; i < argc; ++i) fprintf(flog,"%s ",argv[i]);
	    fprintf(flog,"\t");
	    fprintf(flog,"   %% median_kernel_size %d  gaussian_sigma %g  frames %d\n",
		    median_kernel_size,gaussian_sigma,nimg);
	    fclose(flog);
	}
	return 0;

abort:
	fclose(out);
	remove(OUT_NAME);
	return 1;
}