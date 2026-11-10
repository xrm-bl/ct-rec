/*
 * his_mgf.c - HiPic シーケンスファイル (.his) の各フレームに median filter と
 *             gaussian filter をかけて、別の .his に出力する。
 *
 *   his_mgf <input.his> <output.his> <median_kernel_size> <gaussian_sigma>
 *     median_kernel_size : 0-25      (0, 1 はかけない)
 *     gaussian_sigma     : 0.0-100.0 (0.0 はかけない)
 *   引数は 4 つすべて必須。足りない・多い・数値として読めない・範囲外の
 *   場合はエラー終了する。
 *
 * フィルタ部 (median → gaussian、ミラー境界、double 演算、+0.5 で丸め) は
 * ct-rec の tif_mgf.c のコードをそのまま使っているので、各フレームの結果は
 * tif_mgf で 1 枚ずつ処理した場合とビット単位で一致する。
 *
 * his の扱い (rhis.c と同じ解釈):
 *   各フレーム = 64 バイトヘッダ + コメント (ヘッダ offset 2 の長さ) + 画素
 *   幅 offset 4、高さ offset 6、type offset 12 (0=8bit, 2=16bit)、
 *   総フレーム数は先頭フレームの offset 14 (32bit)。リトルエンディアン。
 *   各フレームのヘッダとコメントはそのままコピーし、画素だけを置き換える。
 *   したがって出力は入力と同じサイズ・同じレイアウトになる。
 *   フレームごとにヘッダを読むので、フレーム間でサイズが違っても処理できる。
 *
 * 出力ファイルが既にある場合 (入力と同じ名前を含む) は上書きせずに終了。
 * 途中でエラーになった場合は書きかけの出力を削除する。
 * 実行内容は tif_mgf と同じ形式で cmd-hst.log に追記する。
 *
 * Build (Linux):   gcc -O3 -march=native -fno-math-errno -fopenmp -o his_mgf his_mgf.c -lm
 * Build (Windows): cl /O2 /openmp his_mgf.c
 *   tif_mgf と同様、-ffast-math / -Ofast は付けないこと (gaussian の加算順が
 *   変わり、tif_mgf と結果が一致しなくなる)。
 */
#define _CRT_SECURE_NO_WARNINGS
#define _FILE_OFFSET_BITS 64
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
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

#define HIS_HDR		64

static unsigned	get16(const unsigned char *p) { return p[0]|(p[1]<<8); }
static uint32_t	get32(const unsigned char *p)
{
	return p[0]|((uint32_t)p[1]<<8)|((uint32_t)p[2]<<16)|((uint32_t)p[3]<<24);
}

/* 整数として完全に読めるか */
static int	parse_int(const char *s,int *v)
{
	char	*end;
	long	x;

	errno=0;
	x=strtol(s,&end,10);
	if (end == s || *end != '\0' || errno != 0) return 0;
	*v=(int)x;
	return 1;
}

/* 実数として完全に読めるか */
static int	parse_double(const char *s,double *v)
{
	char	*end;
	double	x;

	errno=0;
	x=strtod(s,&end);
	if (end == s || *end != '\0' || errno != 0 || x != x) return 0;
	*v=x;
	return 1;
}

static int	file_exists(const char *fn)
{
	FILE	*fp=fopen(fn,"rb");

	if (fp == NULL) return 0;
	fclose(fp);
	return 1;
}

static void	usage(const char *prog)
{
	fprintf(stderr,"Usage: %s <input.his> <output.his> <median_kernel_size> <gaussian_sigma>\n",prog);
	fprintf(stderr,"  median_kernel_size: 0-%d (0 or 1 to skip; 3 and 5 take the fast path)\n",MEDIAN_KMAX);
	fprintf(stderr,"  gaussian_sigma: 0.0-100.0 (0.0 to skip)\n");
}

int	main(int argc,char *argv[])
{
	const char	*in_name,*out_name;
	int		median_kernel_size,i;
	double		gaussian_sigma;
	FILE		*in,*out,*flog;
	unsigned char	hdr[HIS_HDR],*comment=NULL,*pix=NULL;
	size_t		comment_cap=0,pix_cap=0;
	uint32_t	nframes,k;
	unsigned	w0=0,h0=0,type0=0;
	int		mixed=0;

	if (argc != 5) {
	    fprintf(stderr,"Error: 引数が %d 個です (4 個必要)\n",argc-1);
	    usage(argv[0]);
	    return 1;
	}
	in_name =argv[1];
	out_name=argv[2];
	if (!parse_int(argv[3],&median_kernel_size) ||
	    median_kernel_size < 0 || median_kernel_size > MEDIAN_KMAX) {
	    fprintf(stderr,"Error: Invalid median kernel size \"%s\" (must be 0-%d)\n",argv[3],MEDIAN_KMAX);
	    return 1;
	}
	if (!parse_double(argv[4],&gaussian_sigma) ||
	    gaussian_sigma < 0.0 || gaussian_sigma > 100.0) {
	    fprintf(stderr,"Error: Invalid gaussian sigma \"%s\" (must be 0.0-100.0)\n",argv[4]);
	    return 1;
	}

	if ((in=fopen(in_name,"rb")) == NULL) {
	    fprintf(stderr,"Error: %s を開けません\n",in_name);
	    return 1;
	}
	if (fread(hdr,1,HIS_HDR,in) != HIS_HDR || hdr[0] != 'I' || hdr[1] != 'M') {
	    fprintf(stderr,"Error: %s は his ファイルではありません\n",in_name);
	    fclose(in);
	    return 1;
	}
	nframes=get32(hdr+14);
	if (nframes == 0) {
	    fprintf(stderr,"Error: %s のフレーム数 (offset 14) が 0 です\n",in_name);
	    fclose(in);
	    return 1;
	}

	if (file_exists(out_name)) {
	    fprintf(stderr,"Error: %s が既にあります (上書きしません)\n",out_name);
	    fclose(in);
	    return 1;
	}
	if ((out=fopen(out_name,"wb")) == NULL) {
	    fprintf(stderr,"Error: %s を作成できません\n",out_name);
	    fclose(in);
	    return 1;
	}

	fprintf(stderr,"%s: %u frames, median %d, sigma %g\n",
		in_name,nframes,median_kernel_size,gaussian_sigma);

	for (k = 0; k < nframes; k++) {
	    unsigned	cl,w,h,type;
	    size_t	nbytes;
	    int		is_16bit;

	    /* 先頭フレームのヘッダは読み込み済み */
	    if (k > 0 && (fread(hdr,1,HIS_HDR,in) != HIS_HDR || hdr[0] != 'I' || hdr[1] != 'M')) {
		fprintf(stderr,"\nError: frame %u のヘッダが読めません (ファイルが途中で切れている? %u frames のはず)\n",
			k+1,nframes);
		goto abort;
	    }
	    cl=get16(hdr+2); w=get16(hdr+4); h=get16(hdr+6); type=get16(hdr+12);
	    if (type != 0 && type != 2) {
		fprintf(stderr,"\nError: frame %u: type %u には非対応です (0=8bit, 2=16bit のみ)\n",k+1,type);
		goto abort;
	    }
	    if (w == 0 || h == 0) {
		fprintf(stderr,"\nError: frame %u: 画像サイズが 0 です\n",k+1);
		goto abort;
	    }
	    if (k == 0) { w0=w; h0=h; type0=type; }
	    else if (w != w0 || h != h0 || type != type0) mixed=1;

	    is_16bit=(type == 2);
	    nbytes=(size_t)w*h*(is_16bit ? 2 : 1);

	    if (cl > comment_cap) {
		unsigned char *p=(unsigned char *)realloc(comment,cl);

		if (p == NULL) { fprintf(stderr,"\nError: メモリ不足\n"); goto abort; }
		comment=p; comment_cap=cl;
	    }
	    if (nbytes > pix_cap) {
		unsigned char *p=(unsigned char *)realloc(pix,nbytes);

		if (p == NULL) { fprintf(stderr,"\nError: メモリ不足\n"); goto abort; }
		pix=p; pix_cap=nbytes;
	    }
	    if ((cl > 0 && fread(comment,1,cl,in) != cl) || fread(pix,1,nbytes,in) != nbytes) {
		fprintf(stderr,"\nError: frame %u のデータが読めません (ファイルが途中で切れている? %u frames のはず)\n",
			k+1,nframes);
		goto abort;
	    }

	    /* tif_mgf と同じ順序: median → gaussian (画素はリトルエンディアンのまま扱う) */
	    if ((median_kernel_size > 0 &&
		 apply_median(pix,(int)w,(int)h,median_kernel_size,is_16bit) != 0) ||
		(gaussian_sigma > 0.0 &&
		 apply_gaussian(pix,(int)w,(int)h,gaussian_sigma,is_16bit) != 0))
		goto abort;

	    if (fwrite(hdr,1,HIS_HDR,out) != HIS_HDR ||
		(cl > 0 && fwrite(comment,1,cl,out) != cl) ||
		fwrite(pix,1,nbytes,out) != nbytes) {
		fprintf(stderr,"\nError: 書き込みエラー (ディスク容量を確認してください)\n");
		goto abort;
	    }

	    if ((k+1)%100 == 0 || k+1 == nframes) fprintf(stderr,"\r%u / %u",k+1,nframes);
	}
	fprintf(stderr,"\n");

	if (fgetc(in) != EOF)
	    fprintf(stderr,"警告: %s の %u frames の後ろに余分なデータがあります (出力にはコピーしていません)\n",
		    in_name,nframes);
	if (mixed)
	    fprintf(stderr,"注意: フレームごとにサイズまたは type が異なる his でした\n");

	fclose(in);
	free(comment); free(pix);
	if (fclose(out) != 0) {
	    fprintf(stderr,"Error: 書き込みエラー\n");
	    remove(out_name);
	    return 1;
	}

	printf("%s (%u frames, %ux%u, type %u) -> %s  [median %d, sigma %g]\n",
	       in_name,nframes,w0,h0,type0,out_name,median_kernel_size,gaussian_sigma);

	/* append to log file (tif_mgf と同じ形式) */
	if ((flog=fopen("cmd-hst.log","a")) != NULL) {
	    for (i = 0; i < argc; ++i) fprintf(flog,"%s ",argv[i]);
	    fprintf(flog,"\t");
	    fprintf(flog,"   %% median_kernel_size %d  gaussian_sigma %g  frames %u\n",
		    median_kernel_size,gaussian_sigma,nframes);
	    fclose(flog);
	}
	return 0;

abort:
	fclose(in);
	fclose(out);
	remove(out_name);
	free(comment); free(pix);
	return 1;
}