/*
 * tif_mgf5.c - 2D filter (median + gaussian) for a single TIFF,
 *              specialised for small median kernels (3x3 / 5x5).
 *
 * Drop-in replacement for tif_mgf.c:
 *   tif_mgf5 <input_file> <output_file> [median_kernel_size] [gaussian_sigma]
 * Same CLI, same edge mode (mirror), same order (median first, then the
 * separable gaussian), same double-precision arithmetic in the same order,
 * so the output pixels are identical to tif_mgf.c.
 *
 * Why it is faster
 *   1. Median: branch-free sorting networks - 19 compare-exchanges for 3x3
 *      (Smith 1996) and 99 for 5x5 - instead of an insertion sort over the
 *      whole window (~20 / ~150 element moves on average, with a data
 *      dependent branch per move).  The median of a window is by definition
 *      the same value, so the result does not change.
 *   2. Median: the row is split into a 'half'-wide border and an interior.
 *      In the interior (more than 99% of the pixels for k <= 5) the mirror
 *      test cannot trigger, so the four boundary comparisons per window
 *      sample disappear and the window is gathered through row pointers.
 *   3. Median: the malloc/free of the window buffer that tif_mgf.c does
 *      once per image row is gone; the window lives in local variables.
 *   4. Gaussian, vertical pass: tif_mgf.c walks columns (x outer, y inner),
 *      so every single sample touches a different cache line.  Here the
 *      pass is row-blocked (y outer, x inner) with the mirrored source row
 *      pointers computed once per output row.  Each output pixel still
 *      accumulates the same products in the same order, so the arithmetic
 *      is unchanged, but the memory traffic drops by more than an order of
 *      magnitude on large images.
 *   5. Gaussian, both passes: same interior/border split as in 2., and the
 *      integer <-> double conversions are OpenMP-parallel.
 *
 * Median kernel sizes above 5 still work: they fall back to the insertion
 * sort of tif_mgf.c (items 2-5 above still apply), so nothing is lost by
 * using this program everywhere.
 *
 * Even kernel sizes: tif_mgf.c allocates k*k window entries but gathers
 * (2*(k/2)+1)^2 of them, i.e. it overruns its window buffer for k=2 and
 * k=4.  Here an even k is treated as the odd window it actually gathers
 * (k=2 -> 3x3, k=4 -> 5x5), which is what tif_mgf_g.cu does as well.
 *
 * Build (Windows):
 *   cl /DWINDOWS /O2 /openmp /Fetif_mgf5.exe tif_mgf5.c libtiff.lib jpeg.lib lzma.lib zs.lib
 *
 * Build (Linux, recommended - lets gcc vectorise the median networks):
 *   gcc -O3 -march=native -fno-math-errno -o tif_mgf5 tif_mgf5.c -ltiff -lm
 *   add -fopenmp for the multi-threaded build.
 *
 *   Do NOT add -ffast-math (or -Ofast, or -fassociative-math): it would let
 *   gcc reassociate the gaussian tap sum and the output would no longer
 *   match tif_mgf.  -fno-math-errno is safe here (no libm call in the
 *   filter loops; it only frees the compiler from the errno side effect).
 *   -march=native selects pminub/pmaxub / pminuw/pmaxuw; without it the
 *   baseline x86-64 SSE2 forms are used, which are still vectorised.
 *
 *   To inspect what was vectorised:
 *     gcc -O3 -march=native -fno-math-errno -fopt-info-vec-optimized \
 *         -fopt-info-vec-missed -c tif_mgf5.c
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#ifdef _WIN32
#include "tiffio.h"
#else
#include <tiffio.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#define MEDIAN_KMAX	25		/* same limit as tif_mgf.c */
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
/* index reflection (identical to tif_mgf.c EDGE_MIRROR, with a final
   clamp so that a window wider than the image cannot read out of bounds) */

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

/* insertion sort, identical to sort_array_*() of tif_mgf.c (k > 5 path) */
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
/* gaussian kernel (identical math to tif_mgf.c) */

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

/* horizontal pass; per-pixel accumulation order is that of tif_mgf.c */
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
/* integer <-> double conversion (rounding identical to tif_mgf.c) */

/* The pixel count is passed as a plain int: MSVC's /openmp implements
   OpenMP 2.0, whose parallel-for index must be a signed int.  An 8/16-bit
   image of more than 2G pixels cannot be handled by the rest of the
   program (nor by tif_mgf.c) either. */
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

/*----------------------------------------------------------------------*/

static int	process_tiff_file(const char *input_file,const char *output_file,
				  int median_size,double sigma)
{
	TIFF		*in_tiff,*out_tiff;
	uint32_t	width=0,height=0,row;
	uint16_t	bits_per_sample=0;
	void		*data;
	tsize_t		scanline_size;
	char		*desc=NULL;
	int		is_16bit,have_desc;

	if ((in_tiff=TIFFOpen(input_file,"r")) == NULL) {
	    fprintf(stderr,"Error: Cannot open input TIFF file\n");
	    return 1;
	}
	TIFFGetField(in_tiff,TIFFTAG_IMAGEWIDTH,&width);
	TIFFGetField(in_tiff,TIFFTAG_IMAGELENGTH,&height);
	TIFFGetField(in_tiff,TIFFTAG_BITSPERSAMPLE,&bits_per_sample);
	have_desc=TIFFGetField(in_tiff,TIFFTAG_IMAGEDESCRIPTION,&desc);

	if (bits_per_sample != 8 && bits_per_sample != 16) {
	    fprintf(stderr,"Error: Only 8-bit and 16-bit images are supported\n");
	    TIFFClose(in_tiff);
	    return 1;
	}
	is_16bit=(bits_per_sample == 16);
	scanline_size=TIFFScanlineSize(in_tiff);

	if ((data=malloc((size_t)height*scanline_size)) == NULL) {
	    fprintf(stderr,"Error: Cannot allocate memory\n");
	    TIFFClose(in_tiff);
	    return 1;
	}
	for (row = 0; row < height; row++)
	    if (TIFFReadScanline(in_tiff,(char *)data+(size_t)row*scanline_size,
				 row,0) < 0) {
		fprintf(stderr,"Error: Cannot read scanline %u\n",row);
		free(data); TIFFClose(in_tiff);
		return 1;
	    }

	if (have_desc && desc != NULL) {
	    char *copy=(char *)malloc(strlen(desc)+1);

	    if (copy != NULL) { strcpy(copy,desc); desc=copy; }
	    else have_desc=0;
	}
	TIFFClose(in_tiff);

	if (median_size > 0 &&
	    apply_median(data,(int)width,(int)height,median_size,is_16bit) != 0) {
	    free(data); return 1;
	}
	if (sigma > 0.0 &&
	    apply_gaussian(data,(int)width,(int)height,sigma,is_16bit) != 0) {
	    free(data); return 1;
	}

	if ((out_tiff=TIFFOpen(output_file,"w")) == NULL) {
	    fprintf(stderr,"Error: Cannot create output TIFF file\n");
	    free(data); return 1;
	}
	TIFFSetField(out_tiff,TIFFTAG_IMAGEWIDTH,width);
	TIFFSetField(out_tiff,TIFFTAG_IMAGELENGTH,height);
	TIFFSetField(out_tiff,TIFFTAG_BITSPERSAMPLE,bits_per_sample);
	TIFFSetField(out_tiff,TIFFTAG_COMPRESSION,COMPRESSION_NONE);
	TIFFSetField(out_tiff,TIFFTAG_PHOTOMETRIC,PHOTOMETRIC_MINISBLACK);
	TIFFSetField(out_tiff,TIFFTAG_SAMPLESPERPIXEL,1);
	TIFFSetField(out_tiff,TIFFTAG_PLANARCONFIG,PLANARCONFIG_CONTIG);
	TIFFSetField(out_tiff,TIFFTAG_ROWSPERSTRIP,TIFFDefaultStripSize(out_tiff,0));
	if (have_desc) TIFFSetField(out_tiff,TIFFTAG_IMAGEDESCRIPTION,desc);
	TIFFSetField(out_tiff,TIFFTAG_ARTIST,"tif_mgf5_libtiff");

	for (row = 0; row < height; row++)
	    if (TIFFWriteScanline(out_tiff,(char *)data+(size_t)row*scanline_size,
				  row,0) < 0) {
		fprintf(stderr,"Error: Cannot write scanline %u\n",row);
		free(data); TIFFClose(out_tiff);
		return 1;
	    }
	TIFFClose(out_tiff);
	free(data);
	if (have_desc) free(desc);
	return 0;
}

/*----------------------------------------------------------------------*/

int	main(int argc,char *argv[])
{
	const char	*input_file,*output_file;
	int		median_kernel_size=1,result,i;
	double		gaussian_sigma=0.5;
	FILE		*flog;

	if (argc < 3 || argc > 5) {
	    fprintf(stderr,"Usage: %s <input_file> <output_file> [median_kernel_size] [gaussian_sigma]\n",argv[0]);
	    fprintf(stderr,"  median_kernel_size: 0-%d (default: 1, 0 to skip; 3 and 5 take the fast path)\n",MEDIAN_KMAX);
	    fprintf(stderr,"  gaussian_sigma: 0.0-100.0 (default: 0.5, 0.0 to skip)\n");
	    return 1;
	}
	input_file =argv[1];
	output_file=argv[2];

	if (argc > 3) {
	    median_kernel_size=atoi(argv[3]);
	    if (median_kernel_size < 0 || median_kernel_size > MEDIAN_KMAX) {
		fprintf(stderr,"Error: Invalid median kernel size (must be 0-%d)\n",MEDIAN_KMAX);
		return 1;
	    }
	}
	if (argc > 4) {
	    gaussian_sigma=atof(argv[4]);
	    if (gaussian_sigma < 0.0 || gaussian_sigma > 100.0) {
		fprintf(stderr,"Error: Invalid gaussian sigma (must be 0.0-100.0)\n");
		return 1;
	    }
	}

	result=process_tiff_file(input_file,output_file,median_kernel_size,
				 gaussian_sigma);

	/* append to log file (same format as tif_mgf.c) */
	if ((flog=fopen("cmd-hst.log","a")) != NULL) {
	    for (i = 0; i < argc; ++i) fprintf(flog,"%s ",argv[i]);
	    fprintf(flog,"\t");
	    fprintf(flog,"   %% median_kernel_size %d  gaussian_sigma %g\n",
		    median_kernel_size,gaussian_sigma);
	    fclose(flog);
	}
	return result;
}
