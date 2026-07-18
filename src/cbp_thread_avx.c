/*
 * cbp_thread_avx.c - convolution back projection, float + AVX2/FMA
 *
 * cbp_thread_int.c (整数固定小数点BP) の選択制の変種。量子化を行わない
 * ため精度は double 参照 (cbp_thread_nai.c) と実質一致する (実測で
 * int 版の約200倍良い)。速度は CPU の gather/permute 性能に依存し、
 * int 版の増分添字ループに劣る場合がある (2024年代の一般的なデスク
 * トップで実測 3〜5 割遅)。使う場合は makefile の CBP 変数を
 * cbp_thread_avx.c に切り替える。
 *
 * 構成 (cbp_thread_int.c の2フェーズ構成と同じ):
 *   フェーズA: 投影を K スレッドで分担して FFT フィルタ畳み込み
 *   フェーズB: 画像の行帯を K スレッドで分担して逆投影 (AVX2 8画素同時)
 * 逆投影は GPU 版 (cbp.cu BP_GMF) と同じ直接評価:
 *   r = X cosθ + Y sinθ - r0 を毎投影計算 (増分累積の誤差なし)。
 * 出力はスレッド数に依存しない (画素ごとの加算順序は常に m 昇順)。
 *
 * 制約: Float は float であること (-DFloat=double は cbp_thread_nai.c を使用)。
 * GCC では -mavx2 -mfma なしでも #pragma GCC target で有効化する。
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cbp.h"

#if	defined(__GNUC__) && !defined(__AVX2__)
#pragma	GCC target("avx2,fma")
#endif
#include <immintrin.h>

/* Float==float を強制 (double なら cbp_thread_nai.c を使う) */
typedef char	CBP_AVX_requires_float[(sizeof(Float)==sizeof(float))?1:-1];

#ifdef	WINDOWS
#include <windows.h>

#define THREAD_VAR	HANDLE
#define THREAD_FUNC	DWORD WINAPI
#define THREAD_ARG	LPVOID
#define THREAD_NULL	0

#define THREAD_CREATE(var,func,arg) \
	(var=CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)(func),(THREAD_ARG)(size_t)(arg),0,NULL))==0

#define THREAD_JOIN(var) \
	(void)WaitForSingleObject(var,INFINITE); (void)CloseHandle(var)
#else
#include <pthread.h>
#include <unistd.h>

#define THREAD_VAR	pthread_t
#define THREAD_FUNC	void *
#define THREAD_ARG	void *
#define THREAD_NULL	NULL

#define THREAD_CREATE(var,func,arg) \
	pthread_create(&(var),NULL,func,(THREAD_ARG)(size_t)(arg))

#define THREAD_JOIN(var) \
	if (pthread_join(var,NULL)) Error("can not join CBP_thread.")
#endif

/* 既定スレッド数: 実行PCの論理コア数-1 (0 になる場合は 1)。
   環境変数 CBP_THREADS または -DCBP_THREADS=n で上書き可能。 */
static int	DefaultThreads()
{
	int	n;
#ifdef	WINDOWS
	SYSTEM_INFO	si;

	GetSystemInfo(&si);
	n=(int)si.dwNumberOfProcessors-1;
#else
	n=(int)sysconf(_SC_NPROCESSORS_ONLN)-1;
#endif
	return (n<1)?1:n;
}

#ifndef	CBP_THREADS
#define CBP_THREADS	DefaultThreads()
#endif

extern void	Error();

double	Ramachandran(i)
int	i;
{
	return (i==0)?1.0/4.0:(i&1)?-1.0/(M_PI*M_PI*(double)i*(double)i):0.0;
}

double	Shepp(i)
int	i;
{
	return 2.0/(M_PI*M_PI*(1.0-4.0*(double)i*(double)i));
}

double	Chesler(i)
int	i;
{
	if (i==0)
	    return 1.0/8.0-1.0/(2.0*M_PI*M_PI);
	else if (i==1 || i==(-1))
	    return 1.0/16.0-1.0/(2.0*M_PI*M_PI);
	else if (i&1)
	    return -1.0/(2.0*M_PI*M_PI*(double)i*(double)i);
{
	double	i2=(double)i*(double)i,i21=i2-1.0;

	return -(i2+1.0)/(2.0*M_PI*M_PI*i21*i21);
}}

typedef struct {
		double	x,y;
	} Vector;

static void	FFT(L,E,Z,sign)
int		L,sign;
Vector		*E,*Z;
{
	int	j,k,m,n,p,q,L2=L>>1;
	Vector	z,*e,*a,*b;
	double	x,y,w;

	j=L-1;
	for (k=L-2; k>0; k--) {
	    for (m=L2; (j^=m)&m; m>>=1) ;

	    if (j>k) {
		z=Z[j]; Z[j]=Z[k]; Z[k]=z;
	    }
	}
	for (j=2, k=1, m=L2, n=L; n>1; n=m, m>>=1, k=j, j<<=1)
	    for (e=E, p=0; p<k; p++, e+=n) {
		x=e->x; y=e->y*(double)sign;
		for (a=Z+p, q=0; q<m; q++, a+=j) {
		    z=(*(b=a+k)); b->x=a->x-(w=z.x*x-z.y*y); a->x+=w;
				  b->y=a->y-(w=z.x*y+z.y*x); a->y+=w;
		}
	    }
}

static int		K,N,M,*H0,*dH;
static THREAD_VAR	*T;
static Float		**P,**Q,**F;
static Vector		*E,**Z;
static double		*G,*Sd,*Cd,*Ad;

#define POW2(N)		(1<<(int)ceil(log((double)(N))/log(2.0)))
#define ALLOC(type,noe)	(type *)malloc(sizeof(type)*(size_t)(noe))

Float	**InitCBP(Ni,Mi)
int	Ni,Mi;
{
	char	*env;
	Vector	*z;
	int	m,n,k,l,L=POW2(Ni),L2=L<<1,N2=Ni*Ni;
	double	a,da=M_PI/(double)L2;

#ifdef	PTW32_STATIC_LIB
	if (!pthread_win32_process_attach_np())
	    Error("pthreads-win32 initialization failed.");
#endif
	K=((env=getenv("CBP_THREADS"))==NULL)?CBP_THREADS:atoi(env);

	if (K<=0 || K>Mi) Error("bad number of CBP_THREADS.");

	if ((T	    =ALLOC(THREAD_VAR,K    ))==NULL ||
	    (P	    =ALLOC(Float *   ,Mi   ))==NULL ||
	    (P[0]   =ALLOC(Float     ,Mi*Ni))==NULL ||
	    (E	    =ALLOC(Vector    ,L2   ))==NULL ||
	    (Z	    =ALLOC(Vector *  ,K    ))==NULL ||
	    (Z[0]=z =ALLOC(Vector    ,K*L2 ))==NULL ||
	    (G	    =ALLOC(double    ,L2   ))==NULL ||
	    (Q	    =ALLOC(Float *   ,Mi   ))==NULL ||
	    (Q[0]   =ALLOC(Float     ,(size_t)Mi*(Ni+1)+8))==NULL ||
	    (Sd	    =ALLOC(double    ,Mi   ))==NULL ||
	    (Cd	    =ALLOC(double    ,Mi   ))==NULL ||
	    (Ad	    =ALLOC(double    ,Mi   ))==NULL ||
	    (H0	    =ALLOC(int	     ,Ni   ))==NULL ||
	    (dH	    =ALLOC(int	     ,Ni   ))==NULL ||
	    (F	    =ALLOC(Float *   ,Ni   ))==NULL ||
	    (F[0]   =ALLOC(Float     ,N2   ))==NULL) return NULL;

	/* Q[m] は補間用番兵 q[N]=0 を持つ (長さ N+1) */
	for (m=1; m<Mi; m++) {
	    P[m]=P[m-1]+Ni; Q[m]=Q[m-1]+(Ni+1);
	}
	for (n=1; n<Ni; n++) F[n]=F[n-1]+Ni;
	for (k=1; k<K ; k++) Z[k]=Z[k-1]+L2;

	E[0].x=E[L].y=1.0; z[0].x=Filter(L); z[L].x=Filter(0);
	E[0].y=E[L].x=	   z[0].y=	     z[L].y=0.0;
	for (l=L2-1, n=1; n<L; n++, l--) {
	    a=da*(double)n;
	    E[l].x=(-(E[n].x=cos(a))); z[l].x=z[n].x=Filter(L-n);
	    E[l].y=   E[n].y=sin(a)  ; z[l].y=z[n].y=0.0;
	}
	FFT(L2,E,z,-1);

	for (n=0; n<L2; n++) G[n]=z[n].x;

	N=Ni; M=Mi; return P;
}

static double	dr,r0,t0,dt,dtq;
static int	L,L2,N1,v1,v2;

/* フェーズA: 投影を K 本で分担してフィルタ畳み込み。
   スケール dtq を掛けて Q[m] に格納 (逆投影はそのまま加算するだけ)。 */
static THREAD_FUNC	ConvThread(ki)
THREAD_ARG		ki;
{
	size_t	k=(size_t)ki;
	int	m,l,n;
	Float	*p,*q;
	Vector	*z=Z[k];

	for (m=k; m<M; m+=K) {
	    for (	 l=0; l<L;	l++)   z[l].x=	    z[l].y=0.0;
	    for (p=P[m], n=0; n<N; n++, l++) { z[l].x=p[n]; z[l].y=0.0; }
	    for (	    ; n<L; n++, l++)   z[l].x=	    z[l].y=0.0;

	    FFT(L2,E,z,-1);

	    for (n=0; n<L2; n++) {
		z[n].x*=G[n]; z[n].y*=G[n];
	    }
	    FFT(L2,E,z, 1);

	    for (q=Q[m], n=0; n<N; n++) q[n]=(Float)(z[n].x*dtq);
	    q[N]=(Float)0.0;
	}
	return THREAD_NULL;
}

/* 部分ブロック用ストアマスク: mtab+(8-rem) から 8 lane 分ロード */
static const int	mtab[16]={-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0};

/* フェーズB: 画像の行帯 [blo,bhi] を分担し、8画素同時に逆投影。
   r = X cosθ + Y sinθ - r0 → n=floor(r) → q[n]+(q[n+1]-q[n])(r-n) を
   マスク付き gather + FMA で累積する。範囲外 lane は 0 寄与。 */
static THREAD_FUNC	BPThread(ki)
THREAD_ARG		ki;
{
	size_t	k=(size_t)ki;
	int	span=v2-v1+1,
		blo=v1+(int)(((long long)span*(long long)k)/(long long)K),
		bhi=v1+(int)(((long long)span*(long long)(k+1))/(long long)K)-1,
		v,h,m,h1,h2,rem,n0,fast;
	__m256	lane,acc,rh,fl,frac,g0,g1,g1a,okf,va,vb;
	__m256i	Nv,onei,zeroi,seveni,eighti,n,ok,smask,idx,idx1,m8;
	Float	*fr;
	const Float	*q;
	double	vd,bd,r7,rlo;

	if (blo>bhi) return THREAD_NULL;

	lane  =_mm256_set_ps(7.0f,6.0f,5.0f,4.0f,3.0f,2.0f,1.0f,0.0f);
	Nv    =_mm256_set1_epi32(N);
	onei  =_mm256_set1_epi32(1);
	seveni=_mm256_set1_epi32(7);
	eighti=_mm256_set1_epi32(8);
	zeroi =_mm256_setzero_si256();

	for (v=blo; v<=bhi; v++) {
	    fr=F[v]; vd=(double)v;
	    h1=H0[v]; h2=h1+dH[v];

	    for (h=h1; h<=h2; h+=8) {
		rem=h2-h+1;
		smask=_mm256_loadu_si256((const __m256i *)
					 (mtab+((rem>=8)?0:8-rem)));
		acc=_mm256_setzero_ps();

		for (m=0; m<M; m++) {
		    q=Q[m];
		    bd=Ad[m]-vd*Sd[m]+(double)h*Cd[m];
		    r7=bd+7.0*Cd[m];
		    rlo=(bd<r7)?bd:r7;
		    /* 切り捨てはキャストで足りる (真の floor と食い違う負側や
		       丸め誤差は後段のベクトル検証で gather 経路へ落ちる) */
		    n0=(int)rlo;

		    rh  =_mm256_fmadd_ps(lane,_mm256_set1_ps((float)Cd[m]),
					      _mm256_set1_ps((float)bd));
		    fl  =_mm256_floor_ps(rh);
		    n   =_mm256_cvtps_epi32(fl);
		    frac=_mm256_sub_ps(rh,fl);

		    fast=0;
		    if (rem==8 && n0>=0 && n0<=N-8) {
			/* float 丸めで lane の添字が [n0,n0+7] を外れて
			   いないことをベクトル側でも検証してから使う */
			idx=_mm256_sub_epi32(n,_mm256_set1_epi32(n0));
			m8 =_mm256_or_si256(_mm256_cmpgt_epi32(zeroi,idx),
					    _mm256_cmpgt_epi32(idx,seveni));
			fast=_mm256_testz_si256(m8,m8);
		    }
		    if (fast) {
			/* 内部ブロック高速路: 8 lane の添字は [n0,n0+7] の
			   連続9要素に収まるので gather の代わりに
			   連続ロード + permute で引く */
			idx1=_mm256_add_epi32(idx,onei);
			va  =_mm256_loadu_ps(q+n0);
			vb  =_mm256_loadu_ps(q+n0+8);	/* 先頭要素のみ使用 */

			g0 =_mm256_permutevar8x32_ps(va,idx);
			g1a=_mm256_permutevar8x32_ps(va,idx1);
			m8 =_mm256_cmpeq_epi32(idx1,eighti);
			g1 =_mm256_blendv_ps(g1a,
				_mm256_permutevar8x32_ps(vb,zeroi),
				_mm256_castsi256_ps(m8));
		    }
		    else {
			/* 端ブロック/範囲外 lane 混在: マスク付き gather */
			ok  =_mm256_andnot_si256(_mm256_cmpgt_epi32(zeroi,n),
						 _mm256_cmpgt_epi32(Nv,n));
			ok  =_mm256_and_si256(ok,smask);
			okf =_mm256_castsi256_ps(ok);

			g0=_mm256_mask_i32gather_ps(_mm256_setzero_ps(),q,n,
						    okf,4);
			g1=_mm256_mask_i32gather_ps(_mm256_setzero_ps(),q,
						    _mm256_add_epi32(n,onei),
						    okf,4);
		    }
		    acc=_mm256_add_ps(acc,
			_mm256_fmadd_ps(_mm256_sub_ps(g1,g0),frac,g0));
		}
		_mm256_maskstore_ps(fr+h,smask,acc);
	    }
	}
	return THREAD_NULL;
}

Float	**CBP(dri,r0i,t0i)
double	dri,r0i,t0i;
{
	double	N12,R,R2,y,x,t;
	int	v,k,m;

	dr=dri; r0=r0i; t0=t0i;
	L2=(L=POW2(N))<<1;
	dt=M_PI/(double)M;
	dtq=dt/(dr*(double)L2);
	N1=N-1; N12=(double)N1/2.0; R=N12-fabs(N12+r0); R2=R*R;

#ifndef	__INTEL_COMPILER
	v1=(int)ceil(N12-R); v2=(int)(N12+R);
	for (y=N12-(double)(v=v1); v<=v2; v++, y-=1.0) {
	    x=sqrt(R2-y*y); dH[v]=(int)(N12+x)-(H0[v]=(int)ceil(N12-x));
	}
#else
	v1=N; v2=(-1);
	for (y=N12, v=0; v<N; v++, y-=1.0)
	    if (fabs(y)<=R) {
		x=sqrt(R2-y*y); dH[v]=(int)(N12+x)-(H0[v]=(int)ceil(N12-x));

		if (v<v1) v1=v; v2=v;
	    }
#endif
	/* フェーズA: フィルタ畳み込み (投影分割) */
	for (k=1; k<K; k++)
	    if (THREAD_CREATE(T[k],ConvThread,k))
		Error("can not create CBP_thread.");

	(void)ConvThread((THREAD_ARG)0);

	for (k=1; k<K; k++) {
	    THREAD_JOIN(T[k]);
	}
	/* 投影角表 (逆投影の直接評価用; 全て double で一度だけ) */
	for (m=0; m<M; m++) {
	    t=t0+dt*(double)m;
	    Sd[m]=sin(t); Cd[m]=cos(t);
	    Ad[m]=N12*(Sd[m]-Cd[m])-r0;
	}
	(void)memset(F[0],0,sizeof(Float)*(size_t)N*(size_t)N);

	/* フェーズB: 逆投影 (行帯分割, AVX2) */
	for (k=1; k<K; k++)
	    if (THREAD_CREATE(T[k],BPThread,k))
		Error("can not create CBP_thread.");

	(void)BPThread((THREAD_ARG)0);

	for (k=1; k<K; k++) {
	    THREAD_JOIN(T[k]);
	}
	return F;
}

void	TermCBP()
{
	free(F[0]); free(F); free(dH); free(H0);
	free(Ad); free(Cd); free(Sd);
	free(Q[0]); free(Q); free(G); free(Z[0]); free(Z); free(E);
	free(P[0]); free(P);
	free(T);

#ifdef	PTW32_STATIC_LIB
	(void)pthread_win32_process_detach_np();
#endif
}
