
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cbp.h"

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

#ifndef	Int
#define Int	int
#endif

#ifndef	P_BITS
#define P_BITS	12
#endif

#ifndef	R_BITS
#define R_BITS	9
#endif

#ifndef	Q_BITS
#define Q_BITS	(31-P_BITS-R_BITS)
#endif

#ifndef	S_BITS
#define S_BITS	(30-P_BITS-R_BITS)
#endif

#if	defined(__x86_64__) || defined(__ia64__)
#define	__64BIT_DATABUS
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

#ifdef	__64BIT_DATABUS
typedef struct {
		Int	x,y;
	} IPair;
#endif

static int		K,N,M,*H0,*dH;
static THREAD_VAR	*T;
static Float		**P,**Q,**F;
static Vector		*E,**Z;
static double		*G,*Qmax;
#ifdef	__64BIT_DATABUS
static IPair		**QI;
#else
static Int		**QI;
#endif
static Int		*SI,*CI;

#define POW2(N)		(1<<(int)ceil(log((double)(N))/log(2.0)))
#define ALLOC(type,noe)	(type *)malloc(sizeof(type)*(size_t)(noe))

/* ---- Truncation (cupping) pad, controlled by env var PAD_THRESH ----
   When PAD_THRESH (a plain ratio) > 0 and the mean amplitude of the
   sinogram's outermost columns exceeds PAD_THRESH times the overall
   mean (i.e. the sample overfills the field of view), each projection
   is extended on both sides by N/2 samples holding the edge value
   under a cosine decay before filtering.  The step left by the usual
   zero padding, convolved with the ramp kernel's long -1/(2 pi^2 r^2)
   tail, is what produces the bright rim (cupping); the smooth decay
   removes it.  Backprojection stays at N (cost unchanged; only the
   FFT length doubles).  Default (unset/<=0) is fully OFF and
   bit-identical to the previous behaviour. */
static double	PadThr=-1.0;
static int	PadW=0;			/* pad width for this CBP() call (0=off) */

static double	PadThreshold()
{
	char	*e;

	if (PadThr<0.0)
	    PadThr=((e=getenv("PAD_THRESH"))!=NULL && atof(e)>0.0)?atof(e):0.0;
	return PadThr;
}

/* with the pad enabled the FFT must hold N+2*(N/2) input samples */
#define PADPOW2(n)	((PadThreshold()>0.0)?POW2((n)+(n)/2):POW2(n))

static void	DetectPad(Float **Pi,int Ni,int Mi)
{
	static int	said=0;
	double		tot=0.0,edge=0.0;
	int		m,n;

	PadW=0;
	if (PadThreshold()<=0.0) return;
	for (m=0; m<Mi; m++) {
	    for (n=0; n<Ni; n++) tot+=fabs((double)Pi[m][n]);
	    edge+=fabs((double)Pi[m][0])+fabs((double)Pi[m][Ni-1]);
	}
	if (tot>0.0 &&
	    edge/(2.0*(double)Mi)>PadThr*tot/((double)Mi*(double)Ni)) {
	    PadW=Ni/2;
	    if (!said) {
		fprintf(stderr,"PAD_THRESH: truncation pad enabled (W=%d)\n",PadW);
		said=1;
	    }
	}
}

static void	FillPad(Vector *z,Float *pr,int Li,int Ni)
{
	double	w;
	int	n;

	for (n=1; n<=PadW; n++) {
	    w=0.5*(1.0+cos(M_PI*(double)n/(double)PadW));
	    z[Li-n].x	  =pr[0]   *w;
	    z[Li+Ni-1+n].x=pr[Ni-1]*w;
	}
}

Float	**InitCBP(Ni,Mi)
int	Ni,Mi;
{
	char	*env;
	Vector	*z;
	int	m,n,k,l,L=PADPOW2(Ni),L2=L<<1,N2=Ni*Ni;
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
	    (Q[0]   =ALLOC(Float     ,Mi*Ni))==NULL ||
	    (Qmax   =ALLOC(double    ,K    ))==NULL ||
#ifdef	__64BIT_DATABUS
	    (QI	    =ALLOC(IPair *   ,Mi   ))==NULL ||
	    (QI[0]  =ALLOC(IPair     ,(size_t)Mi*(Ni+2)))==NULL ||
#else
	    (QI	    =ALLOC(Int *     ,Mi   ))==NULL ||
	    (QI[0]  =ALLOC(Int	     ,(size_t)Mi*(Ni+2)))==NULL ||
#endif
	    (SI	    =ALLOC(Int	     ,Mi   ))==NULL ||
	    (CI	    =ALLOC(Int	     ,Mi   ))==NULL ||
	    (H0	    =ALLOC(int	     ,Ni   ))==NULL ||
	    (dH	    =ALLOC(int	     ,Ni   ))==NULL ||
	    (F	    =ALLOC(Float *   ,Ni   ))==NULL ||
	    (F[0]   =ALLOC(Float     ,N2   ))==NULL) return NULL;

	/* QI[m] は前後に番兵 (-1 と N) を持つ量子化投影 */
	QI[0]+=1;
	for (m=1; m<Mi; m++) {
	    P[m]=P[m-1]+Ni; Q[m]=Q[m-1]+Ni; QI[m]=QI[m-1]+(Ni+2);
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

static double	dr,r0,t0,dt,qf;
static int	L,L2,N1,v1,v2;
static Int	R0g;

/* フェーズA: 投影を K 本で分担してフィルタ畳み込み (FFT)。
   Q[m] を作り、スレッドごとの最大絶対値を Qmax[k] へ。 */
static THREAD_FUNC	ConvThread(ki)
THREAD_ARG		ki;
{
	size_t	k=(size_t)ki;
	int	m,l,n;
	Float	*p,*q;
	Vector	*z=Z[k];
	double	q0,Q0=0.0;

	for (m=k; m<M; m+=K) {
	    for (	 l=0; l<L;	l++)   z[l].x=	    z[l].y=0.0;
	    for (p=P[m], n=0; n<N; n++, l++) { z[l].x=p[n]; z[l].y=0.0; }
	    for (	    ; n<L; n++, l++)   z[l].x=	    z[l].y=0.0;
	    if (PadW>0) FillPad(z,p,L,N);	/* truncation pad (PAD_THRESH) */

	    FFT(L2,E,z,-1);

	    for (n=0; n<L2; n++) {
		z[n].x*=G[n]; z[n].y*=G[n];
	    }
	    FFT(L2,E,z, 1);

	    for (q=Q[m], n=0; n<N; n++)
		if ((q0=fabs(q[n]=z[n].x))>Q0) Q0=q0;
	}
	Qmax[k]=Q0;

	return THREAD_NULL;
}

/* フェーズB: 画像の行帯 [blo,bhi] を K 本で分担して逆投影。
   全投影を舐め、自分の帯にのみ書くので共有画像1枚で足りる
   (旧実装のスレッド毎の画像複製と直列合算を廃止)。
   量子化スケールはグローバル Q0 なのでスレッド数に依らず結果は同一。 */
static THREAD_FUNC	BPThread(ki)
THREAD_ARG		ki;
{
	size_t	k=(size_t)ki;
	int	span=v2-v1+1,
		blo=v1+(int)(((long long)span*(long long)k)/(long long)K),
		bhi=v1+(int)(((long long)span*(long long)(k+1))/(long long)K)-1,
		v,h,m,n;
	Int	s,c,rv,**fv,rh,*fh,
		**fi=(Int **)F;
#ifndef	__64BIT_DATABUS
	Int	*qi;
	int	n1;
#else
	IPair	*qi;
#endif
	if (blo>bhi) return THREAD_NULL;

	for (m=0; m<M; m++) {
	    qi=QI[m]; s=SI[m]; c=CI[m];

	    rv=(N1*(s-c)+1)/2-blo*s-R0g; fv=fi+blo;
	    for (v=blo; v<=bhi; v++, rv-=s, fv++) {
		h=H0[v]; rh=h*c+rv; fh=(*fv)+h;
		for (h=dH[v]; h>=0; h--, rh+=c, fh++) {
#ifndef	__64BIT_DATABUS
	n=(n1=(int)(rh>>R_BITS+S_BITS))-1;
	*fh+=(qi[n]*(1<<R_BITS)+(qi[n1]-qi[n])*((rh>>S_BITS)&((1<<R_BITS)-1)));
#else
	n=(int)(rh>>R_BITS+S_BITS)-1;
	*fh+=(qi[n].x*(1<<R_BITS)+qi[n].y*((rh>>S_BITS)&((1<<R_BITS)-1)));
#endif
		}
	    }
	}
	/* 自分の帯を int 累積値から Float へ変換 (降順で in-place 安全) */
	for (v=bhi; v>=blo; v--)
	for (h=N1; h>=0; h--) F[v][h]=qf*(double)fi[v][h];

	return THREAD_NULL;
}

Float	**CBP(dri,r0i,t0i)
double	dri,r0i,t0i;
{
	double	N12,R,R2,y,x,q0,Q0,t;
	int	v,k,h,m,n;

	dr=dri; r0=r0i; t0=t0i;
	L2=(L=PADPOW2(N))<<1;
	DetectPad(P,N,M);
	dt=M_PI/(double)M;
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

	for (Q0=Qmax[0], k=1; k<K; k++) {
	    THREAD_JOIN(T[k]);

	    if (Qmax[k]>Q0) Q0=Qmax[k];
	}
	/* グローバル Q0 で量子化 (スレッド数に依らず同一の丸め) */
	q0=ldexp(1.0/Q0,Q_BITS);
	for (m=0; m<M; m++) {
	    Float	*q=Q[m];
#ifndef	__64BIT_DATABUS
	    Int		*qi=QI[m];

	    qi[-1]=qi[N]=0;
	    for (n=0; n<N; n++) qi[n]=(Int)floor(q[n]*q0+0.5);
#else
	    IPair	*qi=QI[m];

	    qi[-1].x=qi[N].x=qi[N].y=0;
	    for (n=N1; n>=0; n--)
		qi[n].y=qi[n+1].x-(qi[n].x=(Int)floor(q[n]*q0+0.5));

	    qi[-1].y=qi[0].x;
#endif
	    t=t0+dt*(double)m;
	    SI[m]=(Int)floor(ldexp(sin(t),R_BITS+S_BITS)+0.5);
	    CI[m]=(Int)floor(ldexp(cos(t),R_BITS+S_BITS)+0.5);
	}
	R0g=(Int)floor(ldexp(r0-1.0,R_BITS+S_BITS)+0.5);
	qf=ldexp(dt*Q0/(dr*(double)L2),-(Q_BITS+R_BITS));

	(void)memset(F[0],0,sizeof(Float)*(size_t)N*(size_t)N);

	/* フェーズB: 逆投影 (行帯分割) */
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
	free(CI); free(SI); free(QI[0]-1); free(QI); free(Qmax);
	free(Q[0]); free(Q); free(G); free(Z[0]); free(Z); free(E);
	free(P[0]); free(P);
	free(T);

#ifdef	PTW32_STATIC_LIB
	(void)pthread_win32_process_detach_np();
#endif
}
