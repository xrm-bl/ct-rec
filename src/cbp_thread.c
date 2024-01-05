
#include <stdlib.h>
#include <math.h>
#include "cbp.h"

#if	defined(WINDOWS) || defined(_WIN32)
#include <windows.h>

#define THREAD_VAR	HANDLE
#define THREAD_FUNC	DWORD WINAPI
#define THREAD_ARG	LPVOID
#define THREAD_NULL	0

#define THREAD_CREATE(var,func,arg)	\
	(var=CreateThread(NULL,0,func,(THREAD_ARG)(arg),0,NULL))==0

#define THREAD_JOIN(var)	\
	WaitForSingleObject(var,INFINITE)==WAIT_FAILED || CloseHandle(var)==0
#else
#include <pthread.h>

#define THREAD_VAR	pthread_t
#define THREAD_FUNC	void *
#define THREAD_ARG	void *
#define THREAD_NULL	NULL

#define THREAD_CREATE(var,func,arg)	\
	pthread_create(&(var),NULL,func,(THREAD_ARG)(arg))

#define THREAD_JOIN(var)	pthread_join(var,NULL)
#endif

#ifndef	CBP_THREADS
#define CBP_THREADS	8
#endif

#if	defined(__x86_64__) || defined(__ia64__)
#define __64BIT_DATABUS
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
static Float		**P,**Q,***F;
static Vector		*E,*Z,*G;

#ifdef	__64BIT_DATABUS
typedef struct {
		Float	x,y;
	} Work;

static Work	**W;
#endif

#define POW2(N)		(1<<(int)ceil(log((double)(N))/log(2.0)))
#define OVER(R,M)	(int)ceil(M_PI*(R)/(double)(M))
#define ALLOC(type,noe)	(type *)malloc(sizeof(type)*(size_t)(noe))

Float	**InitCBP(Ni,Mi)
int	Ni,Mi;
{
	char	*env;
	int	m,n,k,
		L=POW2(Ni),L2=L<<1,MO=Mi*OVER((double)(Ni-1)/2.0,Mi),N2=Ni*Ni;
	double	a,da=M_PI/(double)L2;

#ifdef	PTW32_STATIC_LIB
	if (!pthread_win32_process_attach_np())
	    Error("pthreads-win32 initialization failed.");
#endif
	K=((env=getenv("CBP_THREADS"))==NULL)?CBP_THREADS:atoi(env);

	if (K<=0 || K>Mi) Error("bad number of CBP_THREADS.");

	if ((T  =ALLOC(THREAD_VAR,K        ))==NULL ||
	    (P  =ALLOC(Float *   ,Mi       ))==NULL ||
	    (*P =ALLOC(Float     ,Mi*Ni    ))==NULL ||
	    (E  =ALLOC(Vector    ,L2       ))==NULL ||
	    (Z  =ALLOC(Vector    ,L2       ))==NULL ||
	    (G  =ALLOC(Vector    ,L2       ))==NULL ||
	    (Q  =ALLOC(Float *   ,MO       ))==NULL ||
#ifdef	__64BIT_DATABUS
	    (*Q =ALLOC(Float     ,MO*Ni    ))==NULL ||
	    (W  =ALLOC(Work *    ,K        ))==NULL ||
	    (*W =ALLOC(Work      ,K*Ni     ))==NULL ||
#else
	    (*Q =ALLOC(Float     ,MO*(Ni+1)))==NULL ||
#endif
	    (H0 =ALLOC(int       ,Ni       ))==NULL ||
	    (dH =ALLOC(int       ,Ni       ))==NULL ||
	    (F  =ALLOC(Float **  ,K        ))==NULL ||
	    (*F =ALLOC(Float *   ,K*Ni     ))==NULL ||
	    (**F=ALLOC(Float     ,K*N2     ))==NULL) return NULL;

	for (m=1; m<Mi; m++) P[m]=P[m-1]+Ni;
#ifdef	__64BIT_DATABUS
	for (m=1; m<MO; m++) Q[m]=Q[m-1]+Ni;
	for (k=1; k<K ; k++) W[k]=W[k-1]+Ni;
#else
	Q[0][Ni]=0.0;
	for (m=1; m<MO; m++) {
	    Q[m]=Q[m-1]+Ni+1; Q[m][Ni]=0.0;
	}
#endif
	for (n=1; n<Ni; n++) F[0][n]=F[0][n-1]+Ni;
	for (k=1; k<K ; k++) {
	    F[k]=F[k-1]+Ni; for (n=0; n<Ni; n++) F[k][n]=F[k-1][n]+N2;
	}
	for (n=0; n<L2; n++) {
	    a=da*(double)n; E[n].x=cos(a); E[n].y=sin(a);

	    Z[n].x=Filter(n-L); Z[n].y=0.0;
	}
	FFT(L2,E,Z,-1);

	for (n=0; n<L2; n++) G[n]=Z[n];

	N=Ni; M=Mi; return P;
}

static double	r0,t0,dt,N12;
static int	MO,v1,v2;

static THREAD_FUNC	BP_thread(ki)
THREAD_ARG		ki;
{
	size_t	k=(size_t)ki;
	int	v,h,m,n;
	Float	*q,**fv,*fh,**f=F[k];
	double	t,y,x,rv,rh;
#ifdef	__64BIT_DATABUS
	Work	*w=W[k];
	Float	q2,q1;
#endif
	for (v=0; v<N; v++)
	for (h=0; h<N; h++) f[v][h]=0.0;

	for (m=k; m<MO; m+=K) {
	    q=Q[m];
#ifdef	__64BIT_DATABUS
	    w+=N; q+=N;
	    for (q2=0.0, n=0; n<N; n++, q2=q1) {
		--w; --q; w->y=q2-(q1=w->x=(*q));
	    }
#endif
	    t=t0+dt*(double)m; y=sin(t); x=cos(t);
	    rv=N12*(y-x)-(double)v1*y-r0; fv=f+v1;
	    for (v=v1; v<=v2; v++, rv-=y, fv++) {
		h=H0[v]; rh=(double)h*x+rv; fh=(*fv)+h;
		for (h=dH[v]; h>=0; h--, rh+=x, fh++) {
		    n=(int)rh;
#ifdef	__64BIT_DATABUS
		    *fh+=(w[n].x+w[n].y*(rh-(double)n));
#else
		    *fh+=(q[n]+(q[n+1]-q[n])*(rh-(double)n));
#endif
		}
	    }
	}
	return THREAD_NULL;
}

Float	**CBP(dr,r0i,t0i)
double	dr,r0i,t0i;
{
	double	dtq,R,R2,x,y,r,d;
	int	O,m,l,n,o,v,k,h,L=POW2(N),L2=L<<1;
	Float	*p,*q,*q0,*qM,**Fk,**F0=F[0];
#ifdef	__64BIT_DATABUS
	Float	q2,q1;
#endif
	r0=r0i; t0=t0i;
	N12=(double)(N-1)/2.0; R=N12-fabs(N12+r0); R2=R*R;
	O=OVER(R,M); dtq=(dt=M_PI/(double)(MO=M*O))/(dr*(double)L2);

	for (m=0; m<M; m++) {
	    for (	 l=0; l<L;	l++)   Z[l].x=	    Z[l].y=0.0;
	    for (p=P[m], n=0; n<N; n++, l++) { Z[l].x=p[n]; Z[l].y=0.0; }
	    for (	    ; n<L; n++, l++)   Z[l].x=	    Z[l].y=0.0;

	    FFT(L2,E,Z,-1);

	    for (n=0; n<L2; n++) {
		x=Z[n].x;
		y=Z[n].y; Z[n].x=x*G[n].x-y*G[n].y;
			  Z[n].y=x*G[n].y+y*G[n].x;
	    }
	    FFT(L2,E,Z, 1);

	    for (q=Q[m*O], n=0; n<N; n++) q[n]=Z[n].x*dtq;
	}
	if (O>1) {
#ifdef	__64BIT_DATABUS
	    Z+=N; q0=Q[0]+N;
	    for (q2=0.0, n=0; n<N; n++, q2=q1) {
		--Z; --q0; Z->y=q2-(q1=Z->x=(*q0));
	    }
	    qM=Q[1]; r=(-2.0*r0); l=(int)r; d=r-(double)l;
	    for (n=0; n<N; n++, l--) qM[n]=(l>=0 && l<N)?Z[l].x+Z[l].y*d:0.0;
#else
	    q0=Q[0]; qM=Q[1]; r=(-2.0*r0); l=(int)r; d=r-(double)l;
	    for (n=0; n<N; n++, l--)
		qM[n]=(l>=0 && l<N)?q0[l]+(q0[l+1]-q0[l])*d:0.0;
#endif
	    for (m=MO; (m-=O)>=0; qM=q0) {
		for (q0=Q[m], n=0; n<N; n++)
		    Z[n].y=(qM[n]-(Z[n].x=q0[n]))/(double)O;

		for (o=1; o<O; o++)
		    for (q=Q[m+o], n=0; n<N; n++)
			q[n]=Z[n].x+Z[n].y*(double)o;
	    }
	}

#ifdef	__INTEL_COMPILER
	v1=N; v2=(-1);
	for (y=N12, v=0; v<N; v++, y-=1.0)
	    if (fabs(y)<=R) {
		x=sqrt(R2-y*y); dH[v]=(int)(N12+x)-(H0[v]=(int)ceil(N12-x));

		if (v<v1) v1=v; v2=v;
	    }
#else
	v1=(int)ceil(N12-R); v2=(int)(N12+R);
	for (y=N12-(double)(v=v1); v<=v2; v++, y-=1.0) {
	    x=sqrt(R2-y*y); dH[v]=(int)(N12+x)-(H0[v]=(int)ceil(N12-x));
	}
#endif
	for (k=1; k<K; k++)
	    if (THREAD_CREATE(T[k],BP_thread,k))
		Error("can not create BP_thread.");

	(void)BP_thread((THREAD_ARG)0);

	for (k=1; k<K; k++) {
	    if (THREAD_JOIN(T[k])) Error("can not join BP_thread.");

	    Fk=F[k]; for (v=0; v<N; v++)
		     for (h=0; h<N; h++) F0[v][h]+=Fk[v][h];
	}
	return F0;
}

void	TermCBP()
{
	free(**F); free(*F); free(F); free(dH); free(H0);

#ifdef	__64BIT_DATABUS
	free(*W); free(W);
#endif
	free(*Q); free(Q); free(G); free(Z); free(E);
	free(*P); free(P);
	free(T);

#ifdef	PTW32_STATIC_LIB
	(void)pthread_win32_process_detach_np();
#endif
}
