
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
static Float		**P,***F;
static Vector		*E,**Z,*G;

#define POW2(N)		(1<<(int)ceil(log((double)(N))/log(2.0)))
#define ALLOC(type,noe)	(type *)malloc(sizeof(type)*(size_t)(noe))

Float	**InitCBP(Ni,Mi)
int	Ni,Mi;
{
	char	*env;
	Vector	*z;
	int	m,n,k,L=POW2(Ni),L2=L<<1,N2=Ni*Ni;
	double	a,da=M_PI/(double)L2;

#ifdef	PTW32_STATIC_LIB
	if (!pthread_win32_process_attach_np())
	    Error("pthreads-win32 initialization failed.");
#endif
	K=((env=getenv("CBP_THREADS"))==NULL)?CBP_THREADS:atoi(env);

	if (K<=0 || K>Mi) Error("bad number of CBP_THREADS.");

	if ((T   =ALLOC(THREAD_VAR,K    ))==NULL ||
	    (P   =ALLOC(Float *   ,Mi   ))==NULL ||
	    (*P  =ALLOC(Float     ,Mi*Ni))==NULL ||
	    (E   =ALLOC(Vector    ,L2   ))==NULL ||
	    (Z   =ALLOC(Vector *  ,K    ))==NULL ||
	    (*Z=z=ALLOC(Vector    ,K*L2 ))==NULL ||
	    (G   =ALLOC(Vector    ,L2   ))==NULL ||
	    (H0  =ALLOC(int       ,Ni   ))==NULL ||
	    (dH  =ALLOC(int       ,Ni   ))==NULL ||
	    (F   =ALLOC(Float **  ,K    ))==NULL ||
	    (*F  =ALLOC(Float *   ,K*Ni ))==NULL ||
	    (**F =ALLOC(Float     ,K*N2 ))==NULL) return NULL;

	for (m=1; m<Mi; m++) P[m]=P[m-1]+Ni;
	for (n=1; n<Ni; n++) F[0][n]=F[0][n-1]+Ni;
	for (k=1; k<K ; k++) {
	    Z[k]=Z[k-1]+L2;
	    F[k]=F[k-1]+Ni; for (n=0; n<Ni; n++) F[k][n]=F[k-1][n]+N2;
	}
	for (n=0; n<L2; n++) {
	    a=da*(double)n; E[n].x=cos(a); E[n].y=sin(a);

	    z[n].x=Filter(n-L); z[n].y=0.0;
	}
	FFT(L2,E,z,-1);

	for (n=0; n<L2; n++) G[n]=z[n];

	N=Ni; M=Mi; return P;
}

static double	r0,t0,dtq,dt,N12;
static int	L,L2,v1,v2;

static THREAD_FUNC	CBP_thread(ki)
THREAD_ARG		ki;
{
	size_t	k=(size_t)ki;
	int	v,h,m,l,n;
	Float	*p,**fv,*fh,**f=F[k];
	Vector	*z=Z[k];
	double	x,y,t,rv,rh;
#ifndef	__64BIT_DATABUS
	Float	*q=(Float *)z;
#endif
	for (v=0; v<N; v++)
	for (h=0; h<N; h++) f[v][h]=0.0;

	for (m=k; m<M; m+=K) {
	    for (	 l=0; l<L;	l++)   z[l].x=	    z[l].y=0.0;
	    for (p=P[m], n=0; n<N; n++, l++) { z[l].x=p[n]; z[l].y=0.0; }
	    for (	    ; n<L; n++, l++)   z[l].x=	    z[l].y=0.0;

	    FFT(L2,E,z,-1);

	    for (n=0; n<L2; n++) {
		x=z[n].x;
		y=z[n].y; z[n].x=x*G[n].x-y*G[n].y;
			  z[n].y=x*G[n].y+y*G[n].x;
	    }
	    FFT(L2,E,z, 1);

#ifndef	__64BIT_DATABUS
	    for (n=0; n<N; n++) q[n]=z[n].x*dtq;
	    q[N]=0.0;
#else
	    z+=N; z->x=0.0;
	    for (n=0; n<N; n++) {
		--z; z->y=z[1].x-(z->x*=dtq);
	    }
#endif
	    t=t0+dt*(double)m; y=sin(t); x=cos(t);
	    rv=N12*(y-x)-(double)v1*y-r0; fv=f+v1;
	    for (v=v1; v<=v2; v++, rv-=y, fv++) {
		h=H0[v]; rh=(double)h*x+rv; fh=(*fv)+h;
		for (h=dH[v]; h>=0; h--, rh+=x, fh++) {
		    n=(int)rh;
#ifndef	__64BIT_DATABUS
		    *fh+=(q[n]+(q[n+1]-q[n])*(rh-(double)n));
#else
		    *fh+=(z[n].x+z[n].y*(rh-(double)n));
#endif
		}
	    }
	}
	return THREAD_NULL;
}

Float	**CBP(dr,r0i,t0i)
double	dr,r0i,t0i;
{
	double	R,R2,y,x;
	int	v,k,h;
	Float	**Fk,**F0=F[0];

	r0=r0i; t0=t0i;
	L2=(L=POW2(N))<<1;
	dtq=(dt=M_PI/(double)M)/(dr*(double)L2);
	N12=(double)(N-1)/2.0; R=N12-fabs(N12+r0); R2=R*R;

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
	for (k=1; k<K; k++)
	    if (THREAD_CREATE(T[k],CBP_thread,k))
		Error("can not create CBP_thread.");

	(void)CBP_thread((THREAD_ARG)0);

	for (k=1; k<K; k++) {
	    if (THREAD_JOIN(T[k])) Error("can not join CBP_thread.");

	    Fk=F[k]; for (v=0; v<N; v++)
		     for (h=0; h<N; h++) F0[v][h]+=Fk[v][h];
	}
	return F0;
}

void	TermCBP()
{
	free(**F); free(*F); free(F); free(dH); free(H0);
	free(G); free(*Z); free(Z); free(E);
	free(*P); free(P);
	free(T);

#ifdef	PTW32_STATIC_LIB
	(void)pthread_win32_process_detach_np();
#endif
}
