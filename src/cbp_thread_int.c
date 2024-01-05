
#include <stdlib.h>
#include <math.h>
#include "cbp.h"

#ifdef	WINDOWS
#include <windows.h>

#define THREAD_VAR	HANDLE
#define THREAD_FUNC	DWORD WINAPI
#define THREAD_ARG	LPVOID
#define THREAD_NULL	0

#define THREAD_CREATE(var,func,arg) \
	(var=CreateThread(NULL,0,func,(THREAD_ARG)(arg),0,NULL))==0

#define THREAD_JOIN(var) \
	(void)WaitForSingleObject(var,INFINITE); (void)CloseHandle(var)
#else
#include <pthread.h>

#define THREAD_VAR	pthread_t
#define THREAD_FUNC	void *
#define THREAD_ARG	void *
#define THREAD_NULL	NULL

#define THREAD_CREATE(var,func,arg) \
	pthread_create(&(var),NULL,func,(THREAD_ARG)(arg))

#define THREAD_JOIN(var) \
	if (pthread_join(var,NULL)) Error("can not join CBP_thread.")
#endif

#ifndef	CBP_THREADS
#define CBP_THREADS	8
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

static int		K,N,M,*H0,*dH;
static THREAD_VAR	*T;
static Float		**P,**Q,***F;
static Vector		*E,**Z;
static double		*G;

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
	    (Q[0]   =ALLOC(Float     ,Mi*Ni))==NULL ||
	    (H0	    =ALLOC(int	     ,Ni   ))==NULL ||
	    (dH	    =ALLOC(int	     ,Ni   ))==NULL ||
	    (F	    =ALLOC(Float **  ,K    ))==NULL ||
	    (F[0]   =ALLOC(Float *   ,K*Ni ))==NULL ||
	    (F[0][0]=ALLOC(Float     ,K*N2 ))==NULL) return NULL;

	for (m=1; m<Mi; m++) {
	    P[m]=P[m-1]+Ni; Q[m]=Q[m-1]+Ni;
	}
	for (n=1; n<Ni; n++) F[0][n]=F[0][n-1]+Ni;
	for (k=1; k<K ; k++) {
	    Z[k]=Z[k-1]+L2;
	    F[k]=F[k-1]+Ni; for (n=0; n<Ni; n++) F[k][n]=F[k-1][n]+N2;
	}
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

static double	dr,r0,t0,dt;
static int	L,L2,N1,v1,v2;

static THREAD_FUNC	CBP_thread(ki)
THREAD_ARG		ki;
{
	size_t	k=(size_t)ki;
	int	v,h,m,l,n;
	Float	*p,*q,**f=F[k];
	Int	s,c,rv,**fv,rh,*fh,
		**fi=(Int **)f,R0=(Int)floor(ldexp(r0-1.0,R_BITS+S_BITS)+0.5);
	Vector	*z=Z[k];
	double	q0,t,Q0=0.0;
#ifndef	__64BIT_DATABUS
	Int	*qi=(Int *)z+1;
	int	n1;
#else
	struct Vector {
		Int	x,y;
	} *qi=(struct Vector *)z+1;
#endif
	for (v=0; v<N; v++)
	for (h=0; h<N; h++) fi[v][h]=0;

	for (m=k; m<M; m+=K) {
	    for (	 l=0; l<L;	l++)   z[l].x=	    z[l].y=0.0;
	    for (p=P[m], n=0; n<N; n++, l++) { z[l].x=p[n]; z[l].y=0.0; }
	    for (	    ; n<L; n++, l++)   z[l].x=	    z[l].y=0.0;

	    FFT(L2,E,z,-1);

	    for (n=0; n<L2; n++) {
		z[n].x*=G[n]; z[n].y*=G[n];
	    }
	    FFT(L2,E,z, 1);

	    for (q=Q[m], n=0; n<N; n++)
		if ((q0=fabs(q[n]=z[n].x))>Q0) Q0=q0;
	}
	q0=ldexp(1.0/Q0,Q_BITS);

#ifndef	__64BIT_DATABUS
	qi[-1]=qi[N]=0;
#else
	qi[-1].x=qi[N].x=qi[N].y=0;
#endif
	for (m=k; m<M; m+=K) {
	    q=Q[m];
#ifndef	__64BIT_DATABUS
	    for (n=0; n<N; n++) qi[n]=(Int)floor(q[n]*q0+0.5);
#else
	    for (n=N1; n>=0; n--)
		qi[n].y=qi[n+1].x-(qi[n].x=(Int)floor(q[n]*q0+0.5));

	    qi[-1].y=qi[0].x;
#endif
	    t=t0+dt*(double)m; s=(Int)floor(ldexp(sin(t),R_BITS+S_BITS)+0.5);
			       c=(Int)floor(ldexp(cos(t),R_BITS+S_BITS)+0.5);

	    rv=(N1*(s-c)+1)/2-v1*s-R0; fv=fi+v1;
	    for (v=v1; v<=v2; v++, rv-=s, fv++) {
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
	q0=ldexp(dt*Q0/(dr*(double)L2),-(Q_BITS+R_BITS));
	for (v=N1; v>=0; v--)
	for (h=N1; h>=0; h--) f[v][h]=q0*(double)fi[v][h];

	return THREAD_NULL;
}

Float	**CBP(dri,r0i,t0i)
double	dri,r0i,t0i;
{
	double	N12,R,R2,y,x;
	int	v,k,h;
	Float	**Fk,**F0=F[0];

	dr=dri; r0=r0i; t0=t0i;
	L2=(L=POW2(N))<<1;
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
	for (k=1; k<K; k++)
	    if (THREAD_CREATE(T[k],CBP_thread,k))
		Error("can not create CBP_thread.");

	(void)CBP_thread((THREAD_ARG)0);

	for (k=1; k<K; k++) {
	    THREAD_JOIN(T[k]);

	    Fk=F[k]; for (v=0; v<N; v++)
		     for (h=0; h<N; h++) F0[v][h]+=Fk[v][h];
	}
	return F0;
}

void	TermCBP()
{
	free(F[0][0]); free(F[0]); free(F); free(dH); free(H0);
	free(Q[0]); free(Q); free(G); free(Z[0]); free(Z); free(E);
	free(P[0]); free(P);
	free(T);

#ifdef	PTW32_STATIC_LIB
	(void)pthread_win32_process_detach_np();
#endif
}
