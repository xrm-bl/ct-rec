
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rhp.h"
#include "msd.h"
#include "sif_f.h"

extern void	Error(char *msg);

#ifdef	_WIN32
#include <process.h>
#include <windows.h>

#define THREAD_T	HANDLE
#define FUNCTION_T	unsigned __stdcall
#define RETURN_VALUE	0

#define INIT_THREAD(T,F,A)	\
	(T=(HANDLE)_beginthreadex(NULL,0,F,(void *)(A),0,NULL))==0
#define TERM_THREAD(T)	\
	WaitForSingleObject(T,INFINITE)==WAIT_FAILED || CloseHandle(T)==0
#else
#include <pthread.h>

#define THREAD_T	pthread_t
#define FUNCTION_T	void *
#define RETURN_VALUE	NULL

#define INIT_THREAD(T,F,A)	pthread_create(&(T),NULL,F,(void *)(A))
#define TERM_THREAD(T)		pthread_join(T,NULL)
#endif

#define INIT_MT(T,F,A)	\
if (INIT_THREAD(T,F,A)) Error(#F " : multi-threading initialization error.")

#define TERM_MT(T,F)	\
	if (TERM_THREAD(T)) Error(#F " : multi-threading termination error.")

static int	Threads=40,
		Swapped=0;
static MSD	*msd;
static FOM	***G,***H;

static FUNCTION_T	Compare(void *i)
{
	size_t	j=(size_t)i,
		k=j+Threads*(1-Swapped);

	CalcMSD(msd+j,G[k],H[k]); return RETURN_VALUE;
}

static void	Scan(int Ox1,int Oy1,int Ox2,int Oy2,FOM *D,FOM *S,int m)
{
	int	y,x,
		X=Ox1,
		Y=Oy1;
	double	d=(*D);

	for (y=Oy1; y<=Oy2; y++)
	for (x=Ox1; x<=Ox2; x++) {
	    if (*D<d) {
		d=(*D); X=x; Y=y;
	    }
	    *(S++)+=(*(D++));
	}
	(void)printf("%d\t%d\t%d\t%le\r",m,X,Y,d);
}

static double	Log(double d)
{
	return (d>0.0)?log(d):0.0;
}

/*----------------------------------------------------------------------*/
/* (A) separable Gaussian smoothing of a view (noise robustness).        */
static void Smooth(FOM **F,int Nx,int Ny,double *ker,int r,FOM *tmp)
{
    int    x,y,k,xx,yy;
    double s,w;

    for (y=0; y<Ny; y++)
    for (x=0; x<Nx; x++) {
        s=w=0.0;
        for (k=-r; k<=r; k++){ xx=x+k; if(xx>=0 && xx<Nx){ s+=F[y][xx]*ker[k+r]; w+=ker[k+r]; } }
        tmp[y*Nx+x]=(FOM)(s/w);
    }
    for (y=0; y<Ny; y++) for (x=0; x<Nx; x++) F[y][x]=tmp[y*Nx+x];

    for (x=0; x<Nx; x++)
    for (y=0; y<Ny; y++) {
        s=w=0.0;
        for (k=-r; k<=r; k++){ yy=y+k; if(yy>=0 && yy<Ny){ s+=F[yy][x]*ker[k+r]; w+=ker[k+r]; } }
        tmp[y*Nx+x]=(FOM)(s/w);
    }
    for (y=0; y<Ny; y++) for (x=0; x<Nx; x++) F[y][x]=tmp[y*Nx+x];
}
int	main(int argc,char **argv)
{
	char		*env,
			*mae="memory allocation error.";
	HiPic		hp;
	int		Ox1,Ox2,Oy1,Oy2,Ox,Oy,T,t,y,x,M,m,t0,m0,Y,X;
	THREAD_T	*Thread;
	FOM		**S,**F;
	double		sig;
	double		*ker=NULL;
	FOM		*tmp=NULL;
	int		r=0,kk;
	int		Rc=0,Oyv=0;
	size_t		i;
	double		s;

	if ((argc!=2))
	    	Error("usage : ofct_DO raw/");

	if ((env=getenv("THREADS"))!=NULL &&
	    (Threads=atoi(env))<=0) Error("bad number of THREADS.");

	InitReadHiPic(argv[1],&hp);

	if (hp.Nt%2) (void)fputs("bad number of views (warning).\n",stderr);

		if (argc==6 || argc==7){
			
			if ((Ox1=atoi(argv[2]))<=(-hp.Nx) ||
			    (Ox2=atoi(argv[3]))>=  hp.Nx  || Ox1>Ox2 ||
			    (Oy1=atoi(argv[4]))<=(-hp.Ny) ||
			    (Oy2=atoi(argv[5]))>=  hp.Ny  || Oy1>Oy2)
			    Error("bad range of offsets.");
		} else{
			Ox1=ceil(0.7*hp.Nx); Ox2=hp.Nx-2;
			Oy1=-10; Oy2=10;
		}

	Ox=Ox2-Ox1+1;
	Oy=Oy2-Oy1+1;
		
	T=Threads*2;

	if ((Thread=(THREAD_T *)malloc(sizeof(THREAD_T)*Threads))==NULL ||
	    (msd=(MSD *)malloc(sizeof(MSD)*Threads))==NULL ||
	    (G=(FOM ***)malloc(sizeof(FOM **)*T))==NULL ||
	    (*G=(FOM **)malloc(sizeof(FOM *)*T*hp.Ny))==NULL ||
	    (**G=(FOM *)malloc(sizeof(FOM)*T*hp.Ny*hp.Nx))==NULL ||
	    (H=(FOM ***)malloc(sizeof(FOM **)*T))==NULL ||
	    (*H=(FOM **)malloc(sizeof(FOM *)*T*hp.Ny))==NULL ||
	    (**H=(FOM *)malloc(sizeof(FOM)*T*hp.Ny*hp.Nx))==NULL ||
	    (S=(FOM **)malloc(sizeof(FOM *)*Oy))==NULL ||
	    (*S=(FOM *)malloc(sizeof(FOM)*Oy*Ox))==NULL) Error(mae);

	for (t=0; t<Threads; t++)
	    if (!InitMSD(msd+t,hp.Nx,hp.Ny,
			       hp.Nx,hp.Ny,Ox1,Oy1,Ox2,Oy2)) Error(mae);

	for (t=0; t<T; t++) {
	    if (t>0) {
		G[t]=G[t-1]+hp.Ny; H[t]=H[t-1]+hp.Ny;
	    }
	    for (y=0; y<hp.Ny; y++)
		if (t>0 || y>0) {
		    G[t][y]=G[t][y-1]+hp.Nx; H[t][y]=H[t][y-1]+hp.Nx;
		}
	}
	for (y=0; y<Oy; y++) {
	    if (y>0) S[y]=S[y-1]+Ox;

	    for (x=0; x<Ox; x++) S[y][x]=0.0;
	}
	    /* (A) Gaussian kernel for pre-MSD smoothing (env OFCT_DO_SMOOTH=sigma, 0=off) */
    sig = ((env=getenv("OFCT_DO_SMOOTH"))!=NULL) ? atof(env) : 1.0;
    if (sig>0.0) {
        r=(int)ceil(3.0*sig);
        if ((ker=(double *)malloc(sizeof(double)*(2*r+1)))==NULL ||
            (tmp=(FOM *)malloc(sizeof(FOM)*(size_t)hp.Nx*hp.Ny))==NULL) Error(mae);
        { double sm=0.0; for(kk=-r;kk<=r;kk++){ ker[kk+r]=exp(-0.5*(double)kk*kk/(sig*sig)); sm+=ker[kk+r]; } for(kk=0;kk<2*r+1;kk++) ker[kk]/=sm; }
    }
	M=hp.Nt/2;

	for (m=0; m<M; m+=Threads) {
	    T=(m+Threads<=M)?Threads:M-m;

	    t0=Threads*Swapped;
	    for (t=0; t<T; t++) {
		ReadHiPic(&hp,m+t);
		F=G[t+t0];
		for (y=0; y<hp.Ny; y++)
		for (x=0; x<hp.Nx; x++) F[y][x]=(-Log(hp.T[y][x]));

			if (sig>0.0) Smooth(G[t+t0],hp.Nx,hp.Ny,ker,r,tmp);

		ReadHiPic(&hp,m+t+M);
		F=H[t+t0];
		for (y=0; y<hp.Ny; y++)
		for (x=0; x<hp.Nx; x++) F[y][hp.Nx-1-x]=(-Log(hp.T[y][x]));

			if (sig>0.0) Smooth(H[t+t0],hp.Nx,hp.Ny,ker,r,tmp);
	    }
	    if ((m0=m-Threads)>=0)
		for (t=0; t<Threads; t++) {
		    TERM_MT(Thread[t],Compare);

		    Scan(Ox1,Oy1,Ox2,Oy2,msd[t].D,*S,m0+t);
		}

	    Swapped=1-Swapped; for (i=0; i<T; i++) INIT_MT(Thread[i],Compare,i);
	}
	T=M-(m0=m-Threads);
	for (t=0; t<T; t++) {
	    TERM_MT(Thread[t],Compare);

	    Scan(Ox1,Oy1,Ox2,Oy2,msd[t].D,*S,m0+t);
	}
	(void)printf("\n");
	for (t=0; t<Threads; t++) TermMSD(msd+t);

	s=S[Y=0][X=0];
	for (y=0; y<Oy; y++)
	for (x=0; x<Ox; x++) {
	    if (S[y][x]<s) {
		s=S[y][x]; X=x; Y=y;
	    }
	    S[y][x]/=(double)M;
	}
//	(void)fprintf(stderr,"%d\t%d\t%le\n",+X+Ox1,Y+Oy1,s/(double)M);
//	(void)fprintf(stderr,"%d\t%d\t%le\n",(hp.Nx+X+Ox1)/2,Y+Oy1,s/(double)M);
	{
        double dx=0.0,dy=0.0,a,b,c,den,centerf,oyf;
        if (X>=1 && X<=Ox-2){ a=S[Y][X-1]; b=S[Y][X]; c=S[Y][X+1]; den=a-2.0*b+c; if(den>0.0){ dx=0.5*(a-c)/den; if(dx<-1.0)dx=-1.0; if(dx>1.0)dx=1.0; } }
        if (Y>=1 && Y<=Oy-2){ a=S[Y-1][X]; b=S[Y][X]; c=S[Y+1][X]; den=a-2.0*b+c; if(den>0.0){ dy=0.5*(a-c)/den; if(dy<-1.0)dy=-1.0; if(dy>1.0)dy=1.0; } }
        centerf=((double)(hp.Nx+Ox1+X)+dx)/2.0;
        oyf=(double)(Oy1+Y)+dy;
        Rc=(int)floor(centerf+0.5); Oyv=(int)floor(oyf+0.5);
        (void)fprintf(stderr,"try: ofct_srec_g_c raw %d %d - 1.0 0.0 rec\n",Rc,Oyv);
//        (void)fprintf(stderr,"center=%.3f\tOy=%.3f\t(dx=%.3f dy=%.3f, smooth=%.2f)\n",centerf,oyf,dx,dy,sig);
        (void)fprintf(stderr,"%d\n",Rc);
    }


	free(*S); free(S);
	free(**H); free(*H); free(H);
	free(**G); free(*G); free(G);
	if(ker)free(ker); if(tmp)free(tmp);
	free(msd);
	free(Thread);

	TermReadHiPic(&hp);

	    /* append the executed command to cmd-hst.log in the working directory */
    {
        FILE *logf;
        int  ai;
        if ((logf=fopen("cmd-hst.log","a"))!=NULL){
            for (ai=0; ai<argc; ++ai) fprintf(logf,"%s ",argv[ai]);
            fprintf(logf,"   %% center %d Oy %d smooth %.2f\n",Rc,Oyv,sig);
            fclose(logf);
        }
    }
	return 0;
}
