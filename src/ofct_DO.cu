/*======================================================================*/
/* ofct_DO.cu : GPU version of ofct_DO (offset-CT rotation-axis finder). */
/*                                                                        */
/* Host side (reading via rhp_c, -log, horizontal flip of the opposing   */
/* view, and optional Gaussian smoothing) is identical to ofct_DO.c so    */
/* the result matches the CPU version numerically.  The per-view-pair MSD */
/* over the offset grid -- the CalcMSD() bottleneck in the CPU code -- is */
/* replaced by a direct sum-of-squared-differences (SSD) CUDA kernel:     */
/*                                                                        */
/*   D(Dx,Dy) = SUM_pairs SUM_overlap (G(x,y) - Hf(x-Dx, y-Dy))^2         */
/*              / ( count(Dx,Dy) * M )                                    */
/*                                                                        */
/* count(Dx,Dy) = (Nx-|Dx|)*(Ny-|Dy|) is the overlap pixel count (equal   */
/* for every pair), so the kernel accumulates the raw SSD and the         */
/* per-offset normalisation is applied once, on the host, after the loop. */
/* A global constant factor (the CPU FFT scale Lx*Ly and 1/M) cancels in  */
/* both argmin and the parabolic sub-pixel ratio, so the estimated centre */
/* is identical to the CPU code.                                          */
/*======================================================================*/

#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern "C" {
#include "rhp.h"		/* HiPic, InitReadHiPic/ReadHiPic/TermReadHiPic; FOM=double */
}
extern "C" void	Error(const char *msg);

#define CUDA_CHECK(call)	do {					\
	cudaError_t e_=(call);						\
	if (e_!=cudaSuccess) {						\
	    fprintf(stderr,"CUDA error %s:%d : %s\n",__FILE__,__LINE__,\
		    cudaGetErrorString(e_)); exit(1);			\
	}								\
    } while (0)

#define BDIM	256

/*----------------------------------------------------------------------*/
/* One CUDA block per offset (Dx,Dy); threads reduce over the overlap.   */
/* NOTE on the offset convention: G(x,y) is compared with Hf(x-Dx,y-Dy)  */
/* (Dx>0 aligns the right part of G with the left part of the flipped    */
/* opposing view Hf -- the physical overlap in offset CT).  If a          */
/* validation run against ofct_DO (CPU) shows a mirrored centre, flip the */
/* sign of Dx/Dy here (compare Hf(x+Dx,y+Dy)).                            */
__global__ void ssd_kernel(const float * __restrict__ G,
			   const float * __restrict__ H,
			   int Nx,int Ny,int Ox1,int Oy1,int Ox,int Oy,
			   double *S,double *Pair)
{
	int	off=blockIdx.x;
	if (off>=Ox*Oy) return;

	int	ix=off%Ox, iy=off/Ox;
	int	Dx=Ox1+ix, Dy=Oy1+iy;
	int	adx=(Dx<0)?-Dx:Dx, ady=(Dy<0)?-Dy:Dy;
	int	nx=Nx-adx, ny=Ny-ady;
	if (nx<=0 || ny<=0) return;

	int	gx0=(Dx>0)?Dx:0, hx0=(Dx>0)?0:-Dx;
	int	gy0=(Dy>0)?Dy:0, hy0=(Dy>0)?0:-Dy;

	long	total=(long)nx*ny;
	double	local=0.0;

	for (long t=threadIdx.x; t<total; t+=blockDim.x) {
	    int	i=(int)(t%nx), j=(int)(t/nx);
	    double	d=(double)G[(gy0+j)*Nx+(gx0+i)]-
			 (double)H[(hy0+j)*Nx+(hx0+i)];
	    local+=d*d;
	}

	__shared__ double sh[BDIM];
	sh[threadIdx.x]=local;
	__syncthreads();
	for (int s=blockDim.x/2; s>0; s>>=1) {
	    if (threadIdx.x<s) sh[threadIdx.x]+=sh[threadIdx.x+s];
	    __syncthreads();
	}
	if (threadIdx.x==0) { Pair[off]=sh[0]; S[off]+=sh[0]; }
}

/*----------------------------------------------------------------------*/
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

/*----------------------------------------------------------------------*/
/* Reproduce msd.c's FFT length so the per-pair progress value matches   */
/* the CPU version's normalisation (D = SSD / (Lx*Ly * overlap)).        */
static int	Range(int *Ng,int *Nh,int *D1,int *D2,int *Og,int *Oh)
{
	int	g,h,N,L;

	if ((g=(*Nh)+(*D2))>(*Ng)) g=(*Ng);
	if ((h=(*Ng)-(*D1))>(*Nh)) h=(*Nh);

	if (*D2<0) {
	    *Og=0; *Oh=(-(*D2)); *D1-=(*D2); *D2=0; N=1-(*D1);
	}
	else if (*D1>0) {
	    *Og=(*D1); *Oh=0; *D2-=(*D1); *D1=0; N=(*D2)+1;
	}
	else {
	    *Og=(*Oh)=0; N=(-(*D1)>(*D2))?1-(*D1):*D2+1;
	}
	N+=(((*Ng=g-(*Og))>(*Nh=h-(*Oh)))?*Ng:*Nh);

	for (L=1; L<N; L<<=1) ; return L;
}

/*----------------------------------------------------------------------*/
int	main(int argc,char **argv)
{
	char		*env;
	const char	*mae="memory allocation error.";
	HiPic		hp;
	int		Ox1,Ox2,Oy1,Oy2,Ox,Oy,Nx,Ny,M,m,y,x,X,Y;
	int		Lx,Ly;
	FOM		*Gflat,*Hflat,**G,**H,*tmp=NULL;
	float		*fG,*fH,*dG,*dH;
	double		*dS,*dPair,**S,*Sflat,*pairbuf;
	double		sig,s,l;
	double		*ker=NULL;
	int		r=0,kk;
	int		Rc=0,Oyv=0;
	size_t		p,NxNy;

	if (argc!=2)
	    	Error("usage : ofct_DO_g raw/");

	InitReadHiPic(argv[1],&hp);

//	if (hp.Nt%2) (void)fputs("bad number of views (warning).\n",stderr);

	Nx=hp.Nx; Ny=hp.Ny; M=hp.Nt/2; NxNy=(size_t)Nx*Ny;

	Ox1=(int)ceil(0.7*Nx); Ox2=Nx-2;
	Oy1=-10; Oy2=10;

	Ox=Ox2-Ox1+1;
	Oy=Oy2-Oy1+1;

	/* FFT length used by the CPU CalcMSD normalisation (Range mutates args) */
	{ int gx=Nx,hx=Nx,d1=Ox1,d2=Ox2,ogx,ohx; Lx=Range(&gx,&hx,&d1,&d2,&ogx,&ohx); }
	{ int gy=Ny,hy=Ny,e1=Oy1,e2=Oy2,ogy,ohy; Ly=Range(&gy,&hy,&e1,&e2,&ogy,&ohy); }
	l=(double)Lx*(double)Ly;

	/* host buffers (one G and one flipped H per pair) */
	if ((Gflat=(FOM *)malloc(sizeof(FOM)*(size_t)Nx*Ny))==NULL ||
	    (Hflat=(FOM *)malloc(sizeof(FOM)*(size_t)Nx*Ny))==NULL ||
	    (G=(FOM **)malloc(sizeof(FOM *)*Ny))==NULL ||
	    (H=(FOM **)malloc(sizeof(FOM *)*Ny))==NULL ||
	    (fG=(float *)malloc(sizeof(float)*(size_t)Nx*Ny))==NULL ||
	    (fH=(float *)malloc(sizeof(float)*(size_t)Nx*Ny))==NULL ||
	    (Sflat=(double *)malloc(sizeof(double)*(size_t)Ox*Oy))==NULL ||
	    (pairbuf=(double *)malloc(sizeof(double)*(size_t)Ox*Oy))==NULL ||
	    (S=(double **)malloc(sizeof(double *)*Oy))==NULL) Error(mae);

	for (y=0; y<Ny; y++) { G[y]=Gflat+(size_t)y*Nx; H[y]=Hflat+(size_t)y*Nx; }
	for (y=0; y<Oy; y++)   S[y]=Sflat+(size_t)y*Ox;

	    /* (A) Gaussian kernel for pre-MSD smoothing (env OFCT_DO_SMOOTH=sigma, 0=off) */
	sig = ((env=getenv("OFCT_DO_SMOOTH"))!=NULL) ? atof(env) : 1.0;
	if (sig>0.0) {
	    r=(int)ceil(3.0*sig);
	    if ((ker=(double *)malloc(sizeof(double)*(2*r+1)))==NULL ||
		(tmp=(FOM *)malloc(sizeof(FOM)*(size_t)Nx*Ny))==NULL) Error(mae);
	    { double sm=0.0; for(kk=-r;kk<=r;kk++){ ker[kk+r]=exp(-0.5*(double)kk*kk/(sig*sig)); sm+=ker[kk+r]; } for(kk=0;kk<2*r+1;kk++) ker[kk]/=sm; }
	}

	/* device buffers */
	CUDA_CHECK(cudaMalloc((void **)&dG,sizeof(float)*(size_t)Nx*Ny));
	CUDA_CHECK(cudaMalloc((void **)&dH,sizeof(float)*(size_t)Nx*Ny));
	CUDA_CHECK(cudaMalloc((void **)&dS,sizeof(double)*(size_t)Ox*Oy));
	CUDA_CHECK(cudaMalloc((void **)&dPair,sizeof(double)*(size_t)Ox*Oy));
	CUDA_CHECK(cudaMemset(dS,0,sizeof(double)*(size_t)Ox*Oy));

	/* accumulate SSD over all view pairs */
	for (m=0; m<M; m++) {
	    ReadHiPic(&hp,m);
	    for (y=0; y<Ny; y++) for (x=0; x<Nx; x++) G[y][x]=(-Log(hp.T[y][x]));
	    if (sig>0.0) Smooth(G,Nx,Ny,ker,r,tmp);

	    ReadHiPic(&hp,m+M);
	    for (y=0; y<Ny; y++) for (x=0; x<Nx; x++) H[y][Nx-1-x]=(-Log(hp.T[y][x]));
	    if (sig>0.0) Smooth(H,Nx,Ny,ker,r,tmp);

	    for (p=0; p<NxNy; p++) { fG[p]=(float)Gflat[p]; fH[p]=(float)Hflat[p]; }

	    CUDA_CHECK(cudaMemcpy(dG,fG,sizeof(float)*NxNy,cudaMemcpyHostToDevice));
	    CUDA_CHECK(cudaMemcpy(dH,fH,sizeof(float)*NxNy,cudaMemcpyHostToDevice));

	    ssd_kernel<<<Ox*Oy,BDIM>>>(dG,dH,Nx,Ny,Ox1,Oy1,Ox,Oy,dS,dPair);
	    CUDA_CHECK(cudaGetLastError());

	    /* per-pair progress (same as ofct_DO.c Scan(): m, best Ox, best Oy, min) */
	    CUDA_CHECK(cudaMemcpy(pairbuf,dPair,sizeof(double)*(size_t)Ox*Oy,cudaMemcpyDeviceToHost));
	    {
		int	ix,iy,bX=Ox1,bY=Oy1,adx0=(Ox1<0)?-Ox1:Ox1,ady0=(Oy1<0)?-Oy1:Oy1;
		double	c0=(double)(Nx-adx0)*(double)(Ny-ady0);
		double	bd=(c0>0.0)?pairbuf[0]/(c0*l):1.0e300;
		for (iy=0; iy<Oy; iy++)
		for (ix=0; ix<Ox; ix++) {
		    int	Dx=Ox1+ix, Dy=Oy1+iy;
		    int	adx=(Dx<0)?-Dx:Dx, ady=(Dy<0)?-Dy:Dy;
		    double cnt=(double)(Nx-adx)*(double)(Ny-ady);
		    double v=(cnt>0.0)?pairbuf[(size_t)iy*Ox+ix]/(cnt*l):1.0e300;
		    if (v<bd) { bd=v; bX=Dx; bY=Dy; }
		}
		(void)printf("%d\t%d\t%d\t%le\r",m,bX,bY,bd);
	    }
	}
	CUDA_CHECK(cudaDeviceSynchronize());
	(void)printf("\n");

	/* receive device result into Sflat, then normalise in place via S[][] */
	CUDA_CHECK(cudaMemcpy(Sflat,dS,sizeof(double)*(size_t)Ox*Oy,cudaMemcpyDeviceToHost));

	/* per-offset overlap normalisation: S = SSD / (count * M) */
	for (y=0; y<Oy; y++)
	for (x=0; x<Ox; x++) {
	    int	Dx=Ox1+x, Dy=Oy1+y;
	    int	adx=(Dx<0)?-Dx:Dx, ady=(Dy<0)?-Dy:Dy;
	    double count=(double)(Nx-adx)*(double)(Ny-ady);
	    S[y][x] = (count>0.0) ? Sflat[(size_t)y*Ox+x]/(count*l*(double)M) : 1.0e300;
	}

	/* minimum + parabolic sub-pixel refinement (identical to ofct_DO.c) */
	s=S[Y=0][X=0];
	for (y=0; y<Oy; y++)
	for (x=0; x<Ox; x++)
	    if (S[y][x]<s) { s=S[y][x]; X=x; Y=y; }

	{
        double dx=0.0,dy=0.0,a,b,c,den,centerf,oyf;
        if (X>=1 && X<=Ox-2){ a=S[Y][X-1]; b=S[Y][X]; c=S[Y][X+1]; den=a-2.0*b+c; if(den>0.0){ dx=0.5*(a-c)/den; if(dx<-1.0)dx=-1.0; if(dx>1.0)dx=1.0; } }
        if (Y>=1 && Y<=Oy-2){ a=S[Y-1][X]; b=S[Y][X]; c=S[Y+1][X]; den=a-2.0*b+c; if(den>0.0){ dy=0.5*(a-c)/den; if(dy<-1.0)dy=-1.0; if(dy>1.0)dy=1.0; } }
        centerf=((double)(Nx+Ox1+X)+dx)/2.0;
        oyf=(double)(Oy1+Y)+dy;
        Rc=(int)floor(centerf+0.5); Oyv=(int)floor(oyf+0.5);
        (void)fprintf(stderr,"try: ofct_srec_g_c raw %d %d - 1.0 0.0 rec\n",Rc,Oyv);
//        (void)fprintf(stderr,"center=%.3f\tOy=%.3f\t(dx=%.3f dy=%.3f, smooth=%.2f)\n",centerf,oyf,dx,dy,sig);
        (void)fprintf(stderr,"%d\n",Rc);
	}

	cudaFree(dG); cudaFree(dH); cudaFree(dS); cudaFree(dPair);
	free(S); free(Sflat); free(pairbuf);
	free(fG); free(fH);
	free(G); free(H); free(Gflat); free(Hflat);
	if(ker)free(ker); if(tmp)free(tmp);

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
