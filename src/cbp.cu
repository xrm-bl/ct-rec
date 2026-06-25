
#include <stdlib.h>
#include <math.h>
#include <cufft.h>
#include "cu.h"
#include "cbp.h"
#include "cuda13_compat.h"

#ifndef	THREADS_1D
#define THREADS_1D	256
#endif

#ifndef	THREADS_2D
#define THREADS_2D	16
#endif

#ifndef	CUFFT_LIMIT
#define CUFFT_LIMIT	(1U<<23)
#endif

EXTERN double	Ramachandran(int i)
{
	return (i==0)?1.0/4.0:(i&1)?-1.0/(M_PI*M_PI*(double)i*(double)i):0.0;
}

EXTERN double	Shepp(int i)
{
	return 2.0/(M_PI*M_PI*(1.0-4.0*(double)i*(double)i));
}

EXTERN double	Chesler(int i)
{
	if (i==0)
	    return 1.0/8.0-1.0/(2.0*M_PI*M_PI);
	else if (i==1 || i==(-1))
	    return 1.0/16.0-1.0/(2.0*M_PI*M_PI);
	else if (i&1)
	    return -1.0/(2.0*M_PI*M_PI*(double)i*(double)i);
{
	double	i2=(double)i*(double)i,
		i21=i2-1.0;

	return -(i2+1.0)/(2.0*M_PI*M_PI*i21*i21);
}}

__global__ void	GxPQ(int L1,float2 *G,float PI_MdrL2,float2 *PQ)
{
	int	n=blockIdx.x*blockDim.x+threadIdx.x;

	if (n<L1)
{
	PQ+=((size_t)blockIdx.y*(size_t)L1+(size_t)n);
{
	float2	g=G[n],
		pq=(*PQ);

	PQ->x=(pq.x*g.x-pq.y*g.y)*PI_MdrL2;
	PQ->y=(pq.x*g.y+pq.y*g.x)*PI_MdrL2;
}}}

#ifdef	WAI
static int	NAI(int M,double R)
{
	double	O=M_PI*R/(double)M;
	int	o=(int)O;

	return (O>(double)o)?o+1:o;
}

__global__ void	AI(int N,int M,int L12,int O,float r0,float *Q)
{
	int	l,m,o,
		n=blockIdx.x*blockDim.x+threadIdx.x;
	float	q,q0,dq;

	if (n<N)
{
	float	r=(-r0*2.0f-(float)n);

	if (r<0.0f || r>(float)(N-1))
	    q=0.0f;
	else {
	    l=__float2int_rd(r); q=Q[l]+(Q[l+1]-Q[l])*(r-(float)l);
	}
	Q+=n;
{
	size_t	ML12=(size_t)M*(size_t)L12;
	float	*Q1=Q+ML12,
		*Q2=Q+ML12*(size_t)O;

	for (m=0; m<M; m++) {
	    if (m!=0) *(Q2-=L12)=q=q0;

	    dq=__fdividef(q-(q0=(*(Q1-=L12))),(float)O);

	    for (o=1; o<O; o++) *(Q2-=L12)=(q-=dq);
	}
}}}
#endif

__global__ void	PQ2PQ(int N,int M,int L1,float2 *PQ)
{
	int	n,
		m=blockIdx.x*blockDim.x+threadIdx.x;
	float	q,q1;

	if (m<M)
{
	float2	*P=PQ+((size_t)m*(size_t)L1);
	float	*Q=(float *)P;

	q=P[N].x=P[N].y=0.0f;
	for (n=N-1; n>=0; n--) {
	    q1=q; q=P[n].x=Q[n]; P[n].y=q1-q;
	}
}}

/* 改良版 BP_GMF:
 *   旧版はピクセル座標(x,y)を float で M 回 dθ 増分回転していたため、
 *   投影番号 m と半径に比例して角度・ノルムがドリフトし、CPU(double)版との
 *   差分に半月状バイアス＋外周の縞が出ていた。
 *   ここでは投影角の cos/sin を host 側で一度だけ double 計算し、配列 SC[m]
 *   (=(cosθ_m, sinθ_m)) として渡す。各スレッドは固定ピクセル座標(X,Y)から
 *   毎投影 r = X cosθ_m + Y sinθ_m - r0 を直接求める(累積誤差ゼロ)。
 *   SC[m] はワープ内全スレッドで同一アドレス参照になりブロードキャストされる。 */
__global__ void	BP_GMF(int N,int M,int L1,
		       float xy0,float R2,float r0,
		       const float2 * __restrict__ SC,
		       const float2 * __restrict__ Q,
		       float *F)
{
	ENABLE_SMEM_SPILLING();
	int	h,v,m,n;

	if ((h=blockIdx.x*blockDim.x+threadIdx.x)<N &&
	    (v=blockIdx.y*blockDim.y+threadIdx.y)<N)
{
	float	X=(float)h-xy0,
		Y=xy0-(float)v,
		f=0.0f;

	if (X*X+Y*Y<=R2) {
	    for (m=0; m<M; m++, Q+=L1) {
		float2	sc=__ldg(&SC[m]);		/* (cosθ_m, sinθ_m) */
		float	r=fmaf(Y,sc.y,fmaf(X,sc.x,-r0));/* X cosθ + Y sinθ - r0 */

		n=__float2int_rd(r);			/* floor (負側も正しい) */
		if ((unsigned)n<(unsigned)N) {		/* 上下の範囲外を保護 */
		    float2 q=Q[n];
		    f+=fmaf(q.y,(r-(float)n),q.x);	/* q.x + q.y*(r-n) */
		}
	    }
	}
	F[(size_t)v*(size_t)N+(size_t)h]=f;
}}

/* NOTE: Legacy texture reference BP kernel removed for CUDA 13.x compatibility.
 *       BP_GMF (global memory fetch) kernel is used for all cases. */

static int		N,M,L,L1,L2,batch;
static Float		**p,**f;
static size_t		sof_L2,sof2_L1,sof_N,MOmax;
static float		*gpf,*F;
static float2		*G,*PQ,*B,*SC,*scf;	/* SC:角度表(device) scf:host転送用 */
static cufftHandle	R2C,C2R;

EXTERN Float	**InitCBP(int n,int m)
{
	SETUP_CUDA_GPU();

	N=n; M=m; for (L=1; L<N; L*=2) ; L1=L+1; L2=L*2;

#define ALLOC(type,noe)	(type *)malloc(sizeof(type)*noe)

	if ((p =ALLOC(Float *,(size_t)M		 ))==NULL ||
	    (*p=ALLOC(Float  ,(size_t)M*(size_t)N))==NULL ||
	    (f =ALLOC(Float *,(size_t)N		 ))==NULL ||
	    (*f=ALLOC(Float  ,(size_t)N*(size_t)N))==NULL) return NULL;

	for (m=1; m<M; m++) p[m]=p[m-1]+N;
	for (n=1; n<N; n++) f[n]=f[n-1]+N;

	sof_L2=sizeof(float)*(size_t)L2;
	sof2_L1=sizeof(float2)*(size_t)L1;
	sof_N=sizeof(float)*(size_t)N;

	CUDA_SAFE_CALL(cudaMallocHost((void **)&gpf,sof_L2));

	CUDA_SAFE_CALL(cudaMalloc((void **)&G ,sof2_L1));
#ifdef	WAI
	CUDA_SAFE_CALL(cudaMalloc((void **)&PQ,sof2_L1*(size_t)M
					      *(size_t)NAI(M,(double)(N-1)/2.0)
		      )		 );
#else
	CUDA_SAFE_CALL(cudaMalloc((void **)&PQ,sof2_L1*(size_t)M));
#endif
	CUDA_SAFE_CALL(cudaMalloc((void **)&F, sof_N*(size_t)N));

	/* 投影角表 SC[m]=(cosθ_m,sinθ_m) 用バッファ。角度オーバーサンプリング
	 * (WAI) を使う場合は MO=M*O まで増えるので PQ と同じ最大数で確保する。 */
#ifdef	WAI
	MOmax=(size_t)M*(size_t)NAI(M,(double)(N-1)/2.0);
#else
	MOmax=(size_t)M;
#endif
	CUDA_SAFE_CALL(cudaMalloc((void **)&SC,sizeof(float2)*MOmax));
	CUDA_SAFE_CALL(cudaMallocHost((void **)&scf,sizeof(float2)*MOmax));

	for (n=0; n<L2; n++) gpf[n]=Filter(n-L);

	CUDA_SAFE_CALL(cudaMemcpy(G,gpf,sof_L2,cudaMemcpyHostToDevice));

	CUFFT_SAFE_CALL(cufftPlan1d(&R2C,L2,CUFFT_R2C,1));
	CUFFT_SAFE_CALL(cufftExecR2C(R2C,(cufftReal *)G,(cufftComplex *)G));

	if (batch=(size_t)L2*(size_t)M<=CUFFT_LIMIT) {
	    CUFFT_SAFE_CALL(cufftDestroy(R2C));
	    CUFFT_SAFE_CALL(cufftPlan1d(&R2C,L2,CUFFT_R2C,M));
	    CUFFT_SAFE_CALL(cufftPlan1d(&C2R,L2,CUFFT_C2R,M));
	}
	else {
	    CUFFT_SAFE_CALL(cufftPlan1d(&C2R,L2,CUFFT_C2R,1));
	    CUDA_SAFE_CALL(cudaMalloc((void **)&B,sof2_L1));
	}
	return p;
}

EXTERN void	PrepareCBP()
{
	int	n,m;

	for (n=0; n<L; n++) gpf[n]=0.0;
	for (n=L+N; n<L2; n++) gpf[n]=0.0;

	for (m=0; m<M; m++) {
	    for (n=0; n<N; n++) gpf[L+n]=p[m][n];

	    CUDA_SAFE_CALL(cudaMemcpy(PQ+(size_t)m*(size_t)L1,gpf,sof_L2,
				      cudaMemcpyHostToDevice));
	}
}

static const char	*kef="kernel execution failed.";

#define KEF()	CUT_CHECK_ERROR(kef)

EXTERN void	ExecuteCBP(double dr,double r0,double t0)
{
	float2	*pq;
	int	m;
	double	N1=(double)(N-1),
		R=(2.0*r0+N1>0.0)?-r0:r0+N1;
#ifdef	WAI
	int	O=NAI(M,R),
		MO=M*O;
#else
#define MO	M
#endif
	if (batch) {
	    CUFFT_SAFE_CALL(cufftExecR2C(R2C,(cufftReal *)PQ,
					     (cufftComplex *)PQ));
	}
	else
	    for (pq=PQ, m=0; m<M; m++, pq+=L1) {
		CUDA_SAFE_CALL(cudaMemcpy(B,pq,sof2_L1,
					  cudaMemcpyDeviceToDevice));
		CUFFT_SAFE_CALL(cufftExecR2C(R2C,(cufftReal *)B,
						 (cufftComplex *)B));
		CUDA_SAFE_CALL(cudaMemcpy(pq,B,sof2_L1,
					  cudaMemcpyDeviceToDevice));
	    }
{
	dim3	blocks_1d((L1+THREADS_1D-1)/THREADS_1D,M);

	GxPQ<<<blocks_1d,THREADS_1D>>>
	    (L1,G,(float)(M_PI/((double)MO*dr*(double)L2)),PQ); KEF();
}
	if (batch) {
	    CUFFT_SAFE_CALL(cufftExecC2R(C2R,(cufftComplex *)PQ,
					     (cufftReal *)PQ));
	}
	else
	    for (pq=PQ, m=0; m<M; m++, pq+=L1) {
		CUDA_SAFE_CALL(cudaMemcpy(B,pq,sof2_L1,
					  cudaMemcpyDeviceToDevice));
		CUFFT_SAFE_CALL(cufftExecC2R(C2R,(cufftComplex *)B,
						 (cufftReal *)B));
		CUDA_SAFE_CALL(cudaMemcpy(pq,B,sof2_L1,
					  cudaMemcpyDeviceToDevice));
	    }

#ifdef	WAI
	if (O>1) {
	    AI<<<(N+THREADS_1D-1)/THREADS_1D,THREADS_1D>>>
	      (N,M,L1*2,O,(float)r0,(float *)PQ); KEF();
	}
#endif
	PQ2PQ<<<(MO+THREADS_1D-1)/THREADS_1D,THREADS_1D>>>(N,MO,L1,PQ); KEF();
{
	int	blocks=(N+THREADS_2D-1)/THREADS_2D;
	dim3	blocks_2d(blocks,blocks),
		threads_2d(THREADS_2D,THREADS_2D);
	double	dt=M_PI/(double)MO;

	/* 投影角の cos/sin を double で一度だけ生成して転送(ドリフト除去の要) */
	for (m=0; m<MO; m++) {
	    double th=t0+dt*(double)m;
	    scf[m].x=(float)cos(th); scf[m].y=(float)sin(th);
	}
	CUDA_SAFE_CALL(cudaMemcpy(SC,scf,sizeof(float2)*(size_t)MO,
				  cudaMemcpyHostToDevice));

	    BP_GMF<<<blocks_2d,threads_2d>>>
		  (N,MO,L1,
		   (float)(N1/2.0),(float)(R*R),
		   (float)r0,
		   SC,
		   PQ,
		   F);
}}

EXTERN void	BeginCBP(double dr,double r0,double t0)
{
	PrepareCBP(); ExecuteCBP(dr,r0,t0);
}

EXTERN Float	**EndCBP()
{
	int	y,x;

	KEF();

	for (y=0; y<N; y++) {
	    CUDA_SAFE_CALL(cudaMemcpy(gpf,F+(size_t)y*(size_t)N,sof_N,
				      cudaMemcpyDeviceToHost));

	    for (x=0; x<N; x++) f[y][x]=gpf[x];
	}
	return f;
}

EXTERN Float	**CBP(double dr,double r0,double t0)
{
	BeginCBP(dr,r0,t0); return EndCBP();
}

EXTERN void	TermCBP()
{
	if (!batch) {
	    CUDA_SAFE_CALL(cudaFree(B));
	}
	CUFFT_SAFE_CALL(cufftDestroy(C2R));
	CUFFT_SAFE_CALL(cufftDestroy(R2C));

	CUDA_SAFE_CALL(cudaFree(F));
	CUDA_SAFE_CALL(cudaFree(PQ));
	CUDA_SAFE_CALL(cudaFree(G));
	CUDA_SAFE_CALL(cudaFree(SC));

	CUDA_SAFE_CALL(cudaFreeHost(scf));
	CUDA_SAFE_CALL(cudaFreeHost(gpf));

	free(*f); free(f); free(*p); free(p);
}
