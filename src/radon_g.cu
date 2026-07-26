/* radon_g.cu : GPU implementation of RadonSlice() -- see radon.h.
 *
 * Same geometry and sampling as radon_omp.c (unit steps, bilinear).
 * One thread per (view m, detector bin r); the slice stays resident on
 * the device for the whole transform.  O(N^2*M) samples run in tens of
 * milliseconds on a normal GPU, so this is the "radon is GPU-friendly"
 * half of rec2rec_g (the other half is the existing cbp.cu).
 */
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern "C" {
#include "radon.h"
}

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define RADON_BDIM	256

#define RADON_CHECK(call)	do {					\
	cudaError_t e_=(call);						\
	if (e_!=cudaSuccess) {						\
	    fprintf(stderr,"CUDA error %s:%d : %s\n",__FILE__,__LINE__,	\
		    cudaGetErrorString(e_)); exit(1);			\
	}								\
    } while (0)

__global__ static void radon_kernel(const float * __restrict__ img,int N,
				    float *sino,int M,double t0,int L)
{
	int	r=blockIdx.x*blockDim.x+threadIdx.x,
		m=blockIdx.y;
	if (r>=N || m>=M) return;

	double	cx=(double)(N-1)/2.0,
		th=t0+M_PI*(double)m/(double)M,
		c=cos(th),s=sin(th),
		u=(double)r-cx,
		sum=0.0;

	for (int l=-L; l<=L; l++) {
	    double	h=cx+u*c+(double)l*s,
			v=cx-u*s+(double)l*c;

	    if (h<0.0 || v<0.0) continue;
	    int		h0=(int)h,v0=(int)v;
	    if (h0>=N-1 || v0>=N-1) continue;
	    double	fh=h-(double)h0,fv=v-(double)v0;
	    sum+=(1.0-fv)*((1.0-fh)*(double)__ldg(&img[(size_t)v0*N+h0])
			  +     fh *(double)__ldg(&img[(size_t)v0*N+h0+1]))
		+     fv *((1.0-fh)*(double)__ldg(&img[(size_t)(v0+1)*N+h0])
			  +     fh *(double)__ldg(&img[(size_t)(v0+1)*N+h0+1]));
	}
	sino[(size_t)m*N+r]=(float)sum;
}

extern "C" void	RadonSlice(const float *img,int N,float *sino,int M,double t0)
{
	static float	*dI=NULL,*dS=NULL;	/* reused across slices */
	static int	dN=0,dM=0;
	int		L=(int)ceil((double)N*0.70710678118654752)+1;

	if (dI==NULL || dN!=N || dM!=M) {
	    if (dI!=NULL) { cudaFree(dI); cudaFree(dS); }
	    RADON_CHECK(cudaMalloc((void **)&dI,sizeof(float)*(size_t)N*N));
	    RADON_CHECK(cudaMalloc((void **)&dS,sizeof(float)*(size_t)M*N));
	    dN=N; dM=M;
	}
	RADON_CHECK(cudaMemcpy(dI,img,sizeof(float)*(size_t)N*N,
			       cudaMemcpyHostToDevice));

	dim3	grid((N+RADON_BDIM-1)/RADON_BDIM,M);
	radon_kernel<<<grid,RADON_BDIM>>>(dI,N,dS,M,t0,L);
	RADON_CHECK(cudaGetLastError());

	RADON_CHECK(cudaMemcpy(sino,dS,sizeof(float)*(size_t)M*N,
			       cudaMemcpyDeviceToHost));
}
