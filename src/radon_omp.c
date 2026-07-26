/* radon_omp.c : CPU (OpenMP) implementation of RadonSlice() -- see radon.h.
 *
 * For each view m and detector bin r the line integral is taken with unit
 * (1 px) steps along the ray, sampling the slice by bilinear interpolation.
 * Points outside the image contribute 0.  Views are independent, so the
 * outer loop is parallelised over m.
 */
#include <math.h>
#include <stddef.h>
#include "radon.h"

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

void	RadonSlice(const float *img,int N,float *sino,int M,double t0)
{
	double	cx=(double)(N-1)/2.0;
	int	L=(int)ceil((double)N*0.70710678118654752)+1,	/* half diagonal */
		m;

	#pragma omp parallel for schedule(static)
	for (m=0; m<M; m++) {
	    double	th=t0+M_PI*(double)m/(double)M,
			c=cos(th),s=sin(th);
	    int		r,l;

	    for (r=0; r<N; r++) {
		double	u=(double)r-cx,sum=0.0;

		for (l=-L; l<=L; l++) {
		    /* (h,v) = centre + u*(c,-s) + l*(s,c) ; satisfies
		       (h-cx)*c-(v-cx)*s = u for every l (ray direction) */
		    double	h=cx+u*c+(double)l*s,
				v=cx-u*s+(double)l*c,
				fh,fv;
		    int		h0,v0;

		    if (h<0.0 || v<0.0) continue;
		    h0=(int)h; v0=(int)v;
		    if (h0>=N-1 || v0>=N-1) continue;
		    fh=h-(double)h0; fv=v-(double)v0;
		    sum+=(1.0-fv)*((1.0-fh)*img[(size_t)v0*N+h0]
				  +     fh *img[(size_t)v0*N+h0+1])
			+     fv *((1.0-fh)*img[(size_t)(v0+1)*N+h0]
				  +     fh *img[(size_t)(v0+1)*N+h0+1]);
		}
		sino[(size_t)m*N+r]=(float)sum;
	    }
	}
}
