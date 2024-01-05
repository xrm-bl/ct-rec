
#include <stdlib.h>
#include <math.h>
#include "fft.h"

#ifndef	M_PI
#define	M_PI	3.14159265358979323846
#endif

void	FFT(int S,int L,Complex *E,Complex *C)
{
	int	j,k,m,n,p,q,L2=L>>1;
	Complex	c,*e,*a,*b;
	double	r,i,d;

	for (j=L-1, k=L-2; k>0; k--) {
	    for (m=L2; (j^=m)&m; m>>=1) ;

	    if (j>k) {
		c=C[j]; C[j]=C[k]; C[k]=c;
	    }
	}
	for (j=2, k=1, m=L2, n=L; n>1; n=m, m>>=1, k=j, j<<=1)
	    for (e=E, p=0; p<k; p++, e+=n) {
		r=e->r; i=e->i*(double)S;
		for (a=C+p, q=0; q<m; q++, a+=j) {
		    c=(*(b=a+k)); b->r=a->r-(d=c.r*r-c.i*i); a->r+=d;
				  b->i=a->i-(d=c.r*i+c.i*r); a->i+=d;
		}
	    }
}

int	SetUpFFT(int N,Complex **E)
{
	int	L,l;
	Complex	*e;
	double	d;

	for (L=1; L<N; L<<=1) ;

	if ((*E=e=(Complex *)malloc(sizeof(Complex)*L))==NULL) return 0;

	for (l=0; l<L; l++, e++) {
	    d=M_PI*(double)l/(double)L; e->r=cos(d); e->i=sin(d);
	}
	return L;
}
