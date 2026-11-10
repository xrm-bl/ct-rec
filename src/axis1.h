#ifndef CT_AXIS1_H
#define CT_AXIS1_H
/*
 * axis1.h - rotation axis of a 360-degree offset scan from one layer.
 *
 * A one-dimensional version of ofct_DO's search, shared by his_rec / his_sg
 * (his_main.c, rows read from a.his) and ct_gmem (rows read through rhp.h):
 * for every view pair (m, m+M, 180 degrees apart) the -log transmittances of
 * the layer are smoothed with the same Gaussian (OFCT_DO_SMOOTH, default
 * sigma 10), the second one is mirrored, and the mean squared difference over
 * the overlap is accumulated for every shift Dx in [0.5*Nx, Nx-2] (ofct_DO's
 * range, Dy = 0).  The minimum is refined with a parabola and
 * Rc = (Nx + Dx)/2, as in ofct_DO.
 *
 * The caller supplies the rows: row(ctx, view, out) must fill out[0..Nx-1]
 * with the -log transmittance (0 for a non-positive value, the readers'
 * ERROR_VALUE) of sample view "view" (0 .. 2M-1 in angle order, so that view
 * m and view m+M are 180 degrees apart) at the layer.  Only one row per view
 * is needed, so the estimate costs a few seconds even for thousands of views.
 *
 * Requires <stdio.h>, <stdlib.h>, <string.h> and <math.h>.
 */

typedef void	(*Axis1RowFn)(void *ctx,int view,double *out);

static void	Axis1Fail(const char *msg)
{
	(void)fprintf(stderr,"%s\n",msg); exit(1);
}

/* details of the last estimate, for a caller that prints its own summary
   (ct_gmem); Axis1Quiet=1 suppresses the summary line below (the edge
   warning is always printed) */
static int	Axis1Quiet=0,Axis1Pairs=0,Axis1Shift=0,Axis1Rc=0,Axis1Edge=0;
static double	Axis1Dx=0.0;

/* returns the (fractional) axis position; prints the summary on stderr */
static double	EstimateOffsetAxis1(int Nx,int M,long layer,Axis1RowFn row,void *ctx)
{
	int	m,x,k,r=0,Ox1,Ox2,Ox,X;
	double	*S,*G,*H,*tmp,*ker=NULL,sig,best,result;
	char	*env;

	if (M<1) Axis1Fail("offset axis estimate: fewer than two sample views.");
	if ((G=(double *)malloc(sizeof(double)*(size_t)Nx))==NULL ||
	    (H=(double *)malloc(sizeof(double)*(size_t)Nx))==NULL ||
	    (tmp=(double *)malloc(sizeof(double)*(size_t)Nx))==NULL) Axis1Fail("no memory.");

	/* Gaussian kernel of ofct_DO (OFCT_DO_SMOOTH, default 10; 0 = off) */
	sig=((env=getenv("OFCT_DO_SMOOTH"))!=NULL)?atof(env):10.0;
	if (sig>0.0) {
	    double	sm=0.0;
	    r=(int)ceil(3.0*sig);
	    if ((ker=(double *)malloc(sizeof(double)*(size_t)(2*r+1)))==NULL) Axis1Fail("no memory.");
	    for (k=-r; k<=r; k++) { ker[k+r]=exp(-0.5*(double)k*k/(sig*sig)); sm+=ker[k+r]; }
	    for (k=0; k<2*r+1; k++) ker[k]/=sm;
	}

	Ox1=(int)ceil(0.5*Nx); Ox2=Nx-2; Ox=Ox2-Ox1+1;
	if ((S=(double *)calloc((size_t)Ox,sizeof(double)))==NULL) Axis1Fail("no memory.");

	for (m=0; m<M; m++) {
	    int		p;
	    double	*rw[2]={G,H};

	    for (p=0; p<2; p++) {
		row(ctx,(p==0)?m:m+M,tmp);
		if (p==0) memcpy(G,tmp,sizeof(double)*(size_t)Nx);
		else for (x=0; x<Nx; x++) H[Nx-1-x]=tmp[x];	/* mirrored */
		if (sig>0.0) {	/* same truncated, border-renormalised kernel as ofct_DO */
		    double	*F=rw[p],sum,w;
		    for (x=0; x<Nx; x++) {
			sum=w=0.0;
			for (k=-r; k<=r; k++) if (x+k>=0 && x+k<Nx) { sum+=F[x+k]*ker[k+r]; w+=ker[k+r]; }
			tmp[x]=sum/w;
		    }
		    memcpy(F,tmp,sizeof(double)*(size_t)Nx);
		}
	    }
	    for (k=0; k<Ox; k++) {		/* G(x) vs mirrored H(x-Dx) over the overlap */
		/* keep the smoothing radius r away from both ends of the overlap:
		   there the truncated, re-normalised kernel distorts the profiles
		   asymmetrically and biases the minimum (about +0.6 px with the
		   default sigma 10 on test data; exact without smoothing) */
		int	Dx=Ox1+k,x1=Dx+r,x2=Nx-r;
		double	d,ssd=0.0;
		if (x2-x1<8) { x1=Dx; x2=Nx; }
		for (x=x1; x<x2; x++) { d=G[x]-H[x-Dx]; ssd+=d*d; }
		S[k]+=ssd/(double)(x2-x1);
	    }
	}

	best=S[X=0];
	for (k=1; k<Ox; k++) if (S[k]<best) { best=S[k]; X=k; }
	{
	    double	dx=0.0,a,b,c,den,centerf;
	    int		Rc;
	    if (X>=1 && X<=Ox-2) {
		a=S[X-1]; b=S[X]; c=S[X+1]; den=a-2.0*b+c;
		if (den>0.0) { dx=0.5*(a-c)/den; if (dx<-1.0) dx=-1.0; if (dx>1.0) dx=1.0; }
	    }
	    centerf=((double)(Nx+Ox1+X)+dx)/2.0;
	    Rc=(int)floor(centerf+0.5);
	    Axis1Pairs=M; Axis1Shift=Ox1+X; Axis1Dx=dx; Axis1Rc=Rc; Axis1Edge=(X==0 || X==Ox-1);
	    if (!Axis1Quiet)
		(void)fprintf(stderr,"offset axis estimate at layer %ld: %d view pairs, shift %d%+.3f -> centre %.3f (Rc %d)\n",
			      layer,M,Ox1+X,dx,centerf,Rc);
	    result=centerf;
	    if (X==0 || X==Ox-1)
		(void)fprintf(stderr,"warning: the minimum lies at the edge of the search range (centre %d .. %d, as ofct_DO); "
			      "the estimate is unreliable - give center explicitly.\n",(Nx+Ox1)/2,(Nx+Ox2)/2);
	}

	free(S); free(ker); free(tmp); free(H); free(G);
	return result;
}

#endif /* CT_AXIS1_H */
