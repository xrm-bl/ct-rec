
#include <stdlib.h>
#include <math.h>
#include "msd.h"

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

#ifndef	M_PI
#define	M_PI	3.14159265358979323846
#endif

static Complex	*Table(int L)
{
	Complex	*E,*e;
	int	l;
	double	t;

	if ((E=(Complex *)malloc(sizeof(Complex)*L*2))==NULL) return NULL;

	for (e=E, l=0; l<L; l++, e++) {
	    t=M_PI*(double)l/(double)L; e[L].r=(-(e->r=cos(t)));
					e[L].i=(-(e->i=sin(t)));
	}
	return E;
}

int	InitMSD(MSD *msd,int Ngx,int Ngy,int Nhx,int Nhy,
			 int Dx1,int Dy1,int Dx2,int Dy2)
{
	int	Lx,Ly,y;

	msd->Lx=Lx=Range(&Ngx,&Nhx,&Dx1,&Dx2,&(msd->Ogx),&(msd->Ohx));
	msd->Ly=Ly=Range(&Ngy,&Nhy,&Dy1,&Dy2,&(msd->Ogy),&(msd->Ohy));

	if ((msd->Ex=Table(Lx))==NULL ||
	    (msd->Ey=Table(Ly))==NULL ||
	    (msd->Q=(Complex **)malloc(sizeof(Complex *)*Ly))==NULL ||
	    (*(msd->Q)=(Complex *)malloc(sizeof(Complex)*Ly*Lx))==NULL ||
	    (msd->U=(Complex *)malloc(sizeof(Complex)*Ly))==NULL ||
	    (msd->V=(Complex *)malloc(sizeof(Complex)*Ly))==NULL ||
	    (msd->D=(FOM *)malloc(sizeof(double)*(Dy2-Dy1+1)
						*(Dx2-Dx1+1)))==NULL ||
	    (msd->C=(FOM *)malloc(sizeof(FOM)*(Dy2-Dy1+1)
					     *(Dx2-Dx1+1)))==NULL) return 0;

	for (y=1; y<Ly; y++) msd->Q[y]=msd->Q[y-1]+Lx;

	msd->Ngx=Ngx; msd->Ngy=Ngy; msd->Nhx=Nhx; msd->Nhy=Nhy;
	msd->Dx1=Dx1; msd->Dy1=Dy1; msd->Dx2=Dx2; msd->Dy2=Dy2; return 1;
}

static void	FFT(int S,int L,Complex *E,Complex *C)
{
	int	j,k,m,n,p,q,L2=L>>1;
	Complex	c,*e,*a,*b;
	double	r,i,d,s=(double)S;

	for (j=L-1, k=L-2; k>0; k--) {
	    for (m=L2; (j^=m)&m; m>>=1) ;

	    if (j>k) {
		c=C[j]; C[j]=C[k]; C[k]=c;
	    }
	}
	for (j=2, k=1, m=L2, n=L; n>1; n=m, m>>=1, k=j, j<<=1)
	    for (e=E, p=0; p<k; p++, e+=n) {
		r=e->r; i=e->i*s;
		for (a=C+p, q=0; q<m; q++, a+=j) {
		    c=(*(b=a+k)); b->r=a->r-(d=c.r*r-c.i*i); a->r+=d;
				  b->i=a->i-(d=c.r*i+c.i*r); a->i+=d;
		}
	    }
}

static Complex	Product(double r1,double i1,double r2,double i2)
{
	Complex	c;

	c.r=r1*r2-i1*i2; c.i=r1*i2+i1*r2; return c;
}

void	CalcMSD(MSD *msd,FOM **G,FOM **H)
{
	FOM	**Gy,**Hy,*Gx,*Hx,*D,*C;
	int	y,x,gx,gy,hx,hy,Lxx,z,
		Ngx=msd->Ngx,Ngy=msd->Ngy,Nhx=msd->Nhx,Nhy=msd->Nhy,
		Dx1=msd->Dx1,Dy1=msd->Dy1,Dx2=msd->Dx2,Dy2=msd->Dy2,
		Ogx=msd->Ogx,Ogy=msd->Ogy,Ohx=msd->Ohx,Ohy=msd->Ohy,
		Lx=msd->Lx,Ly=msd->Ly,
		Ngx1=Ngx-1,Ngy1=Ngy-1,Nhx1=Nhx-1,Nhy1=Nhy-1,
		Lx1=Lx-1,Ly1=Ly-1,Lx2=Lx*2,Ly2=Ly*2;
	double	d,gr,gi,hr,hi,t,r,i,*P,l,ly,lyx,
		ngx=(double)Ngx,ngy=(double)Ngy,
		nhx=(double)Nhx,nhy=(double)Nhy,
		ngx2=(double)(Ngx%2),ngy2=(double)(Ngy%2),
		nhx2=(double)(Nhx%2),nhy2=(double)(Nhy%2);
	Complex	*ex,*ey,*q,*u,*v,g,h,c,
		*Ex=msd->Ex,*Ey=msd->Ey,**Q=msd->Q,*U=msd->U,*V=msd->V;

	for (Gy=G+Ogy, Hy=H+Ohy, y=0; y<Ly; y++) {
	    x=0; q=Q[y];
	    if (y<Ngy)
		for (Gx=(*(Gy++))+Ogx; x<Ngx; x++) {
		    d=(*(Gx++)); (q++)->r=d*d;
		}

	    for (; x<Lx; x++) (q++)->r=0.0;

	    x=0; q=Q[y];
	    if (y<Nhy)
		for (Hx=(*(Hy++))+Ohx; x<Nhx; x++) {
		    d=(*(Hx++)); (q++)->i=d*d;
		}

	    for (; x<Lx; x++) (q++)->i=0.0;

	    if (y<Ngy || y<Nhy) FFT(-1,Lx,Ex,Q[y]);
	}
/*
    x=0
*/
	for (u=U, y=0; y<Ly; y++) *(u++)=Q[y][0]; FFT(-1,Ly,Ey,U);

	U->r=U->r*nhx*nhy+ngx*ngy*U->i; U->i=0.0;
	for (ey=Ey+1, gy=Ngy1, hy=Nhy1, u=U+1, v=U+Ly1;
	     u<v;
	     ey++, gy=(gy+Ngy1)%Ly2, hy=(hy+Nhy1)%Ly2, u++, v--) {
	    t=ey->r/ey->i;
	    r=Ey[hy].r; i=Ey[hy].i; d=nhx*(r+i*t);
	    g=Product(u->r+v->r,u->i-v->i,d*r, d*i);
	    r=Ey[gy].r; i=Ey[gy].i; d=ngx*(r+i*t);
	    h=Product(d*r,-d*i,u->i+v->i,u->r-v->r);
	    v->r=u->r=(g.r+h.r)/2.0; v->i=(-(u->i=(g.i+h.i)/2.0));
	}
	u->r=u->r*nhx*nhy2+ngx*ngy2*u->i; u->i=0.0;

	FFT(1,Ly,Ey,U); for (u=U, y=0; y<Ly; y++) Q[y][0]=(*(u++));
/*
    x=1~Lx/2-1
*/
	for (ex=Ex+1, gx=Ngx1, hx=Nhx1, x=1, Lxx=Lx1;
	     x<Lxx;
	     ex++, gx=(gx+Ngx1)%Lx2, hx=(hx+Nhx1)%Lx2, x++, Lxx--) {
	    for (u=U, v=V, y=0; y<Ly; y++) {
		*(u++)=Q[y][x]; *(v++)=Q[y][Lxx];
	    }
	    FFT(-1,Ly,Ey,U); FFT(-1,Ly,Ey,V);

	    t=ex->r/ex->i;
	    r=Ex[gx].r; i=Ex[gx].i; d=r+i*t; gr=d*r; gi=d*i;
	    r=Ex[hx].r; i=Ex[hx].i; d=r+i*t; hr=d*r; hi=d*i;
	    g=Product(U->r+V->r,U->i-V->i,hr*nhy, hi*nhy);
	    h=Product(gr*ngy,-gi*ngy,U->i+V->i,U->r-V->r);
	    V->r=U->r=(g.r+h.r)/2.0; V->i=(-(U->i=(g.i+h.i)/2.0));
	    for (ey=Ey+1, gy=Ngy1, hy=Nhy1, u=U+1, v=V+Ly1, y=1;
		 y<Ly;
		 ey++, gy=(gy+Ngy1)%Ly2, hy=(hy+Nhy1)%Ly2, u++, v--, y++) {
		t=ey->r/ey->i;
		r=Ey[hy].r; i=Ey[hy].i; d=r+i*t; c=Product(hr,hi,d*r,d*i);
		g=Product(u->r+v->r,u->i-v->i,c.r, c.i);
		r=Ey[gy].r; i=Ey[gy].i; d=r+i*t; c=Product(gr,gi,d*r,d*i);
		h=Product(c.r,-c.i,u->i+v->i,u->r-v->r);
		v->r=u->r=(g.r+h.r)/2.0; v->i=(-(u->i=(g.i+h.i)/2.0));
	    }
	    FFT(1,Ly,Ey,U); FFT(1,Ly,Ey,V);
	    for (u=U, v=V, y=0; y<Ly; y++) {
		Q[y][x]=(*(u++)); Q[y][Lxx]=(*(v++));
	    }
	}
/*
    x=Lx/2
*/
	for (u=U, y=0; y<Ly; y++) *(u++)=Q[y][x]; FFT(-1,Ly,Ey,U);

	U->r=U->r*nhx2*nhy+ngx2*ngy*U->i; U->i=0.0;
	for (ey=Ey+1, gy=Ngy1, hy=Nhy1, u=U+1, v=U+Ly1;
	     u<v;
	     ey++, gy=(gy+Ngy1)%Ly2, hy=(hy+Nhy1)%Ly2, u++, v--) {
	    t=ey->r/ey->i;
	    r=Ey[hy].r; i=Ey[hy].i; d=nhx2*(r+i*t);
	    g=Product(u->r+v->r,u->i-v->i,d*r, d*i);
	    r=Ey[gy].r; i=Ey[gy].i; d=ngx2*(r+i*t);
	    h=Product(d*r,-d*i,u->i+v->i,u->r-v->r);
	    v->r=u->r=(g.r+h.r)/2.0; v->i=(-(u->i=(g.i+h.i)/2.0));
	}
	u->r=u->r*nhx2*nhy2+ngx2*ngy2*u->i; u->i=0.0;

	FFT(1,Ly,Ey,U); for (u=U, y=0; y<Ly; y++) Q[y][x]=(*(u++));
/*
    ----
*/
	P=(double *)msd->D;
	for (y=Dy1; y<=Dy2; y++) {
	    FFT(1,Lx,Ex,q=Q[(y<0)?Ly+y:y]);

	    for (x=Dx1; x<=Dx2; x++) *(P++)=q[(x<0)?Lx+x:x].r;
	}
/*
    ----
*/
	for (Gy=G+Ogy, Hy=H+Ohy, y=0; y<Ly; y++) {
	    x=0; q=Q[y];
	    if (y<Ngy) for (Gx=(*(Gy++))+Ogx; x<Ngx; x++) (q++)->r=(*(Gx++));

	    for (; x<Lx; x++) (q++)->r=0.0;

	    x=0; q=Q[y];
	    if (y<Nhy) for (Hx=(*(Hy++))+Ohx; x<Nhx; x++) (q++)->i=(*(Hx++));

	    for (; x<Lx; x++) (q++)->i=0.0;

	    if (y<Ngy || y<Nhy) FFT(-1,Lx,Ex,Q[y]);
	}
/*
    x=0
*/
	for (u=U, y=0; y<Ly; y++) *(u++)=Q[y][0]; FFT(-1,Ly,Ey,U);

	U->r*=U->i; U->i=0.0;
	for (u=U+1, v=U+Ly1; u<v; u++, v--) {
	    *u=Product((u->r+v->r)/2.0,(u->i-v->i)/2.0,
		       (u->i+v->i)/2.0,(u->r-v->r)/2.0);
	    v->r=u->r; v->i=(-u->i);
	}
	u->r*=u->i; u->i=0.0;

	FFT(1,Ly,Ey,U); for (u=U, y=0; y<Ly; y++) Q[y][0]=(*(u++));
/*
    x=1~Lx/2-1
*/
	for (x=1, Lxx=Lx1; x<Lxx; x++, Lxx--) {
	    for (u=U, v=V, y=0; y<Ly; y++) {
		*(u++)=Q[y][x]; *(v++)=Q[y][Lxx];
	    }
	    FFT(-1,Ly,Ey,U); FFT(-1,Ly,Ey,V);

	    *U=Product((U->r+V->r)/2.0,(U->i-V->i)/2.0,
		       (U->i+V->i)/2.0,(U->r-V->r)/2.0);
	    V->r=U->r; V->i=(-U->i);
	    for (u=U+1, v=V+Ly1, y=1; y<Ly; y++, u++, v--) {
		*u=Product((u->r+v->r)/2.0,(u->i-v->i)/2.0,
			   (u->i+v->i)/2.0,(u->r-v->r)/2.0);
		v->r=u->r; v->i=(-u->i);
	    }
	    FFT(1,Ly,Ey,U); FFT(1,Ly,Ey,V);
	    for (u=U, v=V, y=0; y<Ly; y++) {
		Q[y][x]=(*(u++)); Q[y][Lxx]=(*(v++));
	    }
	}
/*
    x=Lx/2
*/
	for (u=U, y=0; y<Ly; y++) *(u++)=Q[y][x]; FFT(-1,Ly,Ey,U);

	U->r*=U->i; U->i=0.0;
	for (u=U+1, v=U+Ly1; u<v; u++, v--) {
	    *u=Product((u->r+v->r)/2.0,(u->i-v->i)/2.0,
		       (u->i+v->i)/2.0,(u->r-v->r)/2.0);
	    v->r=u->r; v->i=(-u->i);
	}
	u->r*=u->i; u->i=0.0;

	FFT(1,Ly,Ey,U); for (u=U, y=0; y<Ly; y++) Q[y][x]=(*(u++));
/*
    ----
*/
	l=(double)Lx*(double)Ly; P=(double *)(D=msd->D); C=msd->C;
	for (y=Dy1; y<=Dy2; y++) {
	    if (y<0) {
		q=Q[Ly+y]; ly=l*(double)(((z=Nhy+y)<Ngy)?z:Ngy);
	    }
	    else {
		q=Q[   y]; ly=l*(double)(((z=Ngy-y)<Nhy)?z:Nhy);
	    }
	    FFT(1,Lx,Ex,q);

	    for (x=Dx1; x<=Dx2; x++) {
		if (x<0)
		    *(C++)=(d=q[Lx+x].r)/
			   (lyx=ly*(double)(((z=Nhx+x)<Ngx)?z:Ngx));
		else
		    *(C++)=(d=q[   x].r)/
			   (lyx=ly*(double)(((z=Ngx-x)<Nhx)?z:Nhx));

		*(D++)=(*(P++)-2.0*d)/lyx;
	    }
	}
}

void	TermMSD(MSD *msd)
{
	free(msd->C); free(msd->D);
	free(msd->V); free(msd->U);
	free(*(msd->Q)); free(msd->Q);
	free(msd->Ey); free(msd->Ex);
}
