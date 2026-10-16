/*
 * rhis_hp.c - the rhp.h reader API (InitReadHiPic / ReadHiPic /
 *             ReadHiPicBand / TermReadHiPic) on top of the HisSet layer of
 *             rhis.c, so that hp_tg_ku.c, ofct_srec.c and ofct_DO.c link
 *             against a.his + conv.bat + output.log instead of rhp_c.c.
 *
 * The arithmetic is copied from rhp_c.c, so the transmittances are
 * bit-identical.  The HiPic struct and ReadHiPicBand depend on FOM, so a
 * program that needs both FOM=float (hp_tg) and FOM=double (ofct_srec)
 * compiles this file twice, once with -DFOM=float -DRHP_NAME_SUFFIX=_f
 * (see rhp.h).  Compile with -DONLY_CT_VIEWS like the callers.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "rhis.h"

#define LEN		2048
#define WORD		unsigned short

#include "rhp.h"
#include "blacklim.h"

#define EPS	1e-9	/* for time interval */

static HisSet	*rhp_his=NULL;	/* one data set per process (as rhp_c.c) */

static void	Error(char *dir,char *name,char *msg)
{
	if (*dir!='\0') (void)fprintf(stderr,"%s/",dir);

	(void)fprintf(stderr,"%s : %s\n",name,msg); exit(1);
}

static int	Compare(OutputLog *OL1,OutputLog *OL2)
{
	int	it=OL2->it-OL1->it;

	return (it)?it:(OL1->it==0)?(OL1->c>OL2->c)?1:(OL1->c<OL2->c)?-1:0:
		       (OL1->it==1)?(OL1->a>OL2->a)?1:(OL1->a<OL2->a)?-1:0:0;
}

static size_t	safe_count(size_t a,size_t b,size_t c)
{
	if (b!=0 && a>SIZE_MAX/b) goto over;
	a*=b;
	if (c!=0 && a>SIZE_MAX/c) goto over;
	return a*c;
over:
	(void)fprintf(stderr,"size overflow in image allocation.\n"); exit(1);
}

#define MALLOC(type,noe)	(type *)malloc(sizeof(type)*(noe))

void	InitReadHiPic(char *dir,HiPic *hp)
{
	char		*env,str[LEN],output_log[LEN]="",
			*bsn="bad sequence number.",
			*dsn="duplicated sequence number.";
	int		l,q,it,i;
	FILE		*file;
	double		c,a,c1=0.0,c2=0.0;
	OutputLog	*OL;

	if ((env=getenv("RHP_O"))!=NULL) {
	    if (strlen(env)>=LEN) Error("","RHP_O","environment value too long.");
	    (void)snprintf(output_log,LEN,"%s",env);
	}
	if (*output_log=='\0')
	    if (snprintf(output_log,LEN,"%s/output.log",dir)>=LEN)
		Error(dir,"output.log","path too long.");

	/* warm the cached values up while still single-threaded:
	   ReadHiPicBand runs on several reader threads. */
	(void)BlackThresh();
	(void)BlackFrac();

	if (rhp_his!=NULL) Error("",dir,"InitReadHiPic called twice.");
	rhp_his=HisOpen(dir);

	hp->Nq=HisNq(rhp_his);
	hp->Nx=HisNx(rhp_his);
	hp->Ny=HisNy(rhp_his);
	hp->q_img=NULL;			/* no files: frames live in the his file */

	if ((hp->dir=strdup(dir))==NULL) Error("",dir,"no memory for directory data.");

	if ((hp->OL=MALLOC(OutputLog,hp->Nq+1))==NULL)
	    Error("",output_log,"no memory for log data.");

	for (q=0; q<=hp->Nq; q++) hp->OL[q].it=(-1);

	hp->Ni=hp->Nt=0;

	if (*output_log=='-')
	    file=stdin;
	else
	    if ((file=fopen(output_log,"r"))==NULL)
		Error("",output_log,"file not open.");

	while (fgets(str,LEN,file)!=NULL)
	    if ((sscanf(str,"%d %lf %lf %d%n",&q,&c,&a,&it,&l)==4)) {
		str[l]='\0';

		if (q<0 || q>=hp->Nq || !HisHasQ(rhp_his,q)) Error("",str,bsn);

		if (it!=0 && it!=1) Error("",str,"unacceptable log data.");

		if (hp->OL[q].it!=(-1)) Error("",str,dsn);

		if (it==0) ++(hp->Ni); else ++(hp->Nt);

		if (hp->Ni+hp->Nt==1) c1=c2=c; else if (c1>c) c1=c;
					       else if (c2<c) c2=c;

		hp->OL[q].q=q; hp->OL[q].c=c; hp->OL[q].a=a; hp->OL[q].it=it;
	    }
	    else
#ifndef	ONLY_CT_VIEWS
		if (file==stdin)
#endif
		break;

	if (file!=stdin) (void)fclose(file);

	if (hp->Ni==0) Error(dir,output_log,"no view referred as I0-image.");
	if (hp->Nt==0) Error(dir,output_log,"no view referred as I-image.");

	if ((hp->D=MALLOC(WORD,safe_count((size_t)hp->Ny,(size_t)hp->Nx,1)))==NULL ||
	    (hp->I=MALLOC(WORD *,hp->Ni+(hp->Ni==1)))==NULL ||
	    (*(hp->I)=MALLOC(WORD,safe_count((size_t)hp->Ni,(size_t)hp->Ny,(size_t)hp->Nx)))==NULL ||
	    (hp->T=MALLOC(FOM *,hp->Ny))==NULL ||
	    (*(hp->T)=MALLOC(FOM,safe_count((size_t)hp->Ny,(size_t)hp->Nx,1)))==NULL)
	    Error("",dir,"no memory for image data.");

	HisReadDarkBand(rhp_his,0,hp->Ny-1,hp->D);

	qsort(hp->OL,hp->Nq+1,sizeof(OutputLog),(int (*)(const void *,const void *))Compare);

	OL=hp->OL+hp->Nt;
	for (i=0; i<hp->Ni; i++) {
	    if (i>0) hp->I[i]=hp->I[i-1]+(size_t)hp->Ny*hp->Nx;
	    HisReadQBand(rhp_his,OL[i].q,0,hp->Ny-1,hp->I[i]);
	}
	if (hp->Ni==1) {
	    hp->Ni=2; OL[1].q=OL[0].q;
		      OL[0].c=c1-EPS; OL[1].c=c2+EPS;
		      OL[1].a=OL[0].a; OL[1].it=OL[0].it; hp->I[1]=hp->I[0];
	}
	for (i=1; i<hp->Ny; i++) hp->T[i]=hp->T[i-1]+hp->Nx;
}

/* ReadHiPicBand: same arithmetic as rhp_c.c, frame comes from the his file.
   raw receives the whole frame (Ny*Nx) because the black-frame decision
   is made over the whole frame, exactly as in rhp_c.c. */
void	ReadHiPicBand(HiPic *hp,int t,int y1,int y2,FOM **dst,unsigned short *raw)
{
	char		str[LEN];
	OutputLog	*ol,*OL;
	int		i,y,x;
	WORD		*D,*I1,*I2,*T;
	double		r1,r2,I_D,T_D;
	double		black_thresh,black_frac;
	int		nclip;
	size_t		o;

	if (t<0 || t>=hp->Nt) {
	    (void)sprintf(str,"%d",t); Error("",str,"bad sequence number.");
	}
	if (y1<0) y1=0;
	if (y2>hp->Ny-1) y2=hp->Ny-1;

	black_thresh=BlackThresh();
	black_frac=BlackFrac();

	ol=hp->OL+t;
	HisReadQBand(rhp_his,ol->q,0,hp->Ny-1,(WORD *)raw);

	OL=hp->OL+hp->Nt;
	for (i=hp->Ni-2; i>0; i--) if (ol->c>=OL[i].c) break;

	r1=1.0-(r2=(ol->c-OL[i].c)/(OL[i+1].c-OL[i].c));

	nclip=0;
	D=hp->D; T=(WORD *)raw;
	for (y=0; y<hp->Ny; y++)
	for (x=0; x<hp->Nx; x++)
	    if ((double)*T++ -(double)*D++ < black_thresh) ++nclip;

	if ((double)nclip > black_frac*(double)hp->Nx*(double)hp->Ny) {
	    (void)fprintf(stderr,
		"Warning\t black\t t=%d clipped=%d/%d (%.1f%%, thresh=%.4g)\n",
		t,nclip,hp->Nx*hp->Ny,
		100.0*(double)nclip/((double)hp->Nx*(double)hp->Ny),black_thresh);
	    BlackCountProjection();
	    for (y=y1; y<=y2; y++)
	    for (x=0; x<hp->Nx; x++)
		dst[y-y1][x]=ERROR_VALUE;

	    return;
	}
	for (y=y1; y<=y2; y++) {
	    o=(size_t)y*(size_t)hp->Nx;
	    D=hp->D+o; I1=hp->I[i]+o; I2=hp->I[i+1]+o; T=(WORD *)raw+o;
	    for (x=0; x<hp->Nx; x++) {
		T_D=(double)T[x]-(double)D[x];
		I_D=r1*(double)I1[x]+r2*(double)I2[x]-(double)D[x];
		if (!(T_D >= black_thresh)) T_D=black_thresh;
		if (!(I_D >= black_thresh)) I_D=black_thresh;
		dst[y-y1][x]=T_D/I_D;
	    }
	}
}

/* ReadHiPic: whole frame into hp->T (same result as rhp_c.c's ReadHiPic) */
void	ReadHiPic(HiPic *hp,int t)
{
	WORD	*raw=(WORD *)malloc(sizeof(WORD)*(size_t)hp->Ny*(size_t)hp->Nx);

	if (raw==NULL) Error("",hp->dir,"no memory for a frame.");
	ReadHiPicBand(hp,t,0,hp->Ny-1,hp->T,raw);
	free(raw);
}

void	TermReadHiPic(HiPic *hp)
{
	free(*(hp->T)); free(hp->T);
	free(*(hp->I)); free(hp->I); free(hp->D);

	free(hp->OL);
	free(hp->dir);

	HisClose(rhp_his); rhp_his=NULL;
}
