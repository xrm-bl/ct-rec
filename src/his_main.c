/*
 * his_main.c - entry point of his_rec / his_srec / his_DO.
 *
 * All three work on the data set as it comes from the beamline: a.his,
 * conv.bat and output.log in one directory (read by rhis.c / rhis_hp.c).
 * The third column of output.log is the rotation angle; when the sample
 * views reach HIS_OFFSET_ANGLE degrees (350, the test ofct_rec has always
 * applied) the scan is a 360-degree offset scan, otherwise a 180-degree
 * scan.  The reconstruction bodies are the unchanged programs, linked
 * side by side with their main() renamed (-DHIS_INPUT):
 *
 *   program    180-degree scan        360-degree offset scan
 *   his_rec    ct_rec_main            ofct_rec_main       (run beside a.his)
 *   his_srec   hp_tg_main             ofct_srec_main      (run one level up)
 *   his_DO     hp2DO_compute          ofct_DO_main        (run one level up)
 *
 * Command lines:
 *   his_rec  layer {center} {pixel size} {offsetangle}      (as ct_rec;
 *            without center the axis is estimated from this layer: a
 *            180-degree scan by ct_rec's own search, an offset scan by a
 *            one-dimensional version of ofct_DO's search, see
 *            EstimateOffsetAxis below)
 *   his_srec HiPic Dr RC RA0 rec                            (as hp_tg)
 *   his_srec HiPic Dr L1 C1 L2 C2 RA0 rec
 *            An offset scan is handed to ofct_srec as
 *            "HiPic RC Oy rangeList Dr RA0 rec" with Oy = 0 (or the
 *            environment variable HIS_SREC_OY) and rangeList "-" or
 *            "L1-L2"; a tilted axis (C1 != C2) is not supported there.
 *   his_DO   HiPic {stride}            (CPU)   his_DO_g HiPic {navg {stride}}
 *            180-degree scan: the hp2DO estimate on the six views named
 *            on the hp2do line of conv.bat (dark, I0 and view at 0 deg,
 *            dark, I0 and view at 180 deg); without such a line the
 *            ct_rec convention q1 q2 q(last+1) q(last+2) is used.  The
 *            extra arguments are ignored in this mode.
 *            360-degree scan: ofct_DO / ofct_DO_g.
 *   Both his_srec and his_DO print a ready-to-run suggestion to stderr.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "rhis.h"

#define LEN			2048
#define HIS_OFFSET_ANGLE	350.0

static void	Fail(const char *msg)
{
	(void)fprintf(stderr,"%s\n",msg); exit(1);
}

/* read <dir>/output.log (or RHP_O) and report the largest angle of a
   sample view and the last q number before the first blank line */
static double	ScanAngles(const char *dir,int *q_last)
{
	char	path[LEN],str[LEN],*env;
	FILE	*f;
	int	q,it,n=0,blank_seen=0,ql=-1;
	double	c,a,amax=-1.0e300;

	if ((env=getenv("RHP_O"))!=NULL && *env!='\0')
	    (void)snprintf(path,LEN,"%s",env);
	else if (snprintf(path,LEN,"%s/output.log",dir)>=LEN)
	    Fail("output.log : path too long.");
	if (!strcmp(path,"-")) Fail("his_*: output.log on stdin (RHP_O=-) is not supported.");
	if ((f=fopen(path,"r"))==NULL) {
	    (void)fprintf(stderr,"%s : file not open.\n",path); exit(1);
	}
	while (fgets(str,LEN,f)!=NULL) {
	    if (sscanf(str,"%d %lf %lf %d",&q,&c,&a,&it)==4) {
		if (it==1 && a>amax) amax=a;
		if (!blank_seen) ql=q;
		n++;
	    }
	    else blank_seen=1;
	}
	(void)fclose(f);
	if (n==0) { (void)fprintf(stderr,"%s : no view found.\n",path); exit(1); }
	if (q_last!=NULL) *q_last=ql;
	return amax;
}

static int	IsOffsetScan(const char *dir,const char *prog)
{
	double	amax=ScanAngles(dir,NULL);
	int	off=(amax>=HIS_OFFSET_ANGLE);

	(void)fprintf(stderr,"%s: sample views reach %.1f degrees -> %s\n",
		      prog,amax,off?"360-degree offset scan":"180-degree scan");
	return off;
}

static const char	*BaseName(const char *p)
{
	const char	*s=strrchr(p,'/'),*b=strrchr(p,'\\');

	if (b!=NULL && (s==NULL || b>s)) s=b;
	return (s!=NULL)?s+1:p;
}

/* ====================================================================== */
#if defined(HIS_REC)

extern int	ct_rec_main(int argc,char **argv);
extern int	ofct_rec_main(int argc,char **argv);

#include "blacklim.h"	/* same per-pixel floor as the readers */

/* ---- rotation axis of an offset scan from one layer ----
   A one-dimensional version of ofct_DO's search: for every view pair
   (m, m+M, 180 degrees apart) the -log transmittances of the layer are
   smoothed with the same Gaussian (OFCT_DO_SMOOTH, default sigma 10),
   the second one is mirrored, and the mean squared difference over the
   overlap is accumulated for every shift Dx in [0.5*Nx, Nx-2] (ofct_DO's
   range, Dy = 0).  The minimum is refined with a parabola and
   Rc = (Nx + Dx)/2, as in ofct_DO.  Only the rows of that layer are read
   from a.his, so this costs a few seconds even for thousands of views. */

typedef struct { int q; double c,a; int it; } View;

static int	CmpI0(const void *a,const void *b)
{
	double	d=((const View *)a)->c-((const View *)b)->c;
	return (d>0.0)?1:(d<0.0)?-1:0;
}
static int	CmpView(const void *a,const void *b)
{
	double	d=((const View *)a)->a-((const View *)b)->a;
	return (d>0.0)?1:(d<0.0)?-1:0;
}

/* output.log up to the first line that is not a view (as the ONLY_CT_VIEWS
   readers do); returns the number of entries */
static int	ReadLog(const char *dir,View **out)
{
	char	path[LEN],str[LEN],*env;
	FILE	*f;
	View	*v=NULL;
	int	n=0,cap=0,q,it;
	double	c,a;

	if ((env=getenv("RHP_O"))!=NULL && *env!='\0') (void)snprintf(path,LEN,"%s",env);
	else (void)snprintf(path,LEN,"%s/output.log",dir);
	if ((f=fopen(path,"r"))==NULL) { (void)fprintf(stderr,"%s : file not open.\n",path); exit(1); }
	while (fgets(str,LEN,f)!=NULL) {
	    if (sscanf(str,"%d %lf %lf %d",&q,&c,&a,&it)!=4) break;
	    if (n==cap) {
		cap=cap?2*cap:256;
		if ((v=(View *)realloc(v,sizeof(View)*(size_t)cap))==NULL) Fail("no memory for output.log.");
	    }
	    v[n].q=q; v[n].c=c; v[n].a=a; v[n].it=it; n++;
	}
	(void)fclose(f);
	*out=v; return n;
}

static double	Log0(double d) { return (d>0.0)?log(d):0.0; }

static int	EstimateOffsetAxis(const char *dir,long layer)
{
	HisSet	*s=HisOpen(dir);
	View	*v,*I0,*V;
	int	n,Ni=0,Nt=0,i,m,M,Nx,x,k,r=0,Ox1,Ox2,Ox,X,Rc;
	unsigned short	*D,*I0row,*T1,*T2;
	double	*S,*G,*H,*tmp,*ker=NULL,thr,c1,c2,sig,best;
	char	*env;

	Nx=HisNx(s);
	if (layer<0 || layer>=HisNy(s)) {
	    (void)fprintf(stderr,"layer %ld is outside the frame (0 .. %d)\n",layer,HisNy(s)-1); exit(1);
	}
	n=ReadLog(dir,&v);
	for (i=0; i<n; i++) { if (v[i].it==0) Ni++; else if (v[i].it==1) Nt++; }
	if (Ni==0 || Nt<2) Fail("offset axis estimate: output.log has no I0 view or fewer than two sample views.");
	if ((I0=(View *)malloc(sizeof(View)*(size_t)(Ni+1)))==NULL ||
	    (V=(View *)malloc(sizeof(View)*(size_t)Nt))==NULL) Fail("no memory.");
	Ni=Nt=0;
	for (i=0; i<n; i++) { if (v[i].it==0) I0[Ni++]=v[i]; else if (v[i].it==1) V[Nt++]=v[i]; }
	qsort(I0,(size_t)Ni,sizeof(View),CmpI0);
	qsort(V,(size_t)Nt,sizeof(View),CmpView);
	M=Nt/2;
	if (Nt%2) (void)fputs("offset axis estimate: odd number of views (warning).\n",stderr);

	/* the layer's row of dark and of every I0 */
	if ((D=(unsigned short *)malloc(sizeof(unsigned short)*(size_t)Nx))==NULL ||
	    (I0row=(unsigned short *)malloc(sizeof(unsigned short)*(size_t)Nx*(size_t)(Ni+1)))==NULL ||
	    (T1=(unsigned short *)malloc(sizeof(unsigned short)*(size_t)Nx))==NULL ||
	    (T2=(unsigned short *)malloc(sizeof(unsigned short)*(size_t)Nx))==NULL ||
	    (G=(double *)malloc(sizeof(double)*(size_t)Nx))==NULL ||
	    (H=(double *)malloc(sizeof(double)*(size_t)Nx))==NULL ||
	    (tmp=(double *)malloc(sizeof(double)*(size_t)Nx))==NULL) Fail("no memory.");
	HisReadDarkBand(s,(int)layer,(int)layer,D);
	for (i=0; i<Ni; i++) {
	    if (!HisHasQ(s,I0[i].q)) { (void)fprintf(stderr,"q%04d is not produced by conv.bat\n",I0[i].q); exit(1); }
	    HisReadQBand(s,I0[i].q,(int)layer,(int)layer,I0row+(size_t)i*Nx);
	}
	if (Ni==1) {		/* as the readers: one I0 -> constant in time */
	    memcpy(I0row+(size_t)Nx,I0row,sizeof(unsigned short)*(size_t)Nx);
	    c1=c2=I0[0].c;
	    for (i=0; i<Nt; i++) { if (V[i].c<c1) c1=V[i].c; if (V[i].c>c2) c2=V[i].c; }
	    I0[1]=I0[0]; I0[0].c=c1-1e-9; I0[1].c=c2+1e-9; Ni=2;
	}
	thr=BlackThresh();

	/* Gaussian kernel of ofct_DO (OFCT_DO_SMOOTH, default 10; 0 = off) */
	sig=((env=getenv("OFCT_DO_SMOOTH"))!=NULL)?atof(env):10.0;
	if (sig>0.0) {
	    double	sm=0.0;
	    r=(int)ceil(3.0*sig);
	    if ((ker=(double *)malloc(sizeof(double)*(size_t)(2*r+1)))==NULL) Fail("no memory.");
	    for (k=-r; k<=r; k++) { ker[k+r]=exp(-0.5*(double)k*k/(sig*sig)); sm+=ker[k+r]; }
	    for (k=0; k<2*r+1; k++) ker[k]/=sm;
	}

	Ox1=(int)ceil(0.5*Nx); Ox2=Nx-2; Ox=Ox2-Ox1+1;
	if ((S=(double *)calloc((size_t)Ox,sizeof(double)))==NULL) Fail("no memory.");

	for (m=0; m<M; m++) {
	    int		j,p;
	    double	r1,r2,t_d,i_d;
	    double	*row[2]={G,H};
	    View	*pv[2]={V+m,V+m+M};

	    for (p=0; p<2; p++) {
		unsigned short	*T=(p==0)?T1:T2,*Ia,*Ib;

		if (!HisHasQ(s,pv[p]->q)) { (void)fprintf(stderr,"q%04d is not produced by conv.bat\n",pv[p]->q); exit(1); }
		HisReadQBand(s,pv[p]->q,(int)layer,(int)layer,T);
		for (j=Ni-2; j>0; j--) if (pv[p]->c>=I0[j].c) break;
		r2=(pv[p]->c-I0[j].c)/(I0[j+1].c-I0[j].c); r1=1.0-r2;
		Ia=I0row+(size_t)j*Nx; Ib=I0row+(size_t)(j+1)*Nx;
		for (x=0; x<Nx; x++) {
		    t_d=(double)T[x]-(double)D[x];
		    i_d=r1*(double)Ia[x]+r2*(double)Ib[x]-(double)D[x];
		    if (!(t_d>=thr)) t_d=thr;
		    if (!(i_d>=thr)) i_d=thr;
		    if (p==0) G[x]=-Log0(t_d/i_d); else H[Nx-1-x]=-Log0(t_d/i_d);
		}
		if (sig>0.0) {	/* same truncated, border-renormalised kernel as ofct_DO */
		    double	*F=row[p],sum,w;
		    for (x=0; x<Nx; x++) {
			sum=w=0.0;
			for (k=-r; k<=r; k++) if (x+k>=0 && x+k<Nx) { sum+=F[x+k]*ker[k+r]; w+=ker[k+r]; }
			tmp[x]=sum/w;
		    }
		    memcpy(F,tmp,sizeof(double)*(size_t)Nx);
		}
	    }
	    for (k=0; k<Ox; k++) {		/* G(x) vs mirrored H(x-Dx) over the overlap */
		int	Dx=Ox1+k;
		double	d,ssd=0.0;
		for (x=Dx; x<Nx; x++) { d=G[x]-H[x-Dx]; ssd+=d*d; }
		S[k]+=ssd/(double)(Nx-Dx);
	    }
	}

	best=S[X=0];
	for (k=1; k<Ox; k++) if (S[k]<best) { best=S[k]; X=k; }
	{
	    double	dx=0.0,a,b,c,den,centerf;
	    if (X>=1 && X<=Ox-2) {
		a=S[X-1]; b=S[X]; c=S[X+1]; den=a-2.0*b+c;
		if (den>0.0) { dx=0.5*(a-c)/den; if (dx<-1.0) dx=-1.0; if (dx>1.0) dx=1.0; }
	    }
	    centerf=((double)(Nx+Ox1+X)+dx)/2.0;
	    Rc=(int)floor(centerf+0.5);
	    (void)fprintf(stderr,"offset axis estimate at layer %ld: %d view pairs, shift %d%+.3f -> centre %.3f, Rc = %d\n",
			  layer,M,Ox1+X,dx,centerf,Rc);
	    if (X==0 || X==Ox-1)
		(void)fprintf(stderr,"warning: the minimum lies at the edge of the search range (centre %d .. %d, as ofct_DO); "
			      "the estimate is unreliable - give center explicitly.\n",(Nx+Ox1)/2,(Nx+Ox2)/2);
	}

	free(S); free(ker); free(tmp); free(H); free(G); free(T2); free(T1); free(I0row); free(D);
	free(V); free(I0); free(v);
	HisClose(s);
	return Rc;
}

int	main(int argc,char **argv)
{
	if (argc<2 || argc>5) {
	    (void)fprintf(stderr,
		"usage : %s layer {center} {pixel size} {offsetangle}\n"
		"        run where a.his, conv.bat and output.log are; a 360-degree\n"
		"        offset scan (angles in output.log) is reconstructed as ofct_rec.\n"
		"        Without center the rotation axis is estimated from this layer\n"
		"        (180-degree scan: ct_rec's search; offset scan: ofct_DO's search\n"
		"        on one layer).\n",argv[0]);
	    return 1;
	}
	if (IsOffsetScan(".",BaseName(argv[0]))) {
	    if (argc<3) {
		char	rc[32],*nargv[4];

		(void)snprintf(rc,sizeof(rc),"%d",EstimateOffsetAxis(".",atol(argv[1])));
		nargv[0]=argv[0]; nargv[1]=argv[1]; nargv[2]=rc; nargv[3]=NULL;
		(void)fprintf(stderr,"%s: offset reconstruction as: ofct_rec %s %s\n",BaseName(argv[0]),argv[1],rc);
		return ofct_rec_main(3,nargv);
	    }
	    return ofct_rec_main(argc,argv);
	}
	return ct_rec_main(argc,argv);
}

/* ====================================================================== */
#elif defined(HIS_SREC)

extern int	hp_tg_main(int argc,char **argv);
extern int	ofct_srec_main(int argc,char **argv);

/* "123" / "123.0" -> 123; anything else is rounded with a warning
   (ofct_srec takes the axis position as an integer) */
static void	IntegerAxis(const char *s,char *out,size_t n)
{
	double	v=atof(s);
	long	r=(long)floor(v+0.5);

	if (fabs(v-(double)r)>1.0e-9)
	    (void)fprintf(stderr,"warning: offset scan: the axis position must be an integer pixel; %s rounded to %ld\n",s,r);
	(void)snprintf(out,n,"%ld",r);
}

int	main(int argc,char **argv)
{
	char	rc[64],rl[128],*oy,*nargv[9];

	if (argc!=6 && argc!=9) {
	    (void)fprintf(stderr,
		"usage : %s HiPic Dr RC RA0 rec\n"
		"usage : %s HiPic Dr L1 C1 L2 C2 RA0 rec\n"
		"        HiPic holds a.his, conv.bat and output.log; run one level above it.\n"
		"        A 360-degree offset scan is reconstructed as ofct_srec with Oy=0\n"
		"        (environment variable HIS_SREC_OY) and RC as the axis position.\n",
		argv[0],argv[0]);
	    return 1;
	}
	if (!IsOffsetScan(argv[1],BaseName(argv[0]))) return hp_tg_main(argc,argv);

	/* ofct_srec: HiPic Rc Oy rangeList Dr RA0 rec */
	if ((oy=getenv("HIS_SREC_OY"))==NULL || *oy=='\0') oy="0";
	if (argc==6) {
	    IntegerAxis(argv[3],rc,sizeof(rc));
	    (void)snprintf(rl,sizeof(rl),"-");
	}
	else {
	    if (atof(argv[4])!=atof(argv[6]))
		Fail("offset scan: a tilted rotation axis (C1 != C2) is not supported by the offset reconstruction; give C1 == C2.");
	    IntegerAxis(argv[4],rc,sizeof(rc));
	    if (atoi(argv[3])<0 || atoi(argv[5])<atoi(argv[3])) Fail("bad layer range L1 L2.");
	    (void)snprintf(rl,sizeof(rl),"%d-%d",atoi(argv[3]),atoi(argv[5]));
	}
	nargv[0]=argv[0]; nargv[1]=argv[1]; nargv[2]=rc; nargv[3]=oy; nargv[4]=rl;
	nargv[5]=argv[2]; nargv[6]=argv[argc-2]; nargv[7]=argv[argc-1]; nargv[8]=NULL;
	(void)fprintf(stderr,"%s: offset reconstruction as: ofct_srec %s %s %s %s %s %s %s\n",
		      BaseName(argv[0]),nargv[1],nargv[2],nargv[3],nargv[4],nargv[5],nargv[6],nargv[7]);
	return ofct_srec_main(8,nargv);
}

/* ====================================================================== */
#elif defined(HIS_DO)

extern int	ofct_DO_main(int argc,char **argv);
extern int	hp2DO_compute(unsigned short **dataD000,unsigned short **dataI000,
			      unsigned short **dataE000,unsigned short **dataD180,
			      unsigned short **dataI180,unsigned short **dataE180,
			      int Nx,int Ny,int z1,int z2,int print_range,
			      double search_percent,const char *out_path,
			      double *A_out,double *B_out);

static unsigned short	**LoadView(HisSet *s,const char *name,int Nx,int Ny)
{
	unsigned short	**rows,*flat;
	int	k=HisResolveName(s,name),y;

	if (k==-2) {
	    (void)fprintf(stderr,"his_DO: %s is not produced by conv.bat.\n",name); exit(1);
	}
	if ((rows=(unsigned short **)malloc(sizeof(unsigned short *)*(size_t)Ny))==NULL ||
	    (flat=(unsigned short *)malloc(sizeof(unsigned short)*(size_t)Ny*(size_t)Nx))==NULL)
	    Fail("his_DO: memory allocation error.");
	for (y=0; y<Ny; y++) rows[y]=flat+(size_t)y*Nx;
	if (k==-1) HisReadDarkBand(s,0,Ny-1,flat); else HisReadQBand(s,k,0,Ny-1,flat);
	return rows;
}

int	main(int argc,char **argv)
{
	HisSet		*s;
	const char	*names[6],*nm[6];
	char		buf[6][32];
	unsigned short	**img[6];
	int		i,Nx,Ny,q_last,n;
	double		A,B,c0,c1;
	const char	*raw;

	if (argc<2 || argc>4) {
	    (void)fprintf(stderr,
		"usage : %s HiPic {stride}            (his_DO,   CPU)\n"
		"        %s HiPic {navg {stride}}     (his_DO_g, GPU)\n"
		"        HiPic holds a.his, conv.bat and output.log; run one level above it.\n"
		"        180-degree scan: hp2DO estimate on the views of the hp2do line of\n"
		"        conv.bat (the extra arguments are ignored); 360-degree offset scan:\n"
		"        ofct_DO.\n",argv[0],argv[0]);
	    return 1;
	}
	raw=argv[1];
	if (IsOffsetScan(raw,BaseName(argv[0]))) return ofct_DO_main(argc,argv);

	/* ---- 180-degree scan: hp2DO on frames from a.his ---- */
	s=HisOpen(raw);
	Nx=HisNx(s); Ny=HisNy(s);
	n=HisHp2doNames(s,names);
	if (n>=6)
	    for (i=0; i<6; i++) nm[i]=names[i];
	else {
	    /* no hp2do line: the convention of ct_rec's centre search */
	    (void)ScanAngles(raw,&q_last);
	    (void)snprintf(buf[0],32,"dark.img");   (void)snprintf(buf[1],32,"q%04d.img",1);
	    (void)snprintf(buf[2],32,"q%04d.img",2); (void)snprintf(buf[3],32,"dark.img");
	    (void)snprintf(buf[4],32,"q%04d.img",q_last+1);
	    (void)snprintf(buf[5],32,"q%04d.img",q_last+2);
	    for (i=0; i<6; i++) nm[i]=buf[i];
	    (void)fprintf(stderr,"his_DO: no hp2do line in conv.bat; using q%d q%d (0 deg) and q%d q%d (180 deg)\n",
			  1,2,q_last+1,q_last+2);
	}
	(void)fprintf(stderr,"his_DO: hp2DO %s %s %s %s %s %s\n",nm[0],nm[1],nm[2],nm[3],nm[4],nm[5]);
	for (i=0; i<6; i++) img[i]=LoadView(s,nm[i],Nx,Ny);

	hp2DO_compute(img[0],img[1],img[2],img[3],img[4],img[5],Nx,Ny,0,Ny-1,0,20.0,NULL,&A,&B);

	/* R0 = A*y + B is minus the axis position: suggest hp_tg-style commands */
	c0=-B; c1=-(A*(double)(Ny-1)+B);
	(void)fprintf(stderr,"axis position: %.3f at layer 0, %.3f at layer %d (tilt %.5f px/layer)\n",c0,c1,Ny-1,-A);
	(void)fprintf(stderr,"try: his_srec_g_c %s 1.0 0 %.1f %d %.1f 0.0 rec   (replace 1.0 by the pixel size Dr)\n",raw,c0,Ny-1,c1);
	(void)fprintf(stderr,"try: his_rec_g_c %d %.1f 1.0 0.0           (in %s)\n",Ny/2,-(A*(double)(Ny/2)+B),raw);

	for (i=0; i<6; i++) { free(img[i][0]); free(img[i]); }
	HisClose(s);

	{   /* record command to cmd-hst.log */
	    FILE *f_log; int i_log;
	    if ((f_log=fopen("cmd-hst.log","a"))!=NULL) {
		for (i_log=0; i_log<argc; ++i_log) fprintf(f_log,"%s ",argv[i_log]);
		fprintf(f_log,"   %% hp2DO axis %.3f / %.3f\n",c0,c1); fclose(f_log);
	    }
	}
	return 0;
}

#else
#error "build with -DHIS_REC, -DHIS_SREC or -DHIS_DO"
#endif
