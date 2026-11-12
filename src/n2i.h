/*
 * n2i.h - the Noise2Inverse front end shared by his_srec_g_F_n2i
 *         (his_main.c, -DHIS_N2I), hp_tg_g_F_n2i and ofct_srec_g_F_n2i
 *         (n2i_main.c).  Static functions, included by one file per
 *         program (as axis1.h).
 *
 * N2iParent() takes the command line of the body program plus a trailing
 * K (2 or 4), splits the sample views listed in <HiPic>/output.log into K
 * interleaved subsets (k, k+K, ...) in angle order and runs the body once
 * per subset into rec/0 .. rec/K-1, RA0 shifted by the angle of the
 * subset's first view so that the K volumes share one geometry.  Every
 * subset is a child process of the same program (the bodies keep static
 * state for one run), started with the environment variables
 *   RHIS_VIEWS=k/K     the reader keeps every K-th sample view starting
 *                      at the k-th (rhp_views.h: rhis_hp.c and rhp_c.c)
 *   RHIS_N2I_CHILD=1   the child skips this front end and runs the body
 * RA0 and rec must be the two arguments before K; the number of sample
 * views must be a multiple of K, on a 360-degree offset scan of 2K (so
 * that the opposing views m and m+Nt/2 fall into the same subset).
 */
#ifndef N2I_H
#define N2I_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef WINDOWS
#include <process.h>
#include <direct.h>
#define N2I_MKDIR(p)	_mkdir(p)
#define N2I_SETENV(n,v)	_putenv_s(n,v)
#else
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#define N2I_MKDIR(p)	mkdir(p,0777)
#define N2I_SETENV(n,v)	setenv(n,v,1)
#endif

#define N2I_LEN			2048
#define N2I_OFFSET_ANGLE	350.0	/* as HIS_OFFSET_ANGLE of his_main.c */

typedef struct {
	const char	*body;		/* name of the body in the usage text (his_srec, hp_tg, ofct_srec) */
	const char	*form1,*form2;	/* the two argument lists of the body, without K */
	int		argc1,argc2;	/* argc of the two forms, K included */
	int		offset;		/* 1: 360-degree offset scan (2K), 0: 180-degree scan (K) */
	int		(*detect)(const char *dir,const char *prog);
					/* if not NULL: decide offset from <dir>/output.log
					   (his_main.c: IsOffsetScan) */
} N2iSpec;

static const char	*N2iBaseName(const char *p)
{
	const char	*s=strrchr(p,'/'),*b=strrchr(p,'\\');

	if (b!=NULL && (s==NULL || b>s)) s=b;
	return (s!=NULL)?s+1:p;
}

static void	N2iFail(const char *msg)
{
	(void)fprintf(stderr,"%s\n",msg); exit(1);
}

/* angles of the sample views listed before the first blank line of
   <dir>/output.log (or RHP_O), sorted ascending; returns their number */
static int	N2iSampleAngles(const char *dir,double **angles)
{
	char	path[N2I_LEN],str[N2I_LEN],*env;
	FILE	*f;
	int	q,it,n=0,cap=1024,i,j;
	double	c,a,*v=(double *)malloc(sizeof(double)*(size_t)cap),t;

	if (v==NULL) N2iFail("no memory.");
	if ((env=getenv("RHP_O"))!=NULL && *env!='\0') (void)snprintf(path,N2I_LEN,"%s",env);
	else if (snprintf(path,N2I_LEN,"%s/output.log",dir)>=N2I_LEN) N2iFail("output.log : path too long.");
	if ((f=fopen(path,"r"))==NULL) { (void)fprintf(stderr,"%s : file not open.\n",path); exit(1); }
	while (fgets(str,N2I_LEN,f)!=NULL) {
	    if (sscanf(str,"%d %lf %lf %d",&q,&c,&a,&it)!=4) break;	/* ONLY_CT_VIEWS */
	    if (it!=1) continue;
	    if (n==cap && (v=(double *)realloc(v,sizeof(double)*(size_t)(cap*=2)))==NULL) N2iFail("no memory.");
	    v[n++]=a;
	}
	(void)fclose(f);
	for (i=1; i<n; i++) { t=v[i]; for (j=i; j>0 && v[j-1]>t; j--) v[j]=v[j-1]; v[j]=t; }
	*angles=v;
	return n;
}

/* run this program once per subset as a child process */
static int	N2iRunChild(char **args)
{
#ifdef WINDOWS
	char	**qa;
	int	i,n,st;

	for (n=0; args[n]!=NULL; n++) ;
	if ((qa=(char **)malloc(sizeof(char *)*(size_t)(n+1)))==NULL) N2iFail("no memory.");
	for (i=0; i<n; i++) {		/* the CRT joins the arguments with blanks and does not quote */
	    if (strchr(args[i],' ')!=NULL) {
		size_t	l=strlen(args[i]);
		if ((qa[i]=(char *)malloc(l+3))==NULL) N2iFail("no memory.");
		qa[i][0]='"'; memcpy(qa[i]+1,args[i],l); qa[i][l+1]='"'; qa[i][l+2]='\0';
	    }
	    else qa[i]=args[i];
	}
	qa[n]=NULL;
	st=(int)_spawnvp(_P_WAIT,args[0],(const char *const *)qa);
	return st;
#else
	pid_t	pid=fork();
	int	st=0;

	if (pid<0) N2iFail("fork failed.");
	if (pid==0) { execvp(args[0],args); perror(args[0]); _exit(127); }
	if (waitpid(pid,&st,0)<0) N2iFail("waitpid failed.");
	return WIFEXITED(st)?WEXITSTATUS(st):128;
#endif
}

static int	N2iParent(int argc,char **argv,const N2iSpec *sp)
{
	char	views[32],ra[32],rec[N2I_LEN],*args[16];
	const char	*prog=N2iBaseName(argv[0]);
	double	*a,ra0;
	int	K,k,Nt,i,n,ira,irec,st,offset;

	if (argc!=sp->argc1 && argc!=sp->argc2) {
	    (void)fprintf(stderr,
		"usage : %s %s K\n"
		"usage : %s %s K\n"
		"        Noise2Inverse: the sample views are split into K (2 or 4)\n"
		"        interleaved subsets and each is reconstructed as %s into\n"
		"        rec/0 .. rec/K-1, RA0 shifted by the angle of the subset's first\n"
		"        view so that the K volumes share one geometry.  The number of\n"
		"        sample views must be a multiple of %s.\n",
		argv[0],sp->form1,argv[0],sp->form2,sp->body,
		(sp->detect!=NULL)?"K (of 2K on a 360-degree offset\n"
		"        scan, so that opposing views fall into the same subset)":
		sp->offset?"2K (the opposing views m and m+Nt/2\n"
		"        of the offset scan must fall into the same subset)":"K");
	    return 1;
	}
	K=atoi(argv[argc-1]);
	if (K!=2 && K!=4) N2iFail("K must be 2 or 4.");
	offset=(sp->detect!=NULL)?(*sp->detect)(argv[1],prog):sp->offset;
	ira=argc-3; irec=argc-2;		/* RA0 and rec are the two arguments before K */
	Nt=N2iSampleAngles(argv[1],&a);
	if (sp->detect==NULL && Nt>0 && (a[Nt-1]>=N2I_OFFSET_ANGLE)!=offset)
	    (void)fprintf(stderr,"%s: warning: the sample views reach %.1f degrees, but %s reconstructs a %s scan.\n",
			  prog,a[Nt-1],sp->body,offset?"360-degree offset":"180-degree");
	if (Nt%K!=0 || (offset && Nt%(2*K)!=0)) {
	    (void)fprintf(stderr,"%s: %d sample views are not a multiple of %d%s.\n",prog,Nt,
			  offset?2*K:K,offset?" (2K: the opposing views of an offset scan must fall into the same subset)":"");
	    return 1;
	}
	ra0=atof(argv[ira]);
	(void)N2I_MKDIR(argv[irec]);
	(void)fprintf(stderr,"%s: %d sample views, K=%d -> %d views per subset, angle step %.5f deg\n",
		      prog,Nt,K,Nt/K,(Nt>1)?(a[1]-a[0])*K:0.0);
	for (k=0; k<K; k++) {
	    (void)snprintf(views,sizeof(views),"%d/%d",k,K);
	    (void)snprintf(ra,sizeof(ra),"%.6f",ra0+(a[k]-a[0]));
	    if (snprintf(rec,N2I_LEN,"%s/%d",argv[irec],k)>=N2I_LEN) N2iFail("rec : path too long.");
	    (void)N2I_MKDIR(rec);
	    N2I_SETENV("RHIS_VIEWS",views); N2I_SETENV("RHIS_N2I_CHILD","1");
	    for (n=0, i=0; i<argc-1; i++) args[n++]=(i==ira)?ra:(i==irec)?rec:argv[i];
	    args[n]=NULL;
	    (void)fprintf(stderr,"==== subset %d / %d: RA0 %s -> %s ====\n",k,K,ra,rec);
	    if ((st=N2iRunChild(args))!=0) {
		(void)fprintf(stderr,"%s: subset %d failed (status %d).\n",prog,k,st);
		free(a); return st;
	    }
	}
	free(a);
	(void)fprintf(stderr,"%s: done, %d subsets in %s/0 .. %s/%d\n",prog,K,argv[irec],argv[irec],K-1);
	return 0;
}

#endif	/* N2I_H */
