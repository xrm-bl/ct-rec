/*
 * rhp_c.c - Combined HiPic(.img) / TIFF(.tif) reader.
 *
 * Drop-in replacement for rhp.c that auto-detects the input format:
 *   - if dark.img exists in the directory, read everything as HiPic .img
 *   - else if dark.tif exists, read dark and every q-image as 16-bit TIFF
 *
 * The public API (InitReadHiPic / ReadHiPic / TermReadHiPic) and the HiPic
 * struct are unchanged, so this links in place of rhp.c. libtiff is required
 * (already linked by the hp_tg_* programs).
 *
 * Environment overrides (same as rhp.c):
 *   RHP_O  output.log path   ("-" = stdin)
 *   RHP_D  dark frame name    (format inferred from its extension)
 *   RHP_Q  q-image prefix char (default 'q')
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef	_WIN32
#include "msdirent.h"
#else
#include <dirent.h>
#endif
#include <stdint.h>   /* SIZE_MAX, uint16_t, uint32_t */
#include "tiffio.h"
#include "rhp.h"

/* selected input format: 0 = HiPic .img, 1 = TIFF .tif. Set in
 * InitReadHiPic, read by ReadHiPic (one dataset per process). */
static int	useTiff = 0;

/* セキュリティ修正: 要素数の積を size_t で計算し，オーバーフロー時は停止する */
static size_t safe_count(size_t a, size_t b, size_t c)
{
	if (b != 0 && a > SIZE_MAX / b) goto over;
	a *= b;
	if (c != 0 && a > SIZE_MAX / c) goto over;
	return a * c;
over:
	(void)fprintf(stderr, "size overflow in image allocation.\n"); exit(1);
}

#define EPS	1e-9	/* for time interval */

static void	Error(char *dir,char *name,char *msg)
{
	if (*dir!='\0') (void)fprintf(stderr,"%s/",dir);

	(void)fprintf(stderr,"%s : %s\n",name,msg); exit(1);
}

static int	StrCmp(char *s1,char *s2)
{
	int	c;

	while ((c=(*s1>='A' && *s1<='Z')?*s1|32:*s1)==(*s2)) {
	    if (c=='\0') return 0;

	    ++s1; ++s2;
	}
	return 1;
}

#define LEN	2048
#define	WORD	unsigned short

/* ---------------------------------------------------------------------- */
/* HiPic .img primitives                                                  */
/* ---------------------------------------------------------------------- */

static int	GetW(FILE *file)
{
	int	lo,hi;

	return ((lo=fgetc(file))==EOF || (hi=fgetc(file))==EOF)?-1:hi*256+lo;
}

static FILE	*OpenImg(char *dir,char *name,int *Nx,int *Ny)
{
	char	path[LEN];
	FILE	*file;
        int	cl;

	if (snprintf(path,LEN,"%s/%s",dir,name) >= LEN) Error(dir,name,"path too long.");
	if ((file=fopen(path,"rb"))==NULL) Error(dir,name,"file not open.");

	if (fgetc(file)!='I' ||
	    fgetc(file)!='M') Error(dir,name,"bad magic number.");

	if ((cl=GetW(file))<0) Error(dir,name,"bad comment length.");

	if ((*Nx=GetW(file))<0 ||
	    (*Ny=GetW(file))<0) Error(dir,name,"bad image size.");

	(void)fseek(file,4L,SEEK_CUR);

	if (GetW(file)!=2) Error(dir,name,"bad image type.");

	(void)fseek(file,50L+cl,SEEK_CUR);

        return file;
}

static void	ReadImg(char *dir,char *name,int Nx,int Ny,WORD *W)
{
	int	x,y,i;
	FILE	*file=OpenImg(dir,name,&x,&y);
	short	s=127;
	char	*c=(char *)&s,
		*dre="data read error.";

	if (Nx!=x || Ny!=y) Error(dir,name,"image size not match.");

	if (*c!=s)
	    for (y=0; y<Ny; y++)
	    for (x=0; x<Nx; x++)
		if ((i=GetW(file))>=0) *(W++)=(WORD)i; else Error(dir,name,dre);
	else
	    if (fread(W,sizeof(WORD)*Ny*Nx,1,file)!=1) Error(dir,name,dre);

	(void)fclose(file);
}

/* ---------------------------------------------------------------------- */
/* TIFF .tif primitives (16-bit monochrome)                               */
/* ---------------------------------------------------------------------- */

static TIFF	*OpenTiff(char *dir,char *name,int *Nx,int *Ny)
{
	char		path[LEN];
	TIFF		*tif;
	uint32_t	width,height;
	uint16_t	bps,spp,photo;

	if (snprintf(path,LEN,"%s/%s",dir,name) >= LEN) Error(dir,name,"path too long.");
	if ((tif=TIFFOpen(path,"r"))==NULL) Error(dir,name,"file not open.");

	TIFFGetField(tif,TIFFTAG_IMAGEWIDTH,&width);
	TIFFGetField(tif,TIFFTAG_IMAGELENGTH,&height);

	*Nx=(int)width;
	*Ny=(int)height;

	TIFFGetField(tif,TIFFTAG_BITSPERSAMPLE,&bps);
	TIFFGetField(tif,TIFFTAG_SAMPLESPERPIXEL,&spp);
	TIFFGetField(tif,TIFFTAG_PHOTOMETRIC,&photo);

	if (bps!=16 || spp!=1 ||
	    (photo!=PHOTOMETRIC_MINISBLACK && photo!=PHOTOMETRIC_MINISWHITE))
	    Error(dir,name,"not a 16-bit monochrome TIFF image.");

	return tif;
}

static void	ReadTiff(char *dir,char *name,int Nx,int Ny,WORD *W)
{
	int		x,y;
	TIFF		*tif=OpenTiff(dir,name,&x,&y);
	char		*dre="data read error.";
	tsize_t		scanline_size;
	tdata_t		buf;

	if (Nx!=x || Ny!=y) Error(dir,name,"image size not match.");

	scanline_size=TIFFScanlineSize(tif);
	if ((buf=_TIFFmalloc(scanline_size))==NULL)
	    Error(dir,name,"no memory for TIFF buffer.");

	for (y=0; y<Ny; y++) {
	    if (TIFFReadScanline(tif,buf,y,0)<0) Error(dir,name,dre);
	    memcpy(W+(size_t)y*Nx,buf,(size_t)Nx*sizeof(WORD));
	}

	_TIFFfree(buf);
	TIFFClose(tif);
}

/* ---------------------------------------------------------------------- */
/* format-dispatching wrappers                                            */
/* ---------------------------------------------------------------------- */

/* read the header only, to obtain the image size */
static void	GetSize(char *dir,char *name,int *Nx,int *Ny)
{
	if (useTiff) {
	    TIFF *tif=OpenTiff(dir,name,Nx,Ny); TIFFClose(tif);
	} else {
	    FILE *file=OpenImg(dir,name,Nx,Ny); (void)fclose(file);
	}
}

static void	Read(char *dir,char *name,int Nx,int Ny,WORD *W)
{
	if (useTiff) ReadTiff(dir,name,Nx,Ny,W);
	else	     ReadImg(dir,name,Nx,Ny,W);
}

/* ---------------------------------------------------------------------- */

static int	Compare(OutputLog *OL1,OutputLog *OL2)
{
	int	it=OL2->it-OL1->it;

	return (it)?it:(OL1->it==0)?(OL1->c>OL2->c)?1:(OL1->c<OL2->c)?-1:0:
		       (OL1->it==1)?(OL1->a>OL2->a)?1:(OL1->a<OL2->a)?-1:0:0;
}

void	InitReadHiPic(char *dir,HiPic *hp)
{
	char		*env,q_img[LEN],str[LEN],
			output_log[LEN]="",
			dark_img[LEN]="",
			dark_tif[LEN]="",
			dark[LEN]="",
			*fnf="file not found.",
			*nmfdd="no memory for directory data.",
			*bsn="bad sequence number.",
			*dsn="duplicated sequence number.";
	DIR		*Dir;
	struct dirent	*sd;
	int		l,q,it,x,y,i,darkFromEnv=0;
	const char	*ext;
	char		qchar;
	FILE		*file;
	double		c,a,c1,c2;
	OutputLog	*OL;

	if ((env=getenv("RHP_O"))!=NULL){	/* 修正: 長さ検査付きコピー */
	    if (strlen(env) >= LEN) Error("","RHP_O","environment value too long.");
	    (void)snprintf(output_log,LEN,"%s",env);
	}
	if ((env=getenv("RHP_D"))!=NULL){
	    if (strlen(env) >= LEN) Error("","RHP_D","environment value too long.");
	    (void)snprintf(dark,LEN,"%s",env);
	    darkFromEnv=1;
	}

	qchar = ((env=getenv("RHP_Q"))!=NULL && *env!='\0' &&
		 (*env|32)>='a' && (*env|32)<='z') ? (char)(*env|32) : 'q';

	hp->Nq=0;

	if ((Dir=opendir(dir))==NULL) Error("",dir,"directory not open.");

#define NAME	sd->d_name

	/* pass 1: locate output.log and the dark frame (.img preferred, .tif fallback) */
	while ((sd=readdir(Dir))!=NULL)
	    if (*output_log=='\0' && !StrCmp(NAME,"output.log")) {
		if (snprintf(output_log,LEN,"%s/%s",dir,NAME) >= LEN) Error(dir,NAME,"path too long.");
	    } else if (!darkFromEnv && *dark_img=='\0' && !StrCmp(NAME,"dark.img"))
		(void)snprintf(dark_img,LEN,"%s",NAME);
	    else if (!darkFromEnv && *dark_tif=='\0' && !StrCmp(NAME,"dark.tif"))
		(void)snprintf(dark_tif,LEN,"%s",NAME);

	if (*output_log=='\0') Error(dir,"output.log",fnf);

	/* decide the input format */
	if (darkFromEnv) {
	    l=(int)strlen(dark);
	    useTiff = (l>=4 && !StrCmp(dark+l-4,".tif")) ||
		      (l>=5 && !StrCmp(dark+l-5,".tiff"));
	} else if (*dark_img!='\0') {
	    useTiff=0; (void)snprintf(dark,LEN,"%s",dark_img);
	} else if (*dark_tif!='\0') {
	    useTiff=1; (void)snprintf(dark,LEN,"%s",dark_tif);
	} else
	    Error(dir,"dark.img/dark.tif",fnf);

	ext = useTiff ? ".tif" : ".img";

	(void)snprintf(q_img,LEN,"%c[0-9]*[0-9]%s",qchar,ext);

	/* pass 2: count q-images with the selected extension */
	rewinddir(Dir);

	while ((sd=readdir(Dir))!=NULL)
	    if ((NAME[0]|32)==qchar &&
		(l=(int)strlen(NAME))>5 &&
		!StrCmp(NAME+l-4,(char *)ext) &&
		sscanf(NAME+1,"%d",&q)==1 &&
		q>=hp->Nq)
		hp->Nq=q+1;

	if (hp->Nq==0) Error(dir,q_img,fnf);

#define MALLOC(type,noe)	(type *)malloc(sizeof(type)*(noe))

	if ((hp->dir=strdup(dir))==NULL ||
	    (hp->q_img=MALLOC(char *,hp->Nq))==NULL) Error("",dir,nmfdd);

	/* pass 3: record q-image file names */
	rewinddir(Dir);

	for (q=0; q<hp->Nq; q++) hp->q_img[q]=NULL;

	while ((sd=readdir(Dir))!=NULL)
	    if ((NAME[0]|32)==qchar &&
		(l=(int)strlen(NAME))>5 &&
		!StrCmp(NAME+l-4,(char *)ext) &&
		sscanf(NAME+1,"%d",&q)==1 &&
		q>=0) {
		if (q>=hp->Nq) Error(dir,NAME,bsn);

                if (hp->q_img[q]!=NULL) Error(dir,NAME,dsn);

                if ((hp->q_img[q]=strdup(NAME))==NULL) Error(dir,NAME,nmfdd);
	    }

	(void)closedir(Dir);

	if ((hp->OL=MALLOC(OutputLog,hp->Nq+1))==NULL)
	    Error("",output_log,"no memory for log data.");

	for (q=0; q<=hp->Nq; q++) hp->OL[q].it=(-1);

	GetSize(dir,dark,&(hp->Nx),&(hp->Ny));

	hp->Ni=hp->Nt=0;

	if (*output_log=='-')
	    file=stdin;
	else
	    if ((file=fopen(output_log,"r"))==NULL)
		Error("",output_log,"file not open.");

	while (fgets(str,LEN,file)!=NULL)
	    if ((sscanf(str,"%d %lf %lf %d%n",&q,&c,&a,&it,&l)==4)) {
		str[l]='\0';

		if (q<0 || q>=hp->Nq || hp->q_img[q]==NULL) Error("",str,bsn);

		if (it!=0 && it!=1) Error("",str,"unacceptable log data.");

		GetSize(dir,hp->q_img[q],&x,&y);

		if (hp->Nx!=x || hp->Ny!=y)
		    Error(dir,hp->q_img[q],"image size not match.");

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

	if (hp->Ni==0) Error(dir,q_img,"no file referred as I0-image.");
	if (hp->Nt==0) Error(dir,q_img,"no file referred as I-image.");

	if ((hp->D=MALLOC(WORD,safe_count((size_t)hp->Ny,(size_t)hp->Nx,1)))==NULL ||
	    (hp->I=MALLOC(WORD *,hp->Ni+(hp->Ni==1)))==NULL ||
	    (*(hp->I)=MALLOC(WORD,safe_count((size_t)hp->Ni,(size_t)hp->Ny,(size_t)hp->Nx)))==NULL ||
	    (hp->T=MALLOC(FOM *,hp->Ny))==NULL ||
	    (*(hp->T)=MALLOC(FOM,safe_count((size_t)hp->Ny,(size_t)hp->Nx,1)))==NULL)
	    Error("",dir,"no memory for image data.");

	Read(dir,dark,hp->Nx,hp->Ny,hp->D);

	qsort(hp->OL,hp->Nq+1,sizeof(OutputLog),(int (*)(const void *,const void *))Compare);

	OL=hp->OL+hp->Nt;
	for (i=0; i<hp->Ni; i++) {
	    if (i>0) {
	    	if (OL[i-1].c+EPS>=OL[i].c){
	    		fprintf(stderr, "%d\t%lf\t%lf\n",i,OL[i-1].c+EPS,OL[i].c );
	//	    	Error("",output_log,"bad time interval.");

	    	}

		hp->I[i]=hp->I[i-1]+(size_t)hp->Ny*hp->Nx;
	    }
	    Read(dir,hp->q_img[OL[i].q],hp->Nx,hp->Ny,hp->I[i]);
	}
	if (hp->Ni==1) {
	    hp->Ni=2; OL[1].q=OL[0].q;
		      OL[0].c=c1-EPS; OL[1].c=c2+EPS;
		      OL[1].a=OL[0].a; OL[1].it=OL[0].it; hp->I[1]=hp->I[0];
	}
	for (y=1; y<hp->Ny; y++) hp->T[y]=hp->T[y-1]+hp->Nx;
}

void	ReadHiPic(HiPic *hp,int t)
{
	char		str[LEN];
	OutputLog	*ol,*OL;
	int		i,y,x;
	WORD		*D,*I1,*I2,*T;
	double		r1,r2,I_D,T_D;
	double		td_sum,black_thresh;
	char		*env_bt;

	if (t<0 || t>=hp->Nt) {
	    (void)sprintf(str,"%d",t); Error("",str,"bad sequence number.");
	}

	env_bt = getenv("CT_REC_BLACK_THRESH");
	black_thresh = (env_bt != NULL) ? atof(env_bt) : 1.0;

	ol=hp->OL+t;
	Read(hp->dir,hp->q_img[ol->q],hp->Nx,hp->Ny,(WORD *)*(hp->T));

	OL=hp->OL+hp->Nt;
	for (i=hp->Ni-2; i>0; i--) if (ol->c>=OL[i].c) break;

	D=hp->D+(size_t)hp->Ny*hp->Nx;
	I1=hp->I[i]+(size_t)hp->Ny*hp->Nx;
	I2=hp->I[i+1]+(size_t)hp->Ny*hp->Nx;
	T=(WORD *)*(hp->T)+(size_t)hp->Ny*hp->Nx;
	r1=1.0-(r2=(ol->c-OL[i].c)/(OL[i+1].c-OL[i].c));

	td_sum = 0.0;
	for (y=hp->Ny-1; y>=0; y--)
	for (x=hp->Nx-1; x>=0; x--) {
	    --D; --I1; --I2; --T;
	    T_D=(double)*T-(double)*D;
	    I_D=r1*(double)*I1+r2*(double)*I2-(double)*D;
	    hp->T[y][x]=(I_D>0.0 && T_D>0.0)?T_D/I_D:ERROR_VALUE;
	    td_sum += T_D;
	}

	/* black check: average of (T-dark) over all pixels */
	if (td_sum / (double)(hp->Nx * hp->Ny) < black_thresh){
	    (void)fprintf(stderr,
	        "Warning\t black\t t=%d avg=%.2f (thresh=%.2f)\n",
	        t, td_sum/(double)(hp->Nx*hp->Ny), black_thresh);
	    for (y=0; y<hp->Ny; y++)
	    for (x=0; x<hp->Nx; x++)
	        hp->T[y][x]=ERROR_VALUE;
	}
}

void	TermReadHiPic(HiPic *hp)
{
	int	q;

	free(*(hp->T)); free(hp->T);
	free(*(hp->I)); free(hp->I); free(hp->D);

	free(hp->OL);

	for (q=0; q<hp->Nq; q++) if (hp->q_img[q]!=NULL) free(hp->q_img[q]);

	free(hp->q_img); free(hp->dir);
}
