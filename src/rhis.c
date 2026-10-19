/*
 * rhis.c - direct reader for a.his + conv.bat + output.log (see rhis.h).
 *
 * This file is the HisSet layer: frame table, conv.bat interpretation,
 * band reads with img_ave-compatible averaging.  The rhp.h API on top of
 * it (InitReadHiPic / ReadHiPic / ReadHiPicBand / TermReadHiPic, used by
 * hp_tg_ku.c, ofct_srec.c and ofct_DO.c) lives in rhis_hp.c, which is
 * compiled once per FOM type (see rhp.h, RHP_NAME_SUFFIX).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdint.h>

#ifdef	_WIN32
#include <windows.h>
typedef HANDLE	HisFd;
#define HIS_BAD_FD	INVALID_HANDLE_VALUE
#else
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
typedef int	HisFd;
#define HIS_BAD_FD	(-1)
#endif

#ifdef	_OPENMP
#include <omp.h>
#endif

#include "rhis.h"

#define LEN		2048
#define HIS_HDR		64
#define WORD		unsigned short

/* ---------------------------------------------------------------------- */
/* data structures                                                        */
/* ---------------------------------------------------------------------- */

enum { HQ_NONE=0, HQ_FRAME=1, HQ_AVG=2 };

typedef struct {
	int	kind;		/* HQ_* */
	int	frame;		/* HQ_FRAME: 0-based frame index */
	int	nlist;		/* HQ_AVG: number of frames */
	int	*list;		/* HQ_AVG: 0-based frame indices */
	/* band cache for HQ_AVG (the dark and the I0 groups are asked for
	   repeatedly; a 200-frame average is worth keeping) */
	int	cy1,cy2;
	WORD	*cache;
} HisQ;

struct HisSet {
	char	*dir,*path,*conv;
	HisFd	fd;
	int	Nx,Ny,nframes,Nq;
	unsigned long long	*off;	/* data offset of every frame */
	HisQ	dark,*Q;
	int	n_avg,n_frame;		/* statistics for the start-up message */
	int	frames_per_scan,nscan,scan,frame0;	/* Z scans, see HisScanLayout */
	char	*hp2do[6];		/* image names of the hp2do line (his_DO) */
	int	n_hp2do;
};

/* ---------------------------------------------------------------------- */
/* errors                                                                 */
/* ---------------------------------------------------------------------- */

static void	HisError(const char *where,const char *msg)
{
	if (where!=NULL && *where!='\0') (void)fprintf(stderr,"%s : ",where);
	(void)fprintf(stderr,"%s\n",msg); exit(1);
}

static void	*HisMalloc(size_t n,const char *what)
{
	void	*p=malloc(n?n:1);

	if (p==NULL) {
	    (void)fprintf(stderr,"rhis: no memory for %s (%.0f bytes).\n",
			  what,(double)n);
	    exit(1);
	}
	return p;
}

/* ---------------------------------------------------------------------- */
/* positional file reads (thread-safe on a shared handle)                 */
/* ---------------------------------------------------------------------- */

static HisFd	HisFileOpen(const char *path)
{
#ifdef	_WIN32
	return CreateFileA(path,GENERIC_READ,FILE_SHARE_READ,NULL,
			   OPEN_EXISTING,FILE_ATTRIBUTE_NORMAL,NULL);
#else
	return open(path,O_RDONLY);
#endif
}

static void	HisFileClose(HisFd fd)
{
#ifdef	_WIN32
	(void)CloseHandle(fd);
#else
	(void)close(fd);
#endif
}

static unsigned long long	HisFileSize(HisFd fd)
{
#ifdef	_WIN32
	LARGE_INTEGER	sz;

	if (!GetFileSizeEx(fd,&sz)) return 0ull;
	return (unsigned long long)sz.QuadPart;
#else
	off_t	sz=lseek(fd,0,SEEK_END);

	return (sz<0)?0ull:(unsigned long long)sz;
#endif
}

/* read exactly n bytes at offset off; returns 0 on success */
static int	HisPread(HisFd fd,void *buf,size_t n,unsigned long long off)
{
	char	*p=(char *)buf;

	while (n>0) {
#ifdef	_WIN32
	    OVERLAPPED	ov;
	    DWORD	want=(n>0x40000000u)?0x40000000u:(DWORD)n,got=0;

	    memset(&ov,0,sizeof(ov));
	    ov.Offset=(DWORD)(off&0xffffffffu);
	    ov.OffsetHigh=(DWORD)(off>>32);
	    if (!ReadFile(fd,p,want,&got,&ov) || got==0) return -1;
#else
	    size_t	want=(n>0x40000000u)?0x40000000u:n;
	    ssize_t	got=pread(fd,p,want,(off_t)off);

	    if (got<=0) return -1;
#endif
	    p+=got; n-=(size_t)got; off+=(unsigned long long)got;
	}
	return 0;
}

/* ---------------------------------------------------------------------- */
/* his file: header scan                                                  */
/* ---------------------------------------------------------------------- */

static unsigned	GetU16(const unsigned char *b)
{
	return (unsigned)b[0]|((unsigned)b[1]<<8);
}

static int	HostIsLittleEndian(void)
{
	WORD	w=1;

	return *(unsigned char *)&w==1;
}

static void	HisScanHeaders(HisSet *s)
{
	unsigned char	hdr[HIS_HDR];
	unsigned long long	pos=0,size=HisFileSize(s->fd),frame_bytes;
	unsigned	cl,w,h,type;
	int	k;

	if (HisPread(s->fd,hdr,HIS_HDR,0)) HisError(s->path,"cannot read the his header.");
	if (hdr[0]!='I' || hdr[1]!='M') HisError(s->path,"bad magic number (not a his file).");

	w=GetU16(hdr+4); h=GetU16(hdr+6); type=GetU16(hdr+12);
	s->nframes=(int)(GetU16(hdr+14)+65536u*GetU16(hdr+16));

	if (w==0 || h==0) HisError(s->path,"bad image size.");
	if (type!=2) HisError(s->path,"only 16-bit (type 2) his files are supported.");
	if (s->nframes<=0) HisError(s->path,"bad number of frames.");

	s->Nx=(int)w; s->Ny=(int)h;
	frame_bytes=(unsigned long long)w*(unsigned long long)h*2ull;

	s->off=(unsigned long long *)HisMalloc(sizeof(unsigned long long)*(size_t)s->nframes,"frame table");

	/* every frame has its own 64-byte header and comment; walk them */
	for (k=0; k<s->nframes; k++) {
	    if (k>0) {
		if (HisPread(s->fd,hdr,HIS_HDR,pos)) {
		    (void)fprintf(stderr,"rhis: %s : frame %d header unreadable "
				  "(file truncated? %d frames expected).\n",
				  s->path,k+1,s->nframes);
		    exit(1);
		}
		if (GetU16(hdr+4)!=w || GetU16(hdr+6)!=h || GetU16(hdr+12)!=type) {
		    (void)fprintf(stderr,"rhis: %s : frame %d has a different "
				  "size/type (%ux%u type %u vs %ux%u type %u).\n",
				  s->path,k+1,GetU16(hdr+4),GetU16(hdr+6),
				  GetU16(hdr+12),w,h,type);
		    exit(1);
		}
	    }
	    cl=GetU16(hdr+2);
	    s->off[k]=pos+HIS_HDR+cl;
	    pos=s->off[k]+frame_bytes;
	    if (size!=0ull && pos>size) {
		(void)fprintf(stderr,"rhis: %s : frame %d runs past the end of the "
			      "file (%d frames expected).\n",s->path,k+1,s->nframes);
		exit(1);
	    }
	}
}

/* rows y1..y2 of frame f -> dst (host byte order) */
static void	HisReadFrameBand(HisSet *s,int f,int y1,int y2,WORD *dst)
{
	size_t	n=(size_t)(y2-y1+1)*(size_t)s->Nx;

	if (HisPread(s->fd,dst,n*sizeof(WORD),
		     s->off[f]+(unsigned long long)y1*(unsigned long long)s->Nx*2ull)) {
	    (void)fprintf(stderr,"rhis: %s : read error in frame %d rows %d-%d.\n",
			  s->path,f+1,y1,y2);
	    exit(1);
	}
	if (!HostIsLittleEndian()) {
	    unsigned char	*b=(unsigned char *)dst,t;
	    size_t		i;

	    for (i=0; i<n; i++) { t=b[2*i]; b[2*i]=b[2*i+1]; b[2*i+1]=t; }
	}
}

/* ---------------------------------------------------------------------- */
/* conv.bat interpretation                                                */
/* ---------------------------------------------------------------------- */

static int	StrCaseEq(const char *a,const char *b)
{
	while (*a && *b && tolower((unsigned char)*a)==tolower((unsigned char)*b)) { a++; b++; }
	return *a=='\0' && *b=='\0';
}

static int	EndsWithImg(const char *tok,size_t *stem_len)
{
	size_t	l=strlen(tok);

	if (l<5 || !StrCaseEq(tok+l-4,".img")) return 0;
	*stem_len=l-4;
	return 1;
}

/* <letter><digits>.img -> letter (lower case) and number; 0 if not */
static int	ParseImgName(const char *tok,char *letter,long *num)
{
	size_t	l,i;
	long	n=0;

	if (!EndsWithImg(tok,&l) || l<2 || !isalpha((unsigned char)tok[0])) return 0;
	for (i=1; i<l; i++) {
	    if (!isdigit((unsigned char)tok[i])) return 0;
	    n=n*10+(tok[i]-'0');
	    if (n>100000000L) return 0;
	}
	*letter=(char)tolower((unsigned char)tok[0]); *num=n;
	return 1;
}

static int	FrameIndex(HisSet *s,const char *tok,int lineno)
{
	char	c;
	long	n;

	if (!ParseImgName(tok,&c,&n) || c=='q') {
	    (void)fprintf(stderr,"rhis: %s line %d : '%s' is not a frame name "
			  "(<letter><digits>.img).\n",s->conv,lineno,tok);
	    exit(1);
	}
	if (n<1 || n-1+s->frame0>=s->nframes) {
	    (void)fprintf(stderr,"rhis: %s line %d : frame %s is outside the his "
			  "file (%d frames%s).\n",s->conv,lineno,tok,s->nframes,
			  (s->nscan>1)?", scan selected by RHIS_SCAN":"");
	    exit(1);
	}
	return (int)(n-1)+s->frame0;
}

/* ---- Z scans: several scans of equal length in one his file ----
   conv.bat names the frames of ONE scan (a0001.img .. a<L>.img, L = the
   highest frame number it uses).  A Z scan (or energy scan) is written
   by the beamline as M consecutive blocks of L frames into one a.his,
   with the same conv.bat and output.log for every block - exactly what
   his_spl_E splits into 001/raw .. 00M/raw.  When the file holds a whole
   multiple M > 1 of L frames, RHIS_SCAN = 1..M selects the block and
   every frame number of conv.bat is offset by (M-1)*L.  Such a file is
   refused without RHIS_SCAN: the first block is never used silently. */
static int	Tokenize(char *line,char **tok,int maxtok);	/* defined below */

/* the highest frame number conv.bat uses (= frames per scan), 0 if none */
static long	HisConvMaxFrame(HisSet *s)
{
	FILE	*f;
	char	line[LEN],*tok[LEN/2],path[LEN],c;
	int	ntok,i;
	long	n,L=0;

	if (snprintf(path,LEN,"%s/%s",s->dir,s->conv)>=LEN) HisError(s->conv,"path too long.");
	if ((f=fopen(path,"r"))==NULL) HisError(path,"file not open.");
	while (fgets(line,LEN,f)!=NULL) {
	    line[strcspn(line,"\r\n")]='\0';
	    if ((ntok=Tokenize(line,tok,LEN/2))<0) HisError(s->conv,"too many words on one line.");
	    for (i=1; i<ntok; i++)
		if (ParseImgName(tok[i],&c,&n) && c!='q' && n>L) L=n;
	}
	(void)fclose(f);
	return L;
}

/* RHIS_SCAN: pick the scan (frames_per_scan and nscan are already set) */
static void	HisScanLayout(HisSet *s)
{
	char	*env;
	long	L=s->frames_per_scan;

	if ((env=getenv("RHIS_SCAN"))!=NULL && *env!='\0') {
	    char	*end;
	    long	k=strtol(env,&end,10);

	    if (*end!='\0' || k<1) HisError(s->path,"RHIS_SCAN must be a scan number (1, 2, ...).");
	    if (L<=0) HisError(s->conv,"RHIS_SCAN is set but conv.bat names no frame.");
	    if (s->nframes%L!=0) {
		(void)fprintf(stderr,"rhis: %s : %d frames is not a whole number of scans of %ld frames "
			      "(the highest frame conv.bat uses); RHIS_SCAN cannot be applied.\n",
			      s->path,s->nframes,L);
		exit(1);
	    }
	    if (k>s->nscan) {
		(void)fprintf(stderr,"rhis: %s : RHIS_SCAN=%ld but the file holds %d scan%s of %ld frames.\n",
			      s->path,k,s->nscan,(s->nscan==1)?"":"s",L);
		exit(1);
	    }
	    s->scan=(int)k; s->frame0=(int)((k-1)*L);
	}
	else if (s->nscan>1) {
	    (void)fprintf(stderr,"rhis: %s : %d frames = %d scans of %ld frames (conv.bat uses frames "
			  "1..%ld).  This looks like a Z scan: set RHIS_SCAN=1..%d to pick the scan "
			  "(RHIS_SCAN=1 if it really is a single scan).\n",
			  s->path,s->nframes,s->nscan,L,L,s->nscan);
	    exit(1);
	}
}

/* q<N>.img -> N, or -1 */
static long	QIndex(const char *tok)
{
	char	c;
	long	n;

	return (ParseImgName(tok,&c,&n) && c=='q')?n:-1L;
}

typedef struct {
	int	cmd;		/* 0 skip, 1 his2img, 2 img_ave, 3 ren, 4 copy */
	int	ntok;
	char	*tok[LEN/2];
} ConvLine;

static int	Tokenize(char *line,char **tok,int maxtok)
{
	int	n=0;

	while (*line) {
	    while (*line && isspace((unsigned char)*line)) line++;
	    if (!*line) break;
	    if (n>=maxtok) return -1;
	    tok[n++]=line;
	    while (*line && !isspace((unsigned char)*line)) line++;
	    if (*line) *line++='\0';
	}
	return n;
}

static void	HisDefineQ(HisSet *s,long q,int lineno,HisQ *def)
{
	HisQ	*Q;

	if (q<0 || q>=s->Nq) {
	    (void)fprintf(stderr,"rhis: %s line %d : bad q number %ld.\n",s->conv,lineno,q);
	    exit(1);
	}
	Q=s->Q+q;
	if (Q->kind!=HQ_NONE) {
	    (void)fprintf(stderr,"rhis: %s line %d : q%04ld.img is defined twice.\n",
			  s->conv,lineno,q);
	    exit(1);
	}
	*Q=*def;
	if (Q->kind==HQ_AVG) {
	    Q->list=(int *)HisMalloc(sizeof(int)*(size_t)Q->nlist,"frame list");
	    memcpy(Q->list,def->list,sizeof(int)*(size_t)Q->nlist);
	    s->n_avg++;
	}
	else s->n_frame++;
	Q->cy1=Q->cy2=(-1); Q->cache=NULL;
}

static void	HisReadConv(HisSet *s)
{
	FILE	*f;
	char	line[LEN],line2[LEN],*tok[LEN/2],path[LEN];
	int	ntok,lineno=0,pass,*list=NULL,nlist,i;
	long	q,maxq=-1;

	if (snprintf(path,LEN,"%s/%s",s->dir,s->conv)>=LEN) HisError(s->conv,"path too long.");
	if ((f=fopen(path,"r"))==NULL) HisError(path,"file not open.");

	/* pass 0: find the highest q number (and the his file name);
	   pass 1: build the table */
	for (pass=0; pass<2; pass++) {
	    rewind(f); lineno=0;
	    while (fgets(line,LEN,f)!=NULL) {
		lineno++;
		line[strcspn(line,"\r\n")]='\0';
		strcpy(line2,line);
		if ((ntok=Tokenize(line2,tok,LEN/2))<0)
		    HisError(s->conv,"too many words on one line.");
		if (ntok==0) continue;

		if (StrCaseEq(tok[0],"rem") || !strncmp(tok[0],"::",2) ||
		    StrCaseEq(tok[0],"@echo") || StrCaseEq(tok[0],"echo") ||
		    StrCaseEq(tok[0],"pause"))
		    continue;
		if (StrCaseEq(tok[0],"hp2do")) {
		    /* not executed; its six image names (D000 I000 E000 D180
		       I180 E180) tell his_DO which views to compare */
		    if (pass==1 && s->n_hp2do==0)
			for (i=1; i<ntok && i<=6 && *tok[i]!='>'; i++) {
			    if ((s->hp2do[i-1]=strdup(tok[i]))==NULL) HisError(s->conv,"no memory.");
			    s->n_hp2do=i;
			}
		    continue;
		}

		if (StrCaseEq(tok[0],"his2img")) continue;	/* handled in HisOpen */
		if (StrCaseEq(tok[0],"img_ave")) {
		    if (ntok<3) {
			(void)fprintf(stderr,"rhis: %s line %d : img_ave needs at least "
				      "one source and one destination.\n",s->conv,lineno);
			exit(1);
		    }
		    if (StrCaseEq(tok[ntok-1],"dark.img")) q=-2; else q=QIndex(tok[ntok-1]);
		    if (q==-1) {
			(void)fprintf(stderr,"rhis: %s line %d : img_ave destination '%s' is "
				      "neither dark.img nor q<N>.img.\n",s->conv,lineno,tok[ntok-1]);
			exit(1);
		    }
		    if (pass==0) { if (q>maxq) maxq=q; continue; }

		    nlist=ntok-2;
		    list=(int *)HisMalloc(sizeof(int)*(size_t)nlist,"frame list");
		    for (i=0; i<nlist; i++) list[i]=FrameIndex(s,tok[1+i],lineno);
		    {
			HisQ	def; memset(&def,0,sizeof(def));
			def.kind=HQ_AVG; def.nlist=nlist; def.list=list;
			if (q==-2) {
			    if (s->dark.kind!=HQ_NONE) HisError(s->conv,"dark.img is defined twice.");
			    s->dark=def; s->dark.cy1=s->dark.cy2=(-1); s->dark.cache=NULL;
			    list=NULL;		/* owned by s->dark now */
			}
			else { HisDefineQ(s,q,lineno,&def); free(list); list=NULL; }
		    }
		    continue;
		}
		if (StrCaseEq(tok[0],"ren") || StrCaseEq(tok[0],"rename") ||
		    StrCaseEq(tok[0],"copy") || StrCaseEq(tok[0],"cp") ||
		    StrCaseEq(tok[0],"mv")) {
		    int	is_copy=StrCaseEq(tok[0],"copy") || StrCaseEq(tok[0],"cp");

		    if (ntok!=3) {
			(void)fprintf(stderr,"rhis: %s line %d : %s needs exactly two names.\n",
				      s->conv,lineno,tok[0]);
			exit(1);
		    }
		    if (StrCaseEq(tok[2],"dark.img")) q=-2; else q=QIndex(tok[2]);
		    if (q==-1) {
			(void)fprintf(stderr,"rhis: %s line %d : destination '%s' is neither "
				      "dark.img nor q<N>.img.\n",s->conv,lineno,tok[2]);
			exit(1);
		    }
		    if (pass==0) { if (q>maxq) maxq=q; continue; }
		    {
			HisQ	def; memset(&def,0,sizeof(def));
			long	src=QIndex(tok[1]);

			if (src>=0) {		/* q<M>.img as source: alias */
			    if (!is_copy) {
				(void)fprintf(stderr,"rhis: %s line %d : renaming a q file is "
					      "not supported.\n",s->conv,lineno);
				exit(1);
			    }
			    if (src>=s->Nq || s->Q[src].kind==HQ_NONE) {
				(void)fprintf(stderr,"rhis: %s line %d : source %s is not "
					      "defined yet.\n",s->conv,lineno,tok[1]);
				exit(1);
			    }
			    def=s->Q[src];
			}
			else if (StrCaseEq(tok[1],"dark.img")) {
			    if (s->dark.kind==HQ_NONE) HisError(s->conv,"dark.img used before it is defined.");
			    def=s->dark;
			}
			else {
			    def.kind=HQ_FRAME; def.frame=FrameIndex(s,tok[1],lineno);
			}
			if (q==-2) {
			    if (s->dark.kind!=HQ_NONE) HisError(s->conv,"dark.img is defined twice.");
			    s->dark=def;
			    if (def.kind==HQ_AVG) {
				s->dark.list=(int *)HisMalloc(sizeof(int)*(size_t)def.nlist,"frame list");
				memcpy(s->dark.list,def.list,sizeof(int)*(size_t)def.nlist);
			    }
			    s->dark.cy1=s->dark.cy2=(-1); s->dark.cache=NULL;
			}
			else HisDefineQ(s,q,lineno,&def);
		    }
		    continue;
		}
		(void)fprintf(stderr,"rhis: %s line %d : unknown command '%s' "
			      "(accepted: his2img img_ave ren copy hp2do rem echo).\n",
			      s->conv,lineno,tok[0]);
		exit(1);
	    }
	    if (pass==0) {
		if (maxq<0) HisError(s->conv,"no q<N>.img is produced by conv.bat.");
		s->Nq=(int)maxq+1;
		s->Q=(HisQ *)HisMalloc(sizeof(HisQ)*(size_t)s->Nq,"q table");
		memset(s->Q,0,sizeof(HisQ)*(size_t)s->Nq);
		for (q=0; q<s->Nq; q++) s->Q[q].cy1=s->Q[q].cy2=(-1);
	    }
	}
	(void)fclose(f);

	if (s->dark.kind==HQ_NONE) HisError(s->conv,"dark.img is not produced by conv.bat.");
}

/* ---------------------------------------------------------------------- */
/* public: open / close                                                   */
/* ---------------------------------------------------------------------- */

/* open the his file, read every frame header and the frame layout of
   conv.bat (frames per scan L, number of scans M); no RHIS_SCAN check */
static HisSet	*HisOpenHeaders(const char *dir)
{
	HisSet	*s=(HisSet *)HisMalloc(sizeof(HisSet),"HisSet");
	char	his_name[LEN]="a.his",path[LEN],*env;

	memset(s,0,sizeof(HisSet));
	s->fd=HIS_BAD_FD;
	s->dark.cy1=s->dark.cy2=(-1);

	if (dir==NULL || *dir=='\0') dir=".";
	s->dir=strdup(dir);
	if (s->dir==NULL) HisError(dir,"no memory.");

	if ((env=getenv("RHIS_C"))!=NULL && *env!='\0') s->conv=strdup(env);
	else s->conv=strdup("conv.bat");
	if (s->conv==NULL) HisError(dir,"no memory.");

	/* his file name: RHIS_H > the his2img line of conv.bat > "a.his".
	   The frame table must exist before conv.bat is interpreted (frame
	   numbers are validated against it), so the name is picked up in a
	   quick pre-scan of conv.bat here. */
	{
	    FILE	*f;
	    char	line[LEN],*tok[LEN/2];
	    int	ntok;

	    if ((env=getenv("RHIS_H"))!=NULL && *env!='\0')
		(void)snprintf(his_name,LEN,"%s",env);
	    else {
		if (snprintf(path,LEN,"%s/%s",dir,s->conv)>=LEN) HisError(s->conv,"path too long.");
		if ((f=fopen(path,"r"))==NULL) HisError(path,"file not open.");
		while (fgets(line,LEN,f)!=NULL) {
		    line[strcspn(line,"\r\n")]='\0';
		    if ((ntok=Tokenize(line,tok,LEN/2))>=2 && StrCaseEq(tok[0],"his2img")) {
			(void)snprintf(his_name,LEN,"%s",tok[1]); break;
		    }
		}
		(void)fclose(f);
	    }
	}

	if (snprintf(path,LEN,"%s/%s",dir,his_name)>=LEN) HisError(his_name,"path too long.");
	s->path=strdup(path);
	if (s->path==NULL) HisError(dir,"no memory.");

	if ((s->fd=HisFileOpen(s->path))==HIS_BAD_FD) HisError(s->path,"file not open.");

	HisScanHeaders(s);
	s->frames_per_scan=(int)HisConvMaxFrame(s);
	s->nscan=1; s->scan=1; s->frame0=0;
	if (s->frames_per_scan>0 && s->nframes>s->frames_per_scan &&
	    s->nframes%s->frames_per_scan==0)
	    s->nscan=s->nframes/s->frames_per_scan;
	return s;
}

HisSet	*HisOpen(const char *dir)
{
	HisSet	*s=HisOpenHeaders(dir);
	int	ndark;

	HisScanLayout(s);
	HisReadConv(s);

	ndark=(s->dark.kind==HQ_AVG)?s->dark.nlist:1;
	(void)fprintf(stderr,
		      "%s: %d x %d, %d frames; %s: dark = %d frame%s, q0..q%d "
		      "(%d averaged, %d single)",
		      s->path,s->Nx,s->Ny,s->nframes,s->conv,ndark,(ndark==1)?"":"s",
		      s->Nq-1,s->n_avg,s->n_frame);
	if (s->nscan>1)
	    (void)fprintf(stderr,"; scan %d / %d (frames %d-%d)",s->scan,s->nscan,
			  s->frame0+1,s->frame0+s->frames_per_scan);
	(void)fputc('\n',stderr);
	return s;
}

int	HisScan(const HisSet *s)	{ return s->scan; }
int	HisNscan(const HisSet *s)	{ return s->nscan; }

void	HisLayout(const char *dir,int *Nx,int *Ny,int *nframes,int *L,int *M)
{
	HisSet	*s=HisOpenHeaders(dir);

	*Nx=s->Nx; *Ny=s->Ny; *nframes=s->nframes; *L=s->frames_per_scan; *M=s->nscan;
	HisClose(s);
}

/* largest angle of a sample view (it == 1) in <dir>/output.log (or RHP_O);
   the his_* programs call a scan with 350 degrees or more an offset scan */
double	HisMaxSampleAngle(const char *dir)
{
	char	path[LEN],str[LEN],*env;
	FILE	*f;
	int	q,it,n=0;
	double	c,a,amax=-1.0e300;

	if ((env=getenv("RHP_O"))!=NULL && *env!='\0') (void)snprintf(path,LEN,"%s",env);
	else if (snprintf(path,LEN,"%s/output.log",dir)>=LEN) HisError(dir,"output.log : path too long.");
	if ((f=fopen(path,"r"))==NULL) HisError(path,"file not open.");
	while (fgets(str,LEN,f)!=NULL)
	    if (sscanf(str,"%d %lf %lf %d",&q,&c,&a,&it)==4) {
		if (it==1 && a>amax) amax=a;
		n++;
	    }
	(void)fclose(f);
	if (n==0) HisError(path,"no view found.");
	return amax;
}

/* "   % RHIS_SCAN=k" for the cmd-hst.log line of the his_* programs, or "" */
const char	*HisScanNote(void)
{
	static char	note[64];
	const char	*env=getenv("RHIS_SCAN");

	if (env==NULL || *env=='\0') return "";
	(void)snprintf(note,sizeof(note),"   %% RHIS_SCAN=%.20s",env);
	return note;
}

void	HisClose(HisSet *s)
{
	int	q;

	if (s==NULL) return;
	if (s->Q!=NULL) {
	    for (q=0; q<s->Nq; q++) { free(s->Q[q].list); free(s->Q[q].cache); }
	    free(s->Q);
	}
	free(s->dark.list); free(s->dark.cache);
	for (q=0; q<s->n_hp2do; q++) free(s->hp2do[q]);
	free(s->off);
	if (s->fd!=HIS_BAD_FD) HisFileClose(s->fd);
	free(s->dir); free(s->path); free(s->conv);
	free(s);
}

int	HisNx(const HisSet *s)		{ return s->Nx; }
int	HisNy(const HisSet *s)		{ return s->Ny; }
int	HisNframes(const HisSet *s)	{ return s->nframes; }
int	HisNq(const HisSet *s)		{ return s->Nq; }
const char *HisPath(const HisSet *s)	{ return s->path; }

int	HisHasQ(const HisSet *s,int q)
{
	return q>=0 && q<s->Nq && s->Q[q].kind!=HQ_NONE;
}

int	HisHp2doNames(const HisSet *s,const char **names)
{
	int	i;

	for (i=0; i<6; i++) names[i]=(i<s->n_hp2do)?s->hp2do[i]:NULL;
	return s->n_hp2do;
}

int	HisResolveName(const HisSet *s,const char *name)
{
	long	q;

	if (StrCaseEq(name,"dark.img")) return -1;
	q=QIndex(name);
	return (q>=0 && q<s->Nq && s->Q[q].kind!=HQ_NONE)?(int)q:-2;
}

/* ---------------------------------------------------------------------- */
/* public: band reads                                                     */
/* ---------------------------------------------------------------------- */

/* img_ave-compatible average of rows y1..y2 over Q->list:
   integer sum, division truncated.  Rows are split over the OpenMP
   threads (each thread owns its rows for every frame), so the result does
   not depend on the thread count. */
static void	HisAverageBand(HisSet *s,const HisQ *Q,int y1,int y2,WORD *dst)
{
	int	rows=y2-y1+1,Nx=s->Nx;

	#pragma omp parallel
	{
	    int	nt=1,id=0,r1,r2,k,i;
	    WORD	*tmp;
	    unsigned long long	*sum;
	    size_t	n;

#ifdef	_OPENMP
	    nt=omp_get_num_threads(); id=omp_get_thread_num();
#endif
	    r1=(int)((long long)rows*id/nt); r2=(int)((long long)rows*(id+1)/nt);
	    if (r2>r1) {
		n=(size_t)(r2-r1)*(size_t)Nx;
		tmp=(WORD *)HisMalloc(sizeof(WORD)*n,"average buffer");
		sum=(unsigned long long *)HisMalloc(sizeof(unsigned long long)*n,"average sum");
		memset(sum,0,sizeof(unsigned long long)*n);
		for (k=0; k<Q->nlist; k++) {
		    HisReadFrameBand(s,Q->list[k],y1+r1,y1+r2-1,tmp);
		    for (i=0; i<(int)n; i++) sum[i]+=tmp[i];
		}
		for (i=0; i<(int)n; i++)
		    dst[(size_t)r1*Nx+i]=(WORD)(sum[i]/(unsigned long long)Q->nlist);
		free(sum); free(tmp);
	    }
	}
}

static void	HisReadEntryBand(HisSet *s,HisQ *Q,int y1,int y2,WORD *dst)
{
	size_t	n;

	if (y1<0 || y2>=s->Ny || y1>y2) {
	    (void)fprintf(stderr,"rhis: bad row range %d-%d (frame height %d).\n",y1,y2,s->Ny);
	    exit(1);
	}
	if (Q->kind==HQ_FRAME) { HisReadFrameBand(s,Q->frame,y1,y2,dst); return; }

	n=(size_t)(y2-y1+1)*(size_t)s->Nx;
	/* averaged entries: serialize so the cache is consistent; the
	   averaging itself is OpenMP-parallel inside */
	#pragma omp critical (rhis_avg_cache)
	{
	    if (Q->cache!=NULL && Q->cy1==y1 && Q->cy2==y2)
		memcpy(dst,Q->cache,sizeof(WORD)*n);
	    else {
		HisAverageBand(s,Q,y1,y2,dst);
		free(Q->cache);
		Q->cache=(WORD *)HisMalloc(sizeof(WORD)*n,"average cache");
		memcpy(Q->cache,dst,sizeof(WORD)*n);
		Q->cy1=y1; Q->cy2=y2;
	    }
	}
}

void	HisReadDarkBand(HisSet *s,int y1,int y2,WORD *dst)
{
	HisReadEntryBand(s,&(s->dark),y1,y2,dst);
}

void	HisReadQBand(HisSet *s,int q,int y1,int y2,WORD *dst)
{
	if (!HisHasQ(s,q)) {
	    (void)fprintf(stderr,"rhis: q%04d.img is not defined by %s.\n",q,s->conv);
	    exit(1);
	}
	HisReadEntryBand(s,s->Q+q,y1,y2,dst);
}
