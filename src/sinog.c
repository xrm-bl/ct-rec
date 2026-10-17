// program ct_sino
// 
// Required files are q???.img, dark.img.
// output file "s????.sin"
//
// Built with -DHIS_INPUT (his_sg) the same program reads a.his directly:
// conv.bat is interpreted in memory by rhis.c, output.log is read as
// usual, and no q????.img file is created.  Run it where a.his, conv.bat
// and output.log are.

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "blacklim.h"
#include "sif_f.h"

#ifdef HIS_INPUT
#include "rhis.h"
static HisSet	*his = NULL;	/* a.his + conv.bat in the current directory */
/* his_sg links this file together with of_sinog.c; his_main.c looks at the
   angles in output.log and calls sinog_main (180-degree scan).  File-scope names
   are static so that the two bodies can share one executable. */
#define main sinog_main
#endif


/*----------------------------------------------------------------------*/
#ifndef M_PI
#define M_PI				3.1415926535897932385
#endif
#define INTEL
#define HiPic_Header_Size	64
#define MAX_SHOT			20010
#define MAXPIXL				20001

/*----------------------------------------------------------------------*/
struct HiPic_Header{
	char			head[2];
	unsigned short			comment_length;
	short			width;
	short			height;
	short			x_offset;
	short			y_offset;
	short			type;
	char			reserved[63 - 13];
	char			*comment;
};
typedef struct HiPic_Header Header;

/*----------------------------------------------------------------------*/

// main data for read transmitted images (data[y][x]) 772
static unsigned short	data[MAXPIXL];

//dark file name and data
static unsigned short	dark[MAXPIXL];
static short	cent_flag;

// image profile from 'output.log'
static float	shottime[MAX_SHOT], shotangle[MAX_SHOT];
static int		Ishot[MAX_SHOT], NST, NI0, II0[MAX_SHOT];

// parameters for interpolation of incident beam
static unsigned short	II01[MAXPIXL], II02[MAXPIXL];     // II01[x], II02[x]
static unsigned short	I[MAXPIXL];                                // I[x]
static double	I0[MAXPIXL];

// file head character  q???.img or q????.img
static char	flhead[]   = "q";
static char	darkfile[] = "dark.img";

// output file name
static char	fn[256];

// flag for q001.img or q0001.img
static int		iFlag;

// projection
static short	N, n_total;

// calculation layer
static long	s_layer, ndiv;

// total absorption
static double	TA, TA1, TA2, TAM, TASD;

static FOM	**fom;


/*----------------------------------------------------------------------*/

static void Error(msg)
char        *msg;
{
	fputs(msg,stderr);
	fputc('\n',stderr);
	exit(1);
}


#ifdef HIS_INPUT
/* his_sg: one layer line of dark (nfq==0) or q<nfq>, straight from a.his.
   Averaged entries (dark, I0 groups) are averaged like img_ave inside
   rhis.c and cached there, so the repeated reads cost nothing.          */
static int read_his_line(int nfq, unsigned short *data, Header *h, long ln)
{
	if(his == NULL) his = HisOpen(".");
	if(ln < 0 || ln >= HisNy(his)){
		fprintf(stderr, "layer %ld is outside the frame (0 .. %d)\n", ln, HisNy(his)-1);
		return((int)ln+1);
	}
	if(nfq == 0) HisReadDarkBand(his, (int)ln, (int)ln, data);
	else{
		if(!HisHasQ(his, nfq)){
			fprintf(stderr, "q%04d.img is not produced by conv.bat\n", nfq);
			return(-10);
		}
		HisReadQBand(his, nfq, (int)ln, (int)ln, data);
	}
	N = HisNx(his);
	h->width = (short)HisNx(his);
	h->height = (short)HisNy(his);
	return(0);
}
#endif

/*----------------------------------------------------------------------*/
static int read_hipic(int nfq, unsigned short *data, Header *h, long ln)
{
#ifdef HIS_INPUT
	return read_his_line(nfq, data, h, ln);
#else
	int		j, i_res;
	FILE	*fi;
	char	fname[1024];

	if(nfq==0){
		snprintf(fname, sizeof(fname),"%s",darkfile);
	}
	else{
		if(iFlag==0) snprintf(fname, sizeof(fname), "%s%03d.img", flhead, nfq);
		if(iFlag==1) snprintf(fname, sizeof(fname), "%s%04d.img", flhead, nfq);
		if(iFlag==2) snprintf(fname, sizeof(fname), "%s%05d.img", flhead, nfq);
	}
	
//open input files
	if((fi = fopen(fname,"rb")) == NULL){
		printf("can not open %s for input\n", fname);
		free(h->comment); 
		return(-10);
	}
// read comment and image data from file1
	if (fread(h, sizeof(char), HiPic_Header_Size, fi) != HiPic_Header_Size){
		fprintf(stderr, "EOF in %s (header)\n", fname);
		free(h->comment); 
		return(-1);
	}
	if (strncmp(h->head, "IM", 2)){
		fclose(fi); 
		free(h->comment); 
		fprintf(stderr, "This is not HiPIc image file.\n"); 
		fprintf(stderr, "File name = %s\n", fname); 
		return(-11);
	}

	h->comment = (char *) malloc(h->comment_length + 1);
	if (h->comment == NULL){
		fclose(fi); 
		free(h->comment); 
		fprintf(stderr, "Connot allocate memory (HiPic_Header.comment)\n"); 
		fprintf(stderr, "File name = %s\n", fname); 
		return(-12);
	}
	if (fread(h->comment, sizeof(char), h->comment_length, fi) != (unsigned short)h->comment_length){ 
		fclose(fi); 
		free(h->comment); 
		fprintf(stderr, "EOF found (comment)\n"); 
		fprintf(stderr, "File name = %s\n", fname); 
		return(-13);
	}
//	x_size = h->width;
//	y_size = h->height;
	h->comment[h->comment_length] = '\0';
	i_res = fseek( fi, h->width*ln*sizeof(unsigned short), SEEK_CUR);
	if( i_res != 0 ){
		fclose(fi);
		free(h->comment); 
		fprintf(stderr, " Error seeking %s.\n", fname);
		fprintf(stderr, "File name = %s\n", fname); 
		return(ln+1);
	}
	else{
		if ((j = fread(data, sizeof(unsigned short), h->width, fi)) != h->width){
			fclose(fi);
			free(h->comment); 
			fprintf(stderr, " Error reading %s at %ld (%d)\n", fname, ln, j);
			fprintf(stderr, "File name = %s\n", fname); 
			return(ln+1);
		}
	}
//	*(data+772) = (*(data+771)+*(data+773))/2;

// check for correct reading of HiPic file
//printf("-------------------------------\n");
//printf("%s\n", fname);
//printf("%s\n", h->head);
//printf("%d\n", h->comment_length);
//printf("%d\n", h->width);
//printf("%d\n", h->height);
//printf("%d\n", h->x_offset);
//printf("%d\n", h->y_offset);
//printf("%d\n", h->type);
//printf("%s\n", h->reserved);
//printf("%s\n", h->comment);
//printf("-------------------------------\n");
	N = h->width;
	free(h->comment);
	fclose(fi);
	return(0);
#endif
}

/*----------------------------------------------------------------------*/

static int read_log()
{
	int			i, j;
	short		nnn;
	FILE		*f;
	char		lne[100], *ss;
	int			flg_I0;
// open parameter file
	if((f = fopen("output.log","r")) == NULL){
		printf("cannot open output.log \n");
		return(-10);
	}

// read parameters
	NI0 = 0;
	NST = 0;
	nnn = 0;
	j   = 0;
	while(j>=0){
		j=j+1;
		lne[0]='\0';
		ss = fgets(lne, 100, f);
		if (i = strlen(lne) > 2){
			sscanf(lne, "%hd %f %f %d", &nnn, &shottime[j], &shotangle[j], &flg_I0);
			if (flg_I0 == 0){
				II0[NI0] = nnn;
				NI0=NI0+1;				// number of I_0
			}
			else{
				Ishot[NST] = nnn;
				NST=NST+1;				// number of I
			}
		}
		else{
			break;
		}
	}
	n_total = (short)(nnn + 2);
	fprintf(stderr, " nshot = %d, NI0 = %d, total = %d \n", NST, NI0, n_total);
	fclose(f);
	if(n_total<1000) iFlag = 0;
	if(n_total>=1000) iFlag = 1;
	if(n_total>=10000) iFlag = 2;
	if(n_total>MAX_SHOT){
		fprintf(stderr, "Too many projections!");
		return(1);
	}
	return(0);
}

/*----------------------------------------------------------------------*/

static int StoreProjection(ln)
long		ln;
{
	int			i, j, k, jx, nshot, *ilp, iplc,x,y;
	long		ll;
	double		t1, t2;
	double		a[MAXPIXL], b[MAXPIXL];
	Header		h;
	char		ffi01[15], ffi02[15], fii[15];
	double		I01, I02;

	double		*po, *poa;
	double		p_sum,p_ave;
	double		black_thresh, black_frac;
	int			nclip, nclip_max;
	FILE		*fi;

	black_thresh = BlackThresh();
	black_frac   = BlackFrac();
	nclip_max    = (int)(black_frac * (double)N);

// p initialization
	po = (double *)malloc(N*NST*sizeof(double));
	TA1 =0.0;
	TA2 =0.0;

	ilp=(int *)malloc(NST*sizeof(int));
	iplc=0;

	poa = (double *)malloc(N*sizeof(double));
	for(j=0;j<N;++j) *(poa+j)=0.0;

// counting for number of projections (NST = nshot)
	nshot = 0;

// loop between I0_1st and I0_2nd
	for ( j = 0; j < NI0-1; ++j){

//IIO[j] and IIO[j+1] are opened
		if ((i = read_hipic(II0[j], &II01[0], &h, ln)) != 0){
			printf("something wrong -- return value is %d(II01)", i);
			return(-1);
		}
		if ((i = read_hipic(II0[j+1], &II02[0], &h, ln)) != 0){
			printf("something wrong -- return value is %d(II02)", i);
			return(-1);
		}
// I0EV
		t1 = shottime[II0[j]];
		t2 = shottime[II0[j+1]];

// 1 layer(ln)
		for (jx = 0; jx < N; ++jx){
			I01   = (double)II01[jx]-dark[jx];
			I02   = (double)II02[jx]-dark[jx];
			a[jx] = (double)(((double)(I02   - I01))    / (t2 - t1));
			b[jx] = (double)(((double)I01*t2 - (double)I02*t1) / (t2 - t1));
		}
		for ( k = II0[j] + 1; k < II0[j+1]; ++k){
			// obtain p(x) from a[jx], b[jx] using shottime[k]
			if ((i = read_hipic(k, &I[0], &h, ln)) != 0){
				printf("something wrong -- return value is %d(II01)", i);
				return(-1);
			}
			*(ilp+nshot)=0;
			/* Per-pixel floor (CT_REC_BLACK_THRESH).  This used to be a
			   hardwired "any pixel below 50 counts discards the whole
			   projection", which throws away every projection of a dense
			   sample and still let log(x/0) through.  Now the pixel is
			   clamped and only a mostly-clipped line is discarded. */
			nclip = 0;
			for (jx = 0; jx < N; ++jx){
				I0[jx]    = a[jx] * shottime[k] + b[jx];
				*(po+N*nshot+jx)= BlackLog((double)(I0[jx]),
					(double)(I[jx]-dark[jx]), black_thresh, &nclip);
			}
			if (nclip > nclip_max){
				printf("Warning \t");
				printf("  %d\t black \t clipped=%d/%d (%.1f%%, thresh=%.4g)\n",
					k, nclip, (int)N, 100.0*(double)nclip/(double)N,
					black_thresh);
				*(ilp+nshot)=1;
				BlackCountProjection();
				for(jx=0;jx<N;++jx){
					*(po+N*nshot+jx)=0.0;
				}
			}

			if(*(ilp+nshot)==0){
// roupe for correction
				p_sum = 0.0;
				p_ave = 0.0;
				/* air reference: outermost 10 px at each end of the line.
				   N-jx would reach po[N*nshot+N], the first pixel of the
				   next (still unwritten) projection - and one element
				   past the buffer for the last one. */
				for (jx=0; jx<10; ++jx){
					p_sum = p_sum + *(po+N*nshot+jx) + *(po+N*nshot+N-1-jx);
				}
				p_ave = p_sum / (10.-0.) / 2.;
				TA=0.0;
				for (jx=0; jx<N; ++jx){
//					*(po+N*nshot+jx) = *(po+N*nshot+jx) - p_ave; /* comment out here */
					TA = TA + *(po+N*nshot+jx);
					*(poa+jx)=*(poa+jx)+*(po+N*nshot+jx);
				}
				TA1=TA1+TA;
				TA2=TA2+TA*TA;
				iplc=iplc+1;
			}
//
			nshot = nshot + 1;
			if((nshot%100)==0) printf("%d\r",nshot);
		} // end of k loop
	} // end of j loop


// correcion for black projection
	for(k=0;k<NST;++k){
		if(*(ilp+k)==1){
//			printf("%d\t%d\n",k,*(ilp+k) );
			for(j=0;j<N;++j){
//				*(po+N*k+j)=*(poa+j)/(double)iplc;
//		(void)printf("%d\t%d\t%lf\n",j,k,*(po+N*k+j));
			}
		}
	}

	
//open output files
//	snprintf(fn, sizeof(fn),"s%04d.sin", ln);
//	if((fi = fopen(fn,"wb")) == NULL){ 
//		printf("failed to open %s for output\n", fn);
//		return(-10);
//	}
//	if(i=fwrite(po,sizeof(double),N*NST,fi)>N*NST*sizeof(double)){
//		fprintf(stderr, " Error writing data to %s at %d\n", fn, i);
//		fclose(fi);
//		return(i+1);
//	}
//
//	fclose(fi);
//		(void)printf("\n%d\n",ln);

	NST=NST/ndiv;
	
	if ((fom=(FOM **)malloc(sizeof(FOM *)*NST))==NULL ||
	    (*fom=(FOM *)malloc(sizeof(FOM)*N*NST))==NULL)
	    Error("memory allocation error for sinograms or tomograms.");
	for (y=1; y<NST; y++) fom[y]=fom[y-1]+N;

//	(void)printf("\n%d\n",ln);
	for (y=0; y<NST; y++){
		for (x=0; x<N; x++){
			fom[y][x]=(float)*(po+N*(y*ndiv)+x);
//		(void)printf("%d\t%d\t%lf\t%lf\n",x,y,fom[y][x],*(po+N*y+x));
		}
	}
//		(void)printf("\n%d\n",ln);

	BlackReport();

	snprintf(fn, sizeof(fn),"s%05ld.tif", ln);
	StoreImageFile_Float(fn,N,NST,fom,SIF_F_desc);

		(void)printf("%s\t%s\n",fn,SIF_F_desc);

	return (0);
}

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC	1000000
extern long clock();
#endif

#define CLOCK()		((double)clock()/(double)CLOCKS_PER_SEC)
#define DESC_LEN	256

int main(argc,argv)
int		argc;
char	**argv;
{
	double		Clock;					// timer setting
	int			i;
	Header		h;

// parameter setting
	ndiv=1;
	if (argc==3){
		ndiv=atoi(argv[2]);
	}
	else if (argc==2){
	}
	else{
		fprintf(stderr, "usage : %s layer {ndiv}\n", argv[0]);
		return(1);
	}
	s_layer = atoi(argv[1]);

// read shot log file
	if ((i = read_log()) != 0){
		fprintf(stderr, "canot read 'output.log' ! (%d)\n", i);
		return(1);
	}

// read dark image
	if ((i = read_hipic(0, &dark[0], &h, s_layer)) != 0){
		fprintf(stderr, "something wrong in dark file (%d)\n", i);
		return(1);
	}
//	fprintf(stderr, "width of projection data = %d \n",  h.width);

	cent_flag = 0;
	Clock=CLOCK();
	if((i=StoreProjection(s_layer)) !=0){
		printf("something wrong in StoreProjection (%d)\n",i);
		return(1);
	}
	
	NST=NST*ndiv;
	TAM=TA1/NST;
	TASD=sqrt(TA2/NST-TAM*TAM);
	fprintf(stderr," TA mean\t%lf\tTASD\t%lf\n", TAM, TASD);
	fprintf(stderr," TASD/TA mean\t%lf\n", TASD/TAM);
//	fprintf(stderr," Store Sinogram \t%lf / sec\n",CLOCK()-Clock);
//	fprintf(stderr, "output file name = %s \n", fn);

// (sino.tmp - width, projections, layer, centre, TA mean and SD for an
//  old batch chain - is no longer written: 2026-10-18)

// append to log file
	FILE		*ff;
	if((ff = fopen("../cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(ff,"%s ",argv[i]);
#ifdef HIS_INPUT
	{ extern const char *HisScanNote(void); fprintf(ff,"%s",HisScanNote()); }
#endif
	fprintf(ff,"\n");
	fclose(ff);

#ifdef HIS_INPUT
	HisClose(his);
#endif
	return 0;
}
