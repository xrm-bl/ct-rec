/* rec2rec.c : re-projection / re-reconstruction of CT slices.
 *
 *   usage : rec2rec rec_in/ rec_out/ [Nt]
 *
 * Every rec*.tif (32bit float CT slice, N x N) found in rec_in/ is
 *   (1) forward Radon-transformed into an Nt-view sinogram (radon.h;
 *       rotation axis = image centre, 180 deg, unit-pixel line integrals),
 *   (2) reconstructed exactly the way sf_rec does it -- sort-filter ring
 *       removal followed by CBP -- and
 *   (3) stored under the SAME file name in rec_out/.
 *
 * No unit scaling is applied anywhere: the input values are already
 * normalised (e.g. cm^-1 from hp_tg/sf_rec), Radon and CBP are linear and
 * mutually inverse for unit pixel spacing, so the output reproduces the
 * input values directly.
 *
 * Nt defaults to the original scan's view count recorded in the input
 * TIFF's ImageDescription ("Dr RC Nt RA0 min max", written by hp_tg /
 * sf_rec / ofct_srec); if that cannot be parsed, Nt=N.  The optional 3rd
 * argument overrides both.  Dr/RA0 from the description are carried over
 * to the output description unchanged.
 *
 * CPU build: radon_omp.c + sort_filter_omp.c + cbp_thread (makefileCPU).
 * GPU build: radon_g.cu  + sort_filter_g.cu  + cbp.cu      (makefileGPU).
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#ifdef	_WIN32
#include "msdirent.h"
#include <direct.h>
#define MKDIR(p)	_mkdir(p)
#else
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#define MKDIR(p)	mkdir(p,0755)
#endif

#include "tiffio.h"
#include "tifwrite.h"
#include "cbp.h"
#include "radon.h"
#ifdef USE_GPU
  #include "sort_filter_g.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_gpu
#else
  #include "sort_filter_omp.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_omp
#endif

extern void	Error(char *);

#define LEN	2048

static int	Nx,Ny,BPS;
static char	desc[512];

/*----------------------------------------------------------------------*/

static void Store32TiffFile(char *wname,int wX,int wY,float *data32,char *wdesc)
{
	TIFF	*image;

	if ((image=TIFFOpen(wname,"w"))==NULL) {
	    fprintf(stderr,"cannot open %s for writing\n",wname); exit(1);
	}
	TIFFSetField(image, TIFFTAG_IMAGEWIDTH, wX);
	TIFFSetField(image, TIFFTAG_IMAGELENGTH, wY);
	TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(image, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
	TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, wdesc);
	TIFFSetField(image, TIFFTAG_ARTIST, "rec2rec");

	ct_write_raw_strips(image,data32,(uint32_t)wX,(uint32_t)wY,sizeof(float));

	TIFFClose(image);
}

/* read a 32bit float TIFF into data32 (iHead=1: header only);
   the ImageDescription is copied into desc[] (empty string if absent). */
static void Read32TiffFile(char *rname,int iHead,float *data32)
{
	TIFF	*image;
	long	i,j;
	float	*rline;
	char	*d=NULL;

	if ((image=TIFFOpen(rname,"r"))==NULL) {
	    fprintf(stderr,"cannot open %s\n",rname); exit(1);
	}
	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &Nx);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, &Ny);
	TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &BPS);
	desc[0]='\0';
	if (TIFFGetField(image, TIFFTAG_IMAGEDESCRIPTION, &d) && d!=NULL) {
	    strncpy(desc,d,sizeof(desc)-1); desc[sizeof(desc)-1]='\0';
	}
	if (iHead==1) { TIFFClose(image); return; }

	if ((rline=(float *)_TIFFmalloc(TIFFScanlineSize(image)))==NULL)
	    Error("cannot allocate memory for line scan.");
	for (i=0; i<Ny; i++) {
	    if (TIFFReadScanline(image,rline,i,0)<0) {
		fprintf(stderr,"cannot get tif line -> %ld\n",i); exit(1);
	    }
	    for (j=0; j<Nx; j++) data32[i*(long)Nx+j]=rline[j];
	}
	_TIFFfree(rline);
	TIFFClose(image);
}

/*----------------------------------------------------------------------*/
/* list rec*.tif in dir, sorted by name */

static int cmpstr(const void *a,const void *b)
{
	return strcmp(*(const char * const *)a,*(const char * const *)b);
}

static int ListRecTif(const char *dir,char ***listp)
{
	DIR		*Dir;
	struct dirent	*sd;
	char		**list=NULL;
	int		n=0,cap=0;
	size_t		l;

	if ((Dir=opendir(dir))==NULL) Error("input directory not open.");
	while ((sd=readdir(Dir))!=NULL) {
	    if (strncmp(sd->d_name,"rec",3)!=0) continue;
	    l=strlen(sd->d_name);
	    if (l<7 || strcmp(sd->d_name+l-4,".tif")!=0) continue;
	    if (n==cap) {
		cap=(cap==0)?256:cap*2;
		if ((list=(char **)realloc(list,sizeof(char *)*cap))==NULL)
		    Error("memory allocation error for file list.");
	    }
	    if ((list[n]=(char *)malloc(l+1))==NULL)
		Error("memory allocation error for file name.");
	    strcpy(list[n],sd->d_name); n++;
	}
	closedir(Dir);
	if (n>1) qsort(list,n,sizeof(char *),cmpstr);
	*listp=list; return n;
}

/*----------------------------------------------------------------------*/

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC	1000000
extern long clock();
#endif
#define CLOCK()		((double)clock()/(double)CLOCKS_PER_SEC)

int	main(int argc,char *argv[])
{
	char	**name,fi[LEN],fo[LEN],comm[300];
	float	*data32=NULL,*sino=NULL,*result_data=NULL;
	Float	**P=NULL,**F;
	double	Dr,RCin,RA0,RC,data_min,data_max,t0,t1,t2;
	int	nfile,fidx,NtArg=0,Nt,NtDesc,N=0,cbpN=0,cbpM=0;
	long	i,j;
	int	kernel_size,num_threads;

	if (argc!=3 && argc!=4)
	    Error("usage : rec2rec rec_in/ rec_out/ [Nt]");
	if (argc==4 && (NtArg=atoi(argv[3]))<2)
	    Error("bad number of views.");

	if ((nfile=ListRecTif(argv[1],&name))==0)
	    Error("no rec*.tif in the input directory.");
	(void)MKDIR(argv[2]);	/* ok if it already exists */

	kernel_size=get_kernel_size_from_env();
	num_threads=get_num_threads_from_env();

	fprintf(stderr,"rec2rec: %d file(s) %s -> %s\n",nfile,argv[1],argv[2]);

	for (fidx=0; fidx<nfile; fidx++) {
	    t0=CLOCK();
	    snprintf(fi,LEN,"%s/%s",argv[1],name[fidx]);
	    snprintf(fo,LEN,"%s/%s",argv[2],name[fidx]);

	    Read32TiffFile(fi,1,NULL);
	    if (Nx!=Ny || BPS!=32) {
		fprintf(stderr,"skip %s (not a square 32bit image: %dx%d/%dbit)\n",
			name[fidx],Nx,Ny,BPS);
		continue;
	    }

	    /* Dr / Nt / RA0 from the description ("Dr RC Nt RA0 min max") */
	    Dr=1.0; RCin=0.0; NtDesc=0; RA0=0.0;
	    (void)sscanf(desc,"%lf%lf%d%lf",&Dr,&RCin,&NtDesc,&RA0);
	    Nt=(NtArg>0)?NtArg:((NtDesc>1)?NtDesc:Nx);

	    if (N!=Nx) {	/* (re)allocate host buffers on size change */
		free(data32);
		N=Nx;
		if ((data32=(float *)malloc(sizeof(float)*(size_t)N*N))==NULL)
		    Error("memory allocation error for slice.");
	    }
	    if (cbpN!=N || cbpM!=Nt) {
		free(sino); free(result_data);
		if ((sino=(float *)malloc(sizeof(float)*(size_t)Nt*N))==NULL ||
		    (result_data=(float *)malloc(sizeof(float)*(size_t)Nt*N))==NULL)
		    Error("memory allocation error for sinogram.");
		if ((P=InitCBP(N,Nt))==NULL)	/* leaks on re-init; size changes are rare */
		    Error("memory allocation error.");
		cbpN=N; cbpM=Nt;
		fprintf(stderr,"%d\t%d\t%d\n",N,N,Nt);
	    }

	    Read32TiffFile(fi,0,data32);

	    /* (1) forward Radon: axis at the image centre, angles RA0+pi*m/Nt */
	    RC=(double)(N-1)/2.0;
	    RadonSlice(data32,N,sino,Nt,RA0);

	    /* (2) ring removal + CBP, exactly as in sf_rec */
	    if (SORT_FILTER_RESTORE(sino,result_data,N,Nt,kernel_size,num_threads)!=0)
		Error("sort filter image processing failed.");
	    for (j=0; j<Nt; j++)
		for (i=0; i<N; i++) P[j][i]=result_data[(size_t)j*N+i];

	    t1=CLOCK();
	    F=CBP(1.0,-RC,RA0);

	    /* (3) store under the same name; values unscaled (see header) */
	    data_max=-32000.; data_min=32000.;
	    for (j=0; j<N; j++)
		for (i=0; i<N; i++) {
		    float v=(float)F[j][i];
		    data32[(size_t)j*N+i]=v;
		    if (data_max<v) data_max=v;
		    if (data_min>v) data_min=v;
		}
	    snprintf(comm,sizeof(comm),"%f\t%f\t%d\t%f\t%f\t%f",
		     Dr,RC,Nt,RA0,(float)data_min,(float)data_max);
	    Store32TiffFile(fo,N,N,data32,comm);
	    t2=CLOCK();

	    fprintf(stderr,"stored: %s\t%d / %d\t%.2f\t%.2f\n",
		    name[fidx],fidx+1,nfile,t1-t0,t2-t1);
	}

	/* append the executed command to cmd-hst.log in the working directory */
	{
	    FILE	*ff;
	    int		ai;

	    if ((ff=fopen("cmd-hst.log","a"))!=NULL) {
		for (ai=0; ai<argc; ai++) fprintf(ff,"%s ",argv[ai]);
		fprintf(ff,"\n");
		fclose(ff);
	    }
	}
	return 0;
}
