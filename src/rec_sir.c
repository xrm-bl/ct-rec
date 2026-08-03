/*
 * rec_sir.c - 3D binning of a 32bit float CT volume (si_sir for rec images).
 *
 *   rec_sir in/ out/ B
 *
 *     in/  : directory holding the 32bit float slices rec?????.tif
 *     out/ : output directory (must exist before the run)
 *     B    : binning factor, applied to x, y and z alike (B >= 1)
 *
 * Every B x B x B block of voxels is replaced by its average.  This is the
 * si_sir operation specialised for the 32bit float reconstructions, so the
 * values stay float and no rounding to an integer takes place.
 *
 * Differences from si_sir, following the intended use here:
 *   - the file name is always rec%05d.tif, on input and on output.
 *   - the first output slice takes the number of the first input slice and
 *     the rest follow one by one, instead of si_sir's renumbering from 0.
 *     With B=3 and a volume rec00003..rec00011 the outputs are rec00003,
 *     rec00004 and rec00005.
 *   - blocks that would be incomplete are dropped (si_sir keeps them and
 *     divides by the number of voxels actually present).  The output size
 *     is therefore Nx/B by Ny/B by Nz/B, every block holding exactly
 *     B*B*B voxels.
 *
 * The ImageDescription tag ("Dr RC Nt RA0 min max", tab separated, written
 * by hp_tg / ofct_srec / sf_rec / rec2rec) is updated for the new grid:
 *   Dr -> Dr*B                    (the voxel is B times larger)
 *   RC -> (RC - (B-1)/2) / B      (same physical axis in the new grid;
 *                                  the (B-1)/2 term is the half-voxel shift
 *                                  of the block centre, so for B=1 nothing
 *                                  changes and for B=2 it is (RC-0.5)/2)
 *   Nt, RA0 are carried over, min and max are recomputed from the result.
 * If the description cannot be parsed it is copied through unchanged and a
 * warning is printed once.
 *
 * Build (Windows):
 *   cl /DWINDOWS /O2 /openmp /Ferec_sir.exe rec_sir.c libtiff.lib jpeg.lib lzma.lib zs.lib
 * Build (Linux):
 *   gcc -O2 -fopenmp -o rec_sir rec_sir.c -ltiff -lm
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include "tiffio.h"
#else
#include <tiffio.h>
#endif
#include "tifwrite.h"

#define LEN		1024
#define DESC_LEN	512
#define MAX_INDEX	99999		/* rec%05d.tif */

static void	Error(char *msg)
{
	fputs(msg,stderr);
	fputc('\n',stderr);
	exit(1);
}

static int	ExistFile(const char *path)
{
	FILE	*fp=fopen(path,"rb");

	if (fp == NULL) return 0;
	fclose(fp);
	return 1;
}

/* Drop trailing separators so that "rec" and "rec/" compare equal. */
static void	TrimSlash(char *s)
{
	size_t	n=strlen(s);

	while (n > 1 && (s[n-1] == '/' || s[n-1] == '\\')) s[--n]='\0';
}

/* The output numbering overlaps the input numbering, so writing into the
   input directory would destroy the slices while they are still being
   read.  This catches the plain "same string" case. */
static int	SamePath(const char *a,const char *b)
{
#ifdef _WIN32
	return _stricmp(a,b) == 0;
#else
	return strcmp(a,b) == 0;
#endif
}

/*----------------------------------------------------------------------*/
/* open a slice, check its geometry, hand back the description */

static TIFF	*OpenSlice(const char *path,int *Nx,int *Ny,char *desc,int dlen)
{
	TIFF		*image;
	uint32_t	w=0,h=0;
	uint16_t	bps=0,spp=1,fmt=SAMPLEFORMAT_IEEEFP;
	char		*d=NULL;

	if ((image=TIFFOpen(path,"r")) == NULL) {
	    fprintf(stderr,"cannot open %s\n",path);
	    exit(1);
	}
	TIFFGetField(image,TIFFTAG_IMAGEWIDTH,&w);
	TIFFGetField(image,TIFFTAG_IMAGELENGTH,&h);
	TIFFGetField(image,TIFFTAG_BITSPERSAMPLE,&bps);
	TIFFGetFieldDefaulted(image,TIFFTAG_SAMPLESPERPIXEL,&spp);
	TIFFGetFieldDefaulted(image,TIFFTAG_SAMPLEFORMAT,&fmt);

	if (bps != 32) {
	    fprintf(stderr,"%s : not a 32bit image (%u bit).\n",path,(unsigned)bps);
	    exit(1);
	}
	if (spp != 1) {		/* the scanline buffer below assumes one sample */
	    fprintf(stderr,"%s : %u samples per pixel, only grayscale is supported.\n",
		    path,(unsigned)spp);
	    exit(1);
	}
	if (fmt != SAMPLEFORMAT_IEEEFP)
	    fprintf(stderr,"%s : containing non-float pixel values (warning).\n",path);

	if (desc != NULL) {
	    desc[0]='\0';
	    if (TIFFGetField(image,TIFFTAG_IMAGEDESCRIPTION,&d) && d != NULL) {
		strncpy(desc,d,(size_t)dlen-1);
		desc[dlen-1]='\0';
	    }
	}
	*Nx=(int)w;
	*Ny=(int)h;
	return image;
}

static void	Store32TiffFile(char *wname,int wX,int wY,float *data32,char *wdesc)
{
	TIFF	*image;

	if ((image=TIFFOpen(wname,"w")) == NULL) {
	    fprintf(stderr,"cannot open %s for writing\n",wname);
	    exit(1);
	}
	TIFFSetField(image,TIFFTAG_IMAGEWIDTH,wX);
	TIFFSetField(image,TIFFTAG_IMAGELENGTH,wY);
	TIFFSetField(image,TIFFTAG_BITSPERSAMPLE,32);
	TIFFSetField(image,TIFFTAG_COMPRESSION,COMPRESSION_NONE);
	TIFFSetField(image,TIFFTAG_PHOTOMETRIC,PHOTOMETRIC_MINISBLACK);
	TIFFSetField(image,TIFFTAG_SAMPLEFORMAT,SAMPLEFORMAT_IEEEFP);
	TIFFSetField(image,TIFFTAG_SAMPLESPERPIXEL,1);
	TIFFSetField(image,TIFFTAG_PLANARCONFIG,PLANARCONFIG_CONTIG);
	TIFFSetField(image,TIFFTAG_IMAGEDESCRIPTION,wdesc);
	TIFFSetField(image,TIFFTAG_ARTIST,"rec_sir");

	ct_write_raw_strips(image,data32,(uint32_t)wX,(uint32_t)wY,sizeof(float));

	TIFFClose(image);
}

/*----------------------------------------------------------------------*/

int	main(int argc,char *argv[])
{
	char	in[LEN],out[LEN],desc[DESC_LEN],wdesc[DESC_LEN+64];
	int	B,i,x,y,ox,oy,dz,bz;
	int	Nx=0,Ny=0,Nz,oNx,oNy,oNz,l_sta=-1,l_dst=-1;
	int	warned=0,haveDesc;
	double	*sum,inv,Dr,RC,RA0;
	int	Nt;
	float	*line,*outbuf;
	FILE	*flog;

	if (argc != 4) {
	    fputs("usage : rec_sir in/ out/ B\n"
		  "        in/  : directory holding rec?????.tif (32bit float)\n"
		  "        out/ : output directory (create it beforehand)\n"
		  "        B    : binning factor for x, y and z\n",stderr);
	    exit(1);
	}
	snprintf(in ,sizeof(in) ,"%s",argv[1]);
	snprintf(out,sizeof(out),"%s",argv[2]);
	TrimSlash(in);
	TrimSlash(out);
	if ((B=atoi(argv[3])) <= 0) Error("bad block size.");
	if (SamePath(in,out))
	    Error("the output directory must differ from the input directory.");

	/* --- locate the slices ------------------------------------------ */
	for (i = 0; i <= MAX_INDEX; i++) {
	    char path[LEN];

	    snprintf(path,sizeof(path),"%s/rec%05d.tif",in,i);
	    if (!ExistFile(path)) continue;
	    if (l_sta < 0) l_sta=i;
	    l_dst=i;
	}
	if (l_sta < 0) {
	    fprintf(stderr,"no rec?????.tif found in %s\n",in);
	    exit(1);
	}
	Nz=l_dst-l_sta+1;

	/* every slice of a block must be present */
	for (i = l_sta; i <= l_dst; i++) {
	    char path[LEN];

	    snprintf(path,sizeof(path),"%s/rec%05d.tif",in,i);
	    if (!ExistFile(path)) {
		fprintf(stderr,"%s is missing (the slice numbers must be contiguous).\n",path);
		exit(1);
	    }
	}

	/* --- geometry ---------------------------------------------------- */
	{
	    char path[LEN];
	    TIFF *image;

	    snprintf(path,sizeof(path),"%s/rec%05d.tif",in,l_sta);
	    image=OpenSlice(path,&Nx,&Ny,desc,DESC_LEN);
	    TIFFClose(image);
	}
	oNx=Nx/B;
	oNy=Ny/B;
	oNz=Nz/B;
	if (oNx < 1 || oNy < 1 || oNz < 1)
	    Error("block size larger than the volume.");

	fprintf(stderr,"rec_sir: %s rec%05d..rec%05d  %dx%dx%d  B=%d  ->  %s  %dx%dx%d\n",
		in,l_sta,l_dst,Nx,Ny,Nz,B,out,oNx,oNy,oNz);
	if (Nx%B || Ny%B || Nz%B)
	    fprintf(stderr,"  (dropping %d column(s), %d row(s), %d slice(s) that do not fill a block)\n",
		    Nx%B,Ny%B,Nz%B);

	/* --- buffers ----------------------------------------------------- */
	if ((sum   =(double *)malloc(sizeof(double)*(size_t)oNx*oNy)) == NULL ||
	    (outbuf=(float  *)malloc(sizeof(float )*(size_t)oNx*oNy)) == NULL ||
	    (line  =(float  *)malloc(sizeof(float )*(size_t)Nx      )) == NULL)
	    Error("no allocatable memory.");

	inv=1.0/((double)B*(double)B*(double)B);

	/* --- one output slice per B input slices ------------------------- */
	for (bz = 0; bz < oNz; bz++) {
	    int	first=l_sta+bz*B;	/* first input slice of this block */
	    char path[LEN];
	    double dmin,dmax;

	    memset(sum,0,sizeof(double)*(size_t)oNx*oNy);

	    for (dz = 0; dz < B; dz++) {
		TIFF	*image;
		int	sNx,sNy;

		snprintf(path,sizeof(path),"%s/rec%05d.tif",in,first+dz);
		image=OpenSlice(path,&sNx,&sNy,(dz==0 && bz==0)?desc:NULL,DESC_LEN);
		if (sNx != Nx || sNy != Ny) {
		    fprintf(stderr,"%s : %dx%d, expected %dx%d.\n",path,sNx,sNy,Nx,Ny);
		    exit(1);
		}
		for (y = 0; y < oNy*B; y++) {
		    double *srow;

		    if (TIFFReadScanline(image,line,(uint32_t)y,0) < 0) {
			fprintf(stderr,"cannot get tif line -> %d of %s\n",y,path);
			exit(1);
		    }
		    srow=sum+(size_t)(y/B)*oNx;
		    for (ox = 0; ox < oNx; ox++) {
			const float *p=line+(size_t)ox*B;
			double s=0.0;

			for (x = 0; x < B; x++) s+=(double)p[x];
			srow[ox]+=s;
		    }
		}
		TIFFClose(image);
	    }

	    dmin= 1e30;
	    dmax=-1e30;
	    for (oy = 0; oy < oNy; oy++)
	    for (ox = 0; ox < oNx; ox++) {
		double v=sum[(size_t)oy*oNx+ox]*inv;

		outbuf[(size_t)oy*oNx+ox]=(float)v;
		if (v < dmin) dmin=v;
		if (v > dmax) dmax=v;
	    }

	    /* description for the new grid */
	    haveDesc=(desc[0] != '\0');
	    if (haveDesc &&
		sscanf(desc,"%lf %lf %d %lf",&Dr,&RC,&Nt,&RA0) == 4)
		snprintf(wdesc,sizeof(wdesc),"%f\t%f\t%d\t%f\t%f\t%f",
			 Dr*(double)B,
			 (RC-((double)B-1.0)/2.0)/(double)B,
			 Nt,RA0,(float)dmin,(float)dmax);
	    else {
		if (haveDesc && !warned) {
		    fprintf(stderr,"cannot parse the image description; copying it unchanged.\n");
		    warned=1;
		}
		snprintf(wdesc,sizeof(wdesc),"%s",desc);
	    }

	    /* output numbering starts at the first input number and then
	       runs consecutively (rec00003, rec00004, ...) */
	    snprintf(path,sizeof(path),"%s/rec%05d.tif",out,l_sta+bz);
	    Store32TiffFile(path,oNx,oNy,outbuf,wdesc);
	    fprintf(stderr,"stored: %s\t%d / %d\t%.4f\t%.4f\r",path,bz+1,oNz,dmin,dmax);
	}
	fputs("\nfinish.\n",stderr);

	free(line);
	free(outbuf);
	free(sum);

	/* append to log file */
	if ((flog=fopen("cmd-hst.log","a")) != NULL) {
	    for (i = 0; i < argc; ++i) fprintf(flog,"%s ",argv[i]);
	    fprintf(flog,"\t   %% binning %d  %dx%dx%d -> %dx%dx%d\n",
		    B,Nx,Ny,Nz,oNx,oNy,oNz);
	    fclose(flog);
	}
	return 0;
}
