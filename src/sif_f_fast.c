
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "sif_f.h"

#ifndef	FOF
#define FOF	float
#endif

FOM	SIF_F_min,SIF_F_max;
char	SIF_F_form[SIF_F_LEN]=SIF_F_FORM,
	SIF_F_desc[SIF_F_LEN]="";

#if	UINT_MAX==4294967295
typedef unsigned int	Word4;
#elif   ULONG_MAX==4294967295
typedef unsigned long	Word4;
#else
#error no suitable type for 32-bit integer.
#endif

#if	UINT_MAX==65535
typedef unsigned int	Word2;
#elif   USHRT_MAX==65535
typedef unsigned short	Word2;
#else
#error no suitable type for 16-bit integer.
#endif

static void	Put2(Word2 word2,FILE *file) { (void)fwrite(&word2,2,1,file); }
static void	Put4(Word4 word4,FILE *file) { (void)fwrite(&word4,4,1,file); }

void	StoreImageFile_Float(char *path,int Nx,int Ny,FOM **cell,char *desc)
{
	FILE	*F;
	Word2	M=0x002a;
	Word4	S,p1,p2,p,L,l;
	int	y,x;
	FOM	fom;
	FOF	fof;

	if (path[0]=='-' && path[1]=='\0')
	    F=stdout;
	else if ((F=fopen(path,"wb"))==NULL) {
	    (void)fprintf(stderr,"%s : file not stored.\n",path); exit(1);
	}
	(void)fwrite((*((unsigned char *)&M)==M)?"II":"MM",2,1,F);
	Put2(M,F);
	Put4(2+2+4+(S=(Word4)Nx*(Word4)Ny*(Word4)sizeof(FOF)),F);

	SIF_F_min=SIF_F_max=cell[0][0]; if (sizeof(FOF)==8) p1=p2=p=0;
	for (y=0; y<Ny; y++)
	for (x=0; x<Nx; x++) {
	    if ((fom=cell[y][x])<SIF_F_min) {
		SIF_F_min=fom; if (sizeof(FOF)==8) p1=p++;
	    }
	    else if (fom>SIF_F_max) {
		SIF_F_max=fom; if (sizeof(FOF)==8) p2=p++;
	    }
	    else if (sizeof(FOF)==8) ++p;

	    fof=fom; (void)fwrite(&fof,sizeof(FOF),1,F);
//		size_t ret = fwrite(*cell, sizeof(FOM), (size_t)(Ny*Nx), F);
	}
	if (desc==NULL) {
	    L=0; Put2(13,F);
	}
	else {
	    if (desc==SIF_F_desc)
		(void)sprintf(SIF_F_desc,SIF_F_form,(double)SIF_F_min,
						    (double)SIF_F_max);

	    L=(Word4)strlen(desc)+1; Put2(14,F);
	}

#define PUTD(tag,type,count)	Put2(tag,F); Put2(type,F); Put4(count,F)
#define PUTD3(tag,word2)	PUTD(tag,3,1); Put2(word2,F); Put2(0,F)
#define PUTD4(tag,word4)	PUTD(tag,4,1); Put4(word4,F)

	PUTD4(0x100,Nx);	/* 0: image width */
	PUTD4(0x101,Ny);	/* 1: image length */
	PUTD3(0x102,(Word2)sizeof(FOF)*8);
				/* 2: bits per sample */
	PUTD3(0x103,1);		/* 3: image compression */
	PUTD3(0x106,1);		/* 4: photometric interpretation */
	if (L!=0) {
	    PUTD(0x10e,2,L);	/* 5: image description */
	    if (L>4)
		Put4(8+S+2+12*14+4,F);
	    else {
		(void)fwrite(desc,1,L,F); for (l=L; l<4; l++) (void)fputc(0,F);
	    }
	}
	PUTD4(0x111,8);		/* 5/6: strip offset */
	PUTD3(0x115,1);		/* 6/7: samples per pixel */
	PUTD4(0x116,Ny);	/* 7/8: rows per strip */
	PUTD4(0x117,S);		/* 8/9: strip byte counts */
	PUTD3(0x11c,1);		/* 9/10: planar configuration */
	PUTD3(0x153,3);		/* 10/11: sample format */
	if (sizeof(FOF)==4) {
	    PUTD(0x154,11,1); fof=SIF_F_min; (void)fwrite(&fof,4,1,F);
				/* 11/12: minimum sample value */
	    PUTD(0x155,11,1); fof=SIF_F_max; (void)fwrite(&fof,4,1,F);
				/* 12/13: maximum sample value */
	}
	else {
	    PUTD(0x154,12,1); Put4(8+p1*(Word4)sizeof(FOF),F);
				/* 11/12: minimum sample value */
	    PUTD(0x155,12,1); Put4(8+p2*(Word4)sizeof(FOF),F);
				/* 12/13: maximum sample value */
	}
	Put4(0,F);

	if (L>4) (void)fwrite(desc,1,L,F);

	if (F!=stdout) (void)fclose(F);
}
