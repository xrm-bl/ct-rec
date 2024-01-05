
#include <stdio.h>
#include <stdlib.h>
#include <setjmp.h>

static jmp_buf	env;

static int	FGetC(file)
FILE		*file;
{
	int	code;

	if ((code=fgetc(file))==EOF) longjmp(env,1); return code;
}

static void	FSeek(file,offset,origin)
FILE		*file;
long		offset;
int		origin;
{
	if (fseek(file,offset,origin)) longjmp(env,2);
}

typedef unsigned short	Word2;
typedef unsigned long	Word4;

static Word2	Get2I(tiff)
FILE		*tiff;
{
	Word2	lo=(Word2)FGetC(tiff);

	return lo|((Word2)FGetC(tiff)<<8);
}

static Word4	Get4I(tiff)
FILE		*tiff;
{
	Word4	lo=(Word4)Get2I(tiff);

	return lo|((Word4)Get2I(tiff)<<16);
}

static Word2	Get2M(tiff)
FILE		*tiff;
{
	Word2	hi=(Word2)FGetC(tiff);

	return (hi<<8)|(Word2)FGetC(tiff);
}

static Word4	Get4M(tiff)
FILE		*tiff;
{
	Word4	hi=(Word4)Get2M(tiff);

	return (hi<<16)|(Word4)Get2M(tiff);
}

#define ERROR(msg)	fprintf(stderr,"%s : %s\n",path,msg)

static void	PrintImageDescription(path)
char		*path;
{
	FILE	*tiff;
	int	code;
	Word2	(*Get2)(),cnt;
	Word4	(*Get4)(),ifd,pos,len;

	if ((tiff=fopen(path,"rb"))==NULL) {
	    ERROR("file not found."); return;
	}
	if (setjmp(env)) {
	    fclose(tiff); ERROR("unexpected end of file."); return;
	}
	switch(code=FGetC(tiff)) {
	    case 'I' : Get2=Get2I; Get4=Get4I; break;
	    case 'M' : Get2=Get2M; Get4=Get4M; break;
	    default  : code=0; break;
	}
	if (!code || FGetC(tiff)!=code || (*Get2)(tiff)!=0x002a) {
	    fclose(tiff); ERROR("bad magic number."); return;
	}
	while (ifd=(*Get4)(tiff)) {
	    FSeek(tiff,ifd,0);
	    for (cnt=(*Get2)(tiff); cnt; cnt--)
		if ((*Get2)(tiff)==0x10e)
		    if ((*Get2)(tiff)==2) {
			pos=ftell(tiff);

			if ((len=(*Get4)(tiff))>4L) FSeek(tiff,(*Get4)(tiff),0);
			while (len-- && (code=FGetC(tiff))!='\0') putchar(code);
			putchar('\n');

			FSeek(tiff,pos+8L,0);
		    }
		    else
			FSeek(tiff,8L,1);
		else
		    FSeek(tiff,10L,1);
	}
	fclose(tiff);
}

int	main(argc,argv)
int	argc;
char	**argv;
{
	if (argc==1) {
	    fputs("usage : pid TIFF {TIFF ...}\n",stderr); exit(1);
	}
	while (--argc) PrintImageDescription(*(++argv)); return 0;
}
