
#include <stdio.h>
#include <stdlib.h>
#include "cell.h"
#include "rif.h"

#define BITS	(8*(int)sizeof(Cell))

typedef unsigned int	Word4;
typedef unsigned short	Word2;
typedef unsigned char	Byte;

static char	*Path;

static void	Error(msg)
char		*msg;
{
	(void)fprintf(stderr,"%s : %s\n",Path,msg); exit(1);
}

#define WARNING(msg)	(void)fprintf(stderr,"%s : %s (warning).\n",Path,msg)

static FILE	*file;

static int	GetChar()
{
	int	data=fgetc(file);

	if (data==EOF) Error("unexpected end of file.");

	return data;
}

static Word4	bytes;

static Byte	bitN[256]={
		  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
		 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
		 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
		 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
		 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
		 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
		 96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
		112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,
		128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
		144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
		160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
		176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
		192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
		208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
		224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
		240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255
		},
		bitR[256]={
		  0,128, 64,192, 32,160, 96,224, 16,144, 80,208, 48,176,112,240,
		  8,136, 72,200, 40,168,104,232, 24,152, 88,216, 56,184,120,248,
		  4,132, 68,196, 36,164,100,228, 20,148, 84,212, 52,180,116,244,
		 12,140, 76,204, 44,172,108,236, 28,156, 92,220, 60,188,124,252,
		  2,130, 66,194, 34,162, 98,226, 18,146, 82,210, 50,178,114,242,
		 10,138, 74,202, 42,170,106,234, 26,154, 90,218, 58,186,122,250,
		  6,134, 70,198, 38,166,102,230, 22,150, 86,214, 54,182,118,246,
		 14,142, 78,206, 46,174,110,238, 30,158, 94,222, 62,190,126,254,
		  1,129, 65,193, 33,161, 97,225, 17,145, 81,209, 49,177,113,241,
		  9,137, 73,201, 41,169,105,233, 25,153, 89,217, 57,185,121,249,
		  5,133, 69,197, 37,165,101,229, 21,149, 85,213, 53,181,117,245,
		 13,141, 77,205, 45,173,109,237, 29,157, 93,221, 61,189,125,253,
		  3,131, 67,195, 35,163, 99,227, 19,147, 83,211, 51,179,115,243,
		 11,139, 75,203, 43,171,107,235, 27,155, 91,219, 59,187,123,251,
		  7,135, 71,199, 39,167,103,231, 23,151, 87,215, 55,183,119,247,
		 15,143, 79,207, 47,175,111,239, 31,159, 95,223, 63,191,127,255
		},
		*byte;

#define BYTE()	byte[GetChar()]

static void	InitNICDecode(ofs,cnt)
Word4		ofs,cnt;
{
	(void)fseek(file,(long)ofs,0); bytes=cnt;
}

static int	NICDecode()
{
	if (bytes==0) Error("too few codes.");

	--bytes; return BYTE();
}

static void	TermNICDecode()
{
	if (bytes!=0) WARNING("too many codes");
}

static int	LZWnum,LZWlen;
static Word4	LZWbuf,LZWmask;

static Word2	LZWGetCode()
{
	for (LZWnum-=LZWlen; LZWnum<0; LZWnum+=8) {
	    if (bytes==0) Error("LZW too few codes.");

	    LZWbuf=(LZWbuf<<8)|BYTE(); --bytes;
	}
	return (Word2)((LZWbuf>>LZWnum)&LZWmask);
}

static int	LZWrem;
static Word2	LZWold,LZWcode;

static void	InitLZWDecode(ofs,cnt)
Word4		ofs,cnt;
{
	(void)fseek(file,(long)ofs,0); bytes=cnt;

	LZWnum=0; LZWlen=9; LZWmask=511; LZWrem=0; LZWold=256; LZWcode=258;
}

static int	LZWbase,LZWsize;
static Byte	LZWring[3839],LZWTail[3838],*LZWtail=LZWTail-258;
static Word2	LZWlast,LZWHead[3838],*LZWhead=LZWHead-258;

static int	LZWDecode()
{
	Word2	new;

	while (LZWrem==0) {
	    if ((new=LZWGetCode())==257) Error("LZW unexpected end of code.");

	    if (new==256) {
		LZWlen=9; LZWmask=511;
	    }
	    else if (LZWold==256) {
		if (new>=256) Error("LZW bad code after clear code.");

		LZWring[LZWbase=0]=(Byte)(LZWlast=new); LZWrem=LZWsize=1;
		LZWcode=258;
	    }
	    else {
		if (LZWcode==4096) Error("LZW missing clear code.");

		if (new>LZWcode) Error("LZW bad code.");

		if (new<LZWcode) {
		    LZWbase=LZWsize=0;
		    for (LZWlast=new; LZWlast>=258; LZWlast=LZWhead[LZWlast])
			LZWring[LZWsize++]=LZWtail[LZWlast];

		    LZWring[LZWsize++]=(Byte)LZWlast;
		}
		else {
		    LZWring[LZWbase=(LZWbase+3838)%3839]=(Byte)LZWlast;
		    ++LZWsize;
		}
		LZWrem=LZWsize;

		LZWhead[LZWcode]=LZWold; LZWtail[LZWcode]=(Byte)LZWlast;
		if (++LZWcode==LZWmask && LZWlen<12) {
		    ++LZWlen; LZWmask=(LZWmask<<1)|1;
		}
	    }
	    LZWold=new;
	}
	return (int)LZWring[(LZWbase+(--LZWrem))%3839];
}

static void	TermLZWDecode()
{
	if (LZWrem>0)
	    WARNING("LZW too many data");
	else {
	    for (LZWnum-=LZWlen; LZWnum<0 && bytes>0; LZWnum+=8) {
		LZWbuf=(LZWbuf<<8)|BYTE(); --bytes;
	    }
	    if (bytes!=0)
		WARNING("LZW too many code");
	    else
		if (!(LZWnum>=0    && ((LZWbuf>>LZWnum  )& LZWmask    )==257) &&
		    !(LZWnum>=(-1) && ((LZWbuf>>LZWnum+1)&(LZWmask>>1))==257))
		    WARNING("LZW missing end of code");
	}
}

static int	MPBGetCode()
{
	if (bytes==0) Error("MPB too few codes.");

	--bytes; return BYTE();
}

static int	MPBlen;

static void	InitMPBDecode(ofs,cnt)
Word4		ofs,cnt;
{
	(void)fseek(file,(long)ofs,0); bytes=cnt;

	MPBlen=0;
}

static int	MPBrep;

static int	MPBDecode()
{
	if (MPBlen==0) {
	    while ((MPBlen=MPBGetCode())==128) ;

	    if (MPBlen>128) {
		MPBlen-=256; return MPBrep=MPBGetCode();
	    }
	}
	else if (MPBlen<0) {
	    ++MPBlen; return MPBrep;
	}
	else
	    --MPBlen;

	return MPBGetCode();
}

static void	TermMPBDecode()
{
	if (MPBlen!=0)
	    WARNING("MPB too many data");
	else if (bytes!=0)
	    WARNING("MPB too many code");
}

static int	(*Get1)();

static Word2	Get2I()
{
	Word2	lo=(Word2)(*Get1)();

	return lo|((Word2)(*Get1)()<<8);
}

static Word2	Get2M()
{
	Word2	hi=(Word2)(*Get1)();

	return (hi<<8)|(Word2)(*Get1)();
}

static Word4	Get4I()
{
	Word4	lo=(Word4)Get2I();

	return lo|((Word4)Get2I()<<16);
}

static Word4	Get4M()
{
	Word4	hi=(Word4)Get2M();

	return (hi<<16)|Get2M();
}

static Word2	(*Get2)();
static Word4	(*Get4)();

static int	GetIFD(ifd,tag)
Word4		ifd;
int		tag;
{
	Word2	cnt;

	(void)fseek(file,(long)ifd,0);
	for (cnt=(*Get2)(); cnt!=0; cnt--) {
	    if ((*Get2)()==tag) return (int)(*Get2)();

	    (void)fseek(file,10L,1);
	}
	return 0;
}

static Word4	GetIFD4(ifd,tag,def)
Word4		ifd,def;
int		tag;
{
	switch(GetIFD(ifd,tag)) {
	    case 3 : if ((*Get4)()!=1) break; return (Word4)(*Get2)();
	    case 4 : if ((*Get4)()!=1) break; return (*Get4)();
	}
	return def;
}

static Word2	GetIFD3(ifd,tag,def)
Word4		ifd;
int		tag;
Word2		def;
{
	return (GetIFD(ifd,tag)==3 && (*Get4)()==1)?(*Get2)():def;
}

#define ALLOC(type,size)	(type *)malloc(sizeof(type)*(size_t)(size))

#define GETSTRIP_MSG_LEN	256
#define GETSTRIP_ERROR(format)	(void)sprintf(msg,format,title); Error(msg)

static Word4	*GetStrip(ifd,tag,spi,title)
Word4		ifd,spi;
int		tag;
char		*title;
{
	Word4	idx,*strip=ALLOC(Word4,spi);
	char	msg[GETSTRIP_MSG_LEN],*bns="bad number of strip %s.";

	if (strip==NULL) {
	    GETSTRIP_ERROR("no allocatable memory for strip %s.");
	}
	switch(GetIFD(ifd,tag)) {
	    case 3 :
		if ((*Get4)()!=spi) {
		    GETSTRIP_ERROR(bns);
		}
		if (spi>2) (void)fseek(file,(long)(*Get4)(),0);
		for (idx=0; idx<spi; idx++) strip[idx]=(*Get2)();
		break;
	    case 4 :
		if ((*Get4)()!=spi) {
		    GETSTRIP_ERROR(bns);
		}
		if (spi>1) (void)fseek(file,(long)(*Get4)(),0);
		for (idx=0; idx<spi; idx++) strip[idx]=(*Get4)();
		break;
	    default :
		GETSTRIP_ERROR("missing strip %s.");
	}
	return strip;
}

void	ReadImageFile(path,Nx,Ny,BPS,cell,desc)
char	*path,**desc;
int	*Nx,*Ny,*BPS;
#ifdef	FUNCTION_CELL
void	(*cell)();
#else
Cell	***cell;
#endif
{
	int	code,bps,(*Decode)(),dp,y,idx,x,bits;
	Word4	ifd,iw,il,len,rps,spi,*ofs,*cnt;
	void	(*InitDecode)(),(*TermDecode)();
	Cell	mask,*line,carry,data=0;
#ifndef	FUNCTION_CELL
	char	*cmae="no allocatable memory for cell.";
#endif

	Path=path;

	if ((file=fopen(path,"rb"))==NULL) Error("file not found.");

	switch(code=GetChar()) {
	    case 'I' : Get4=Get4I; Get2=Get2I; break;
	    case 'M' : Get4=Get4M; Get2=Get2M; break;
	    default  : code=0; break;
	}
	Get1=GetChar;

	if (code==0 || GetChar()!=code || (*Get2)()!=0x002a)
	    Error("bad magic number.");

	if ((ifd=(*Get4)())==0) Error("no image file directory.");

	if (Nx!=NULL || Ny!=NULL || cell!=NULL) {
	    if ((iw=GetIFD4(ifd,0x100,(Word4)0))==0 ||
		(il=GetIFD4(ifd,0x101,(Word4)0))==0)
		Error("missing image size.");

	    if ((iw>>8*sizeof(int)-1)!=0 ||
		(il>>8*sizeof(int)-1)!=0) Error("too large image size.");

	    if (Nx!=NULL) *Nx=(int)iw;
	    if (Ny!=NULL) *Ny=(int)il;
	}
	if (BPS!=NULL || cell!=NULL) {
	    if ((bps=(int)GetIFD3(ifd,0x102,(Word2)1))<=0 || bps>BITS)
		Error("bad bits per sample.");

	    if (BPS!=NULL) *BPS=bps;
	}
	if (desc!=NULL)
	    if (GetIFD(ifd,0x10e)!=2)
		*desc=NULL;
	    else {
		if ((len=(*Get4)())>4) (void)fseek(file,(long)(*Get4)(),0);

		if ((*desc=ALLOC(char,len+1))==NULL)
		    Error("no allocatable memory for image description.");

		if (fread(*desc,(size_t)1,(size_t)len,file)!=len)
		    Error("missing image description.");

		(*desc)[len]='\0';
	    }

	if (cell!=NULL) {
	    switch(GetIFD3(ifd,0x103,(Word2)1)) {
		case 1 :
		    InitDecode=InitNICDecode;
		    TermDecode=TermNICDecode; Decode=NICDecode; break;
		case 5 :
		    InitDecode=InitLZWDecode;
		    TermDecode=TermLZWDecode; Decode=LZWDecode; break;
		case 32773 :
		    InitDecode=InitMPBDecode;
		    TermDecode=TermMPBDecode; Decode=MPBDecode; break;
		default :
		    Error("unknown image compression method.");
	    }
	    switch(GetIFD3(ifd,0x106,(Word2)0xffff)) {
		case 0  :
		case 1  :
		case 3  :
		case 4  : break;
		default : Error("bad photometric interpretation.");
	    }
	    switch(GetIFD3(ifd,0x10a,(Word2)1)) {
		case 1  : byte=bitN; break;
		case 2  : byte=bitR; break;
		default : Error("unknown fill order.");
	    }
	    if (GetIFD3(ifd,0x115,(Word2)1)!=1) Error("bad samples per pixel.");

	    if ((rps=GetIFD4(ifd,0x116,il))>il)
		rps=il;
	    else if (rps==0)
		Error("bad rows per strip.");

	    if (GetIFD3(ifd,0x11c,(Word2)1)!=1)
		Error("bad planar configuration.");

	    switch(GetIFD3(ifd,0x13d,(Word2)1)) {
		case 1  : dp=0; break;
		case 2  : dp=1; break;
		default : Error("unknown differencing predictor.");
	    }

#ifdef	FUNCTION_CELL
	    if ((line=ALLOC(Cell,iw))==NULL)
		Error("no allocatable memory for scanline.");
#else
	    if ((*cell=ALLOC(Cell *,il))==NULL) Error(cmae);

	    for (y=0; y<il; y++)
		if (((*cell)[y]=ALLOC(Cell,iw))==NULL) Error(cmae);
#endif
	    spi=(il+rps-1)/rps;
	    ofs=GetStrip(ifd,0x111,spi,"offset");
	    cnt=GetStrip(ifd,0x117,spi,"bytes count");

	    mask=(Cell)((Cell)(~0)>>BITS-bps);

	    Get1=Decode;

	    for (y=0; y<il; y++) {
#ifndef	FUNCTION_CELL
		line=(*cell)[y];
#endif
		if (y%rps==0) {
		    idx=y/rps; (*InitDecode)(ofs[idx],cnt[idx]);
		}
		if (bps==8)
		    for (x=0; x<iw; x++) line[x]=(Cell)(*Get1)();
		else if (bps==16)
		    for (x=0; x<iw; x++) line[x]=(Cell)(*Get2)();
		else if (bps==32)
		    for (x=0; x<iw; x++) line[x]=(Cell)(*Get4)();
		else
		    for (bits=(-bps), x=0; x<iw; x++, bits-=bps) {
			carry=(Cell)((bits+BITS>=8)?0:data<<-bits);

			for (; bits<0; bits+=8)
			    data=(Cell)((data<<8)|(*Get1)());

			line[x]=(carry|(data>>bits))&mask;
		    }

		if (y+1==il || (y+1)%rps==0) (*TermDecode)();

		if (dp!=0)
		    for (x=1; x<iw; x++) line[x]=(line[x-1]+line[x])&mask;

#ifdef	FUNCTION_CELL
		for (x=0; x<iw; x++) (*cell)(x,y,line[x]);
#endif
	    }
	    free(cnt); free(ofs);
#ifdef	FUNCTION_CELL
	    free(line);
#endif
	}
	(void)fclose(file);
}
