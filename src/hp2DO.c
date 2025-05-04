
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "cell.h"
#include "sif.h"

static void	Abort(path,msg)
char		*path,*msg;
{
	fprintf(stderr,"%s : %s\n",path,msg); exit(1);
}

#define GetByte(file)	fgetc(file)

static int	GetWord(file)
FILE		*file;
{
	int	lo,hi;

	return ((lo=GetByte(file))==EOF ||
		(hi=GetByte(file))==EOF)?-1:(hi<<8)|lo;
}

static FILE	*OpenFile(path,Nx,Ny,firstFile)
char		*path;
int		*Nx,*Ny,firstFile;
{
	FILE	*file=fopen(path,"rb");
	int	len;

	if (file==NULL) Abort(path,"file not open.");

	if (GetByte(file)!='I' ||
	    GetByte(file)!='M') Abort(path,"bad magic number.");

	if ((len=GetWord(file))<0) Abort(path,"bad comment length.");

	if (firstFile) {
	    if ((*Nx=GetWord(file))<=0 ||
		(*Ny=GetWord(file))<=0) Abort(path,"bad image size.");
	}
	else
	    if (*Nx!=GetWord(file) ||
		*Ny!=GetWord(file)) Abort(path,"image size not match.");

	fseek(file,4L,1);		/* image offset */

	if (GetWord(file)!=2) Abort(path,"bad image type.");

	fseek(file,(long)(len+50),1);	/* reserved data and comment */

	return file;
}

static void	ReadLine(path,file,Nx,line)
char		*path;
FILE		*file;
int		Nx,*line;
{
	int	x,word;

	for (x=0; x<Nx; x++) {
	    if ((word=GetWord(file))<0) Abort(path,"unexpected end of file.");

	    line[x]=word;
	}
}

static double	Log(I0,I)
double		I0,I;
{
	return (I0<=0.0 || I<=0.0)?0.0:log(I0/I);
}

static void	ReadProjection(pathD,fileD,
			       pathI,fileI,
			       pathE,fileE,Nx,dark,work,p)
char		*pathD,*pathI,*pathE;
FILE		*fileD,*fileI,*fileE;
int		Nx,*dark,*work;
double		*p;
{
	int	x;

	ReadLine(pathD,fileD,Nx,dark);
	ReadLine(pathI,fileI,Nx,work);
	for (x=0; x<Nx; x++) p[x]=(double)work[x]-(double)dark[x];

	ReadLine(pathE,fileE,Nx,work);
	for (x=0; x<Nx; x++) p[x]=Log(p[x],(double)work[x]-(double)dark[x]);
}

static double	Sqrt(d)
double		d;
{
	return (d<0.0)?0.0:sqrt(d);
}

static double	CalculateRMS(N1,dx,x1,x2,p000,p180)
int		N1,dx,x1,x2;
double		*p000,*p180;
{
	int	x;
	double	d,sum=0.0;

	for (x=x1; x<x2; x++) {
	    d=p000[x+dx]-p180[N1-x]; sum+=(d*d);
	}
	return Sqrt(sum/(double)(x2-x1));
}

static void	Error(msg)
char		*msg;
{
	fputs(msg,stderr); fputc('\n',stderr); exit(1);
}

#define MA(cnt,ptr)	malloc((size_t)(cnt)*sizeof(*(ptr)))

#define PathD000	argv[1]
#define PathI000	argv[2]
#define PathE000	argv[3]
#define PathD180	argv[4]
#define PathI180	argv[5]
#define PathE180	argv[6]

#define BPS	(8*sizeof(Cell))

#define DESC_LEN	256

int	main(argc,argv)
int	argc;
char	**argv;
{
	FILE	*fileD000,*fileI000,*fileE000,
		*fileD180,*fileI180,*fileE180;
	int	Nx,Ny,z1,z2,*dark,*work,N1,N2,Nc,y,dx0,x1,x2,dx,x;
	double	*p000,*p180,**array,RMS0,RMS,R0,CC,D,A,B,RMS1,RMS2,dRMS,
		sCC=0.0,sCCy=0.0,sCCR0=0.0,sCCy2=0.0,sCCR02=0.0,sCCyR0=0.0;
	char	desc[DESC_LEN],
		*cmae="cell memory allocation error.";
	Cell	**cell;

	if (argc<7 || argc>10)
	Error("usage : hp2DO D000 I000 E000 D180 I180 E180 {y1 y2} {TIFF}");

	fileD000=OpenFile(PathD000,&Nx,&Ny,1);
	fileI000=OpenFile(PathI000,&Nx,&Ny,0);
	fileE000=OpenFile(PathE000,&Nx,&Ny,0);
	fileD180=OpenFile(PathD180,&Nx,&Ny,0);
	fileI180=OpenFile(PathI180,&Nx,&Ny,0);
	fileE180=OpenFile(PathE180,&Nx,&Ny,0);

	if (argc==7 || argc==8) {
	    z1=0; z2=Ny-1;
	}
	else
	    if ((z1=atoi(argv[7]))<0 ||
		(z2=atoi(argv[8]))>=Ny || z1>=z2) Error("bad slice range.");

	printf("%d\t%d\n",Nx,Ny);

	if ((dark=(int    *)MA(Nx,dark))==NULL ||
	    (work=(int    *)MA(Nx,work))==NULL ||
	    (p000=(double *)MA(Nx,p000))==NULL ||
	    (p180=(double *)MA(Nx,p180))==NULL)
	    Error("line memory allocation error.");

	N1=Nx-1; N2=Nx/2;

	if (argc==8 || argc==10) {
	    if ((array=(double **)MA(Ny,array))==NULL) Error(cmae);

	    Nc=N2*2+1;
	    for (y=0; y<Ny; y++)
		if ((array[y]=(double *)MA(Nc,*array))==NULL) Error(cmae);
	}

    /* Set number of threads */
    omp_set_num_threads(16);
	#pragma omp for
	for (y=0; y<Ny; y++) {
	    ReadProjection(PathD000,fileD000,
			   PathI000,fileI000,
			   PathE000,fileE000,Nx,dark,work,p000);

	    ReadProjection(PathD180,fileD180,
			   PathI180,fileI180,
			   PathE180,fileE180,Nx,dark,work,p180);

	    RMS0=CalculateRMS(N1,dx0=0,x1=0,x2=Nx,p000,p180);
	    if (argc==8 || argc==10) array[y][N2]=RMS0;

	    for (dx=1; dx<=N2; dx++) {
		if ((RMS=CalculateRMS(N1,-dx,dx,Nx,p000,p180))<RMS0) {
		    RMS0=RMS; dx0=(-dx); x1=dx; x2=Nx;
		}
		if (argc==8 || argc==10) array[y][N2+dx]=RMS;

		if ((RMS=CalculateRMS(N1,dx,0,Nx-dx,p000,p180))<RMS0) {
		    RMS0=RMS; dx0=dx; x1=0; x2=Nx-dx;
		}
		if (argc==8 || argc==10) array[y][N2-dx]=RMS;
	    }
	    R0=(-(double)(N1+dx0)/2.0);

	    CC=0.0; for (x=x1; x<x2; x++) CC+=(p000[x+dx0]*p180[N1-x]);

	    printf("%d\t%lg\t%lf\t%lf\n",y,R0,RMS0,CC);

	    if (y>=z1 && y<=z2) {
		sCC+=CC;
		sCCy+=(CC*(double)y);
		sCCR0+=(CC*R0);
		sCCy2+=(CC*(double)y*(double)y);
		sCCR02+=(CC*R0*R0);
		sCCyR0+=(CC*(double)y*R0);
	    }
	}
	D=sCCy2*sCC-sCCy*sCCy;
	A=(sCC*sCCyR0-sCCy*sCCR0)/D;
	B=(sCCy2*sCCR0-sCCy*sCCyR0)/D;
	printf("%lf\t%lf\t%lf",A,B,Sqrt((sCCR02-A*sCCyR0-B*sCCR0)/sCC));
	if (argc==9 || argc==10) printf("\t%d\t%d",z1,z2);
	putchar('\n');

	if (argc==8 || argc==10) {
	    RMS1=RMS2=RMS0;
	    for (y=0; y<Ny; y++) for (dx=0; dx<Nc; dx++) {
		if (array[y][dx]<RMS1) RMS1=array[y][dx];
		if (array[y][dx]>RMS2) RMS2=array[y][dx];
	    }
	    dRMS=(RMS2-RMS1)/((double)(1<<BPS)-1.0);

	    cell=(Cell **)array;
	    for (y=0; y<Ny; y++) for (dx=0; dx<Nc; dx++)
		cell[y][dx]=(Cell)floor((array[y][dx]-RMS1)/dRMS+0.5);

	    sprintf(desc,"%d\t%lf\t%lf\n",Nx,RMS1,RMS2);
	    StoreImageFile(argv[argc-1],Nc,Ny,BPS,cell,desc);
	}
	return 0;
}
