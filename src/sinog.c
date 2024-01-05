// program ct_sino
// 
// Required files are q???.img, dark.img.
// output file "s????.sin"

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "sif_f.h"


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
	short			comment_length;
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
unsigned short	data[MAXPIXL];

//dark file name and data
unsigned short	dark[MAXPIXL];
short	cent_flag;

// image profile from 'output.log'
float	shottime[MAX_SHOT], shotangle[MAX_SHOT];
int		Ishot[MAX_SHOT], NST, NI0, II0[MAX_SHOT];

// parameters for interpolation of incident beam
unsigned short	II01[MAXPIXL], II02[MAXPIXL];     // II01[x], II02[x]
unsigned short	I[MAXPIXL];                                // I[x]
double	I0[MAXPIXL];

// file head character  q???.img or q????.img
char	flhead[]   = "q";
char	darkfile[] = "dark.img";

// output file name
char	fn[15];

// flag for q001.img or q0001.img
int		iFlag;

// projection
short	N, n_total;

// calculation layer
long	s_layer, ndiv;

// total absorption
double	TA, TA1, TA2, TAM, TASD;

static FOM	**fom;


/*----------------------------------------------------------------------*/

static void Error(msg)
char        *msg;
{
	fputs(msg,stderr);
	fputc('\n',stderr);
	exit(1);
}


/*----------------------------------------------------------------------*/
int read_hipic(int nfq, unsigned short *data, Header *h, long ln)
{
	int		j, i_res;
	FILE	*fi;
	char	fname[20];

	if(nfq==0){
		sprintf(fname,"%s",darkfile);
	}
	else{
		if(iFlag==0) sprintf(fname, "%s%03d.img", flhead, nfq);
		if(iFlag==1) sprintf(fname, "%s%04d.img", flhead, nfq);
		if(iFlag==2) sprintf(fname, "%s%05d.img", flhead, nfq);
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
}

/*----------------------------------------------------------------------*/

int read_log()
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
	FILE		*fi;


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
			for (jx = 0; jx < N; ++jx){
				I0[jx]    = a[jx] * shottime[k] + b[jx];
				if ((I[jx]-dark[jx]) < 50){
					if(*(ilp+nshot)==0){
						printf("Warning \t");
//						printf("  jx = %d, I0 = %f, I = %d, dark = %d, ln =%d \n", jx, I0[jx], I[jx], dark[jx], ln);
//						printf("  II01 = %d, II02 = %d\n", II01[jx], II02[jx]);
//						printf("  t1   = %f, t2   = %f \n", t1, t2);
//						printf("  a = %f,  b = %f\n", a[jx], b[jx]);
						printf("  %d\t black\n", k);
						*(ilp+nshot)=1;
					}
				}
				*(po+N*nshot+jx)= log((double)(I0[jx])/(double)(I[jx]-dark[jx]));
			}
			if(*(ilp+nshot)==1){
				for(jx=0;jx<N;++jx){
					*(po+N*nshot+jx)=0.0;
				}
			}

			if(*(ilp+nshot)==0){
// roupe for correction
				p_sum = 0.0;
				p_ave = 0.0;
				for (jx=0; jx<10; ++jx){
					p_sum = p_sum + *(po+N*nshot+jx) + *(po+N*nshot+N-jx);
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
//	sprintf(fn,"s%04d.sin", ln);
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

	sprintf(fn,"s%05ld.tif", ln);
	StoreImageFile_Float(fn,N,NST,fom,SIF_F_desc);

		(void)printf("%s\t%s\n",fn,SIF_F_desc);

	return (0);
}

/*----------------------------------------------------------------------*/

float calc_center(){
int		i;
long	j, dx, dx0, N1;
double	mm1, mx1, mm2, mx2;
double	*p000, *p180, rsum, rd, rmsd0, rmsd;
float	center;
Header	h;

// file name and data for calculation of rotation center (data[y][x])
// odd = I_0, even = transmitted
char	c_file[13];
unsigned short	cent_1[MAXPIXL], cent_2[MAXPIXL];
unsigned short	cent_3[MAXPIXL], cent_4[MAXPIXL];

// read center file
//I_0 file
	if ((i = read_hipic(1, &cent_1[0], &h, s_layer)) != 0){
		printf("something wrong -- return value is %d", i);
		return(-1);
	}

// zero degree file(transmitted)
	if ((i = read_hipic(2, &cent_2[0], &h, s_layer)) != 0){
		printf("something wrong -- return value is %d", i);
		return(-1);
	}

//I_0 file
	if ((i = read_hipic(n_total-1, &cent_3[0], &h, s_layer)) != 0){
		printf("something wrong -- return value is %d", i);
		return(-1);
	}

//  180 degree file(transmitted)
	if ((i = read_hipic(n_total, &cent_4[0], &h, s_layer)) != 0){
		printf("something wrong -- return value is %d", i);
		return(-1);
	}

// calculate center with RMSD
		p000=(double *)malloc(N*sizeof(double));
		p180=(double *)malloc(N*sizeof(double));
		for(j=0;j<N;++j){
			*(p000+j)=log((double)(cent_2[j]-dark[j]) / (double)(cent_1[j]-dark[j]));
			*(p180+j)=log((double)(cent_4[j]-dark[j]) / (double)(cent_3[j]-dark[j]));
		}
		
		N1=N-1;
		rsum=0.0;
		for(j=0;j<N;++j){
			rd=*(p000+j)-*(p180+N1-j);
			rsum=rsum+rd*rd;
		}
		rmsd0=sqrt(rsum/(double)N);
		dx0=0;

		for(dx=1;dx<=N/2;++dx){
			rsum=0.0;
			for(j=0;j<N1-dx;++j){
				rd=*(p000+dx+j)-*(p180+N1-j);
				rsum=rsum+rd*rd;
			}
			rmsd=sqrt(rsum/(double)(N-dx));
			if(rmsd0>rmsd){
				rmsd0=rmsd;
				dx0=dx;
//				printf("%d\t%lf\t%lf\n", dx,rmsd0,rmsd);
			}

			rsum=0.0;
			for(j=0;j<N1-dx;++j){
				rd=*(p000+j)-*(p180-dx+N1-j);
				rsum=rsum+rd*rd;
			}
			rmsd=sqrt(rsum/(double)(N-dx));
			if(rmsd0>rmsd){
				rmsd0=rmsd;
				dx0=-dx;
//				printf("%d\t%lf\t%lf\n", -dx,rmsd0,rmsd);
			}
		}
		center=(N1+dx0)/2.0;
//		printf("%d\t%lf\t%lf\n", dx0,rmsd0,center);

	return(center);
}


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
	FILE		*fo;

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

//open output files
	if((fo = fopen("sino.tmp","w")) == NULL){ 
		printf("can not open %s for output\n", fn);
		return(1);
	}
	fprintf(fo, "%d\n", h.width);
	fprintf(fo, "%d\n", NST);
	fprintf(fo, "%04ld\n", s_layer);
	fprintf(fo, "%f\n", calc_center());
	fprintf(fo, "%lf\n", TAM); // Mean of Total Absorption from sinogram
	fprintf(fo, "%lf\n", TASD); // SD of Total Absorption from sinogram
	fclose(fo);

	return 0;
}
