// program ct_rec
// 
// Required files are q???.img, dark.img.
// output file "rec????.tif"

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cbp.h"
#include <math.h>
#include "tiffio.h"
#include "sort_filter_omp.h"

/*----------------------------------------------------------------------*/
#ifndef M_PI
#define M_PI				3.1415926535897932385
#endif
#define INTEL
#define MAX_SHOT			10010
#define MAXPIXL				18001

/*----------------------------------------------------------------------*/

// main data for read transmitted images (data[y][x])
unsigned short	data[MAXPIXL];

//dark file name and data
unsigned short	dark[MAXPIXL];
unsigned short	cent_flag;

// image profile from 'output.log'
float	shottime[MAX_SHOT], shotangle[MAX_SHOT];
int		Ishot[MAX_SHOT], NI0, II0[MAX_SHOT];

// parameters for interpolation of incident beam
unsigned short	II01[MAXPIXL], II02[MAXPIXL];     // II01[x], II02[x]
unsigned short	I[MAXPIXL];                                // I[x]
double	I0[MAXPIXL];

double	*po;

// file head character  q???.img or q????.img
char	flhead[]   = "q";
char	darkfile[] = "dark.tif";

// output file name
char	fn[15];

// flag for q001.img or q0001.img
int		iFlag, nfq, BPS;

// projection
unsigned short	N, M, n_total;

// calculation layer
long	s_layer;

// total absorption
double	TA, TA1, TA2, TAM, TASD;


/*----------------------------------------------------------------------*/

int Read16TiffFileLine(char* rname, int iHead)
{
	TIFF* image;
	long i, j;
	unsigned short spp, com, pm, pc, rps;
	unsigned short* rline;
	
	char	*desc;
	int		Nx,Ny;

	image = TIFFOpen(rname, "r");

	TIFFGetField(image, TIFFTAG_IMAGEWIDTH, &Nx);
	TIFFGetField(image, TIFFTAG_IMAGELENGTH, &Ny);
	TIFFGetField(image, TIFFTAG_BITSPERSAMPLE, &BPS);
//	TIFFGetField(image, TIFFTAG_COMPRESSION, &com);
//	TIFFGetField(image, TIFFTAG_PHOTOMETRIC, &pm);
//	TIFFGetField(image, TIFFTAG_SAMPLESPERPIXEL, &spp);
//	TIFFGetField(image, TIFFTAG_ROWSPERSTRIP, &rps);
//	TIFFGetField(image, TIFFTAG_PLANARCONFIG, &pc);
	TIFFGetField(image, TIFFTAG_IMAGEDESCRIPTION, &desc);

	N=Nx;
	
	if(iHead==1){
		TIFFClose(image);
		return 1;
	}
	if ((rline = (unsigned short*)_TIFFmalloc(TIFFScanlineSize(image))) == NULL) {
		printf("cannot allocate memory for line scan\n");
		return 1;
	}

		if (TIFFReadScanline(image, rline, (int)s_layer, 0) < 0) {
			printf("cannot get tif line in %s at line %ld\n", rname, s_layer);
			return 1;
		}
		for (j=0;j<N;j++) data[j] = *(rline+j);

	_TIFFfree(rline);
	TIFFClose(image);

	//printf("%s\t%d\t%d\t%d\t%d\t",rname,Nx,Ny,data[N/3],data[2*N/3]);

	return 0;
}

/*----------------------------------------------------------------------*/
void Store32TiffFile(char *wname, int wX, int wY, int wBPS, float *data32, char *wdesc)
{
	TIFF *image;
	long i;

	image = TIFFOpen(wname, "w");

	TIFFSetField(image, TIFFTAG_IMAGEWIDTH, wX);
	TIFFSetField(image, TIFFTAG_IMAGELENGTH, wY);
	TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, 32);
	TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(image, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP );
	TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, 1);
	TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, wdesc);
	TIFFSetField(image, TIFFTAG_ARTIST, "ct_rec_tif");
//	TIFFSetField(image, TIFFTAG_MINSAMPLEVALUE, mmmin );
//	TIFFSetField(image, TIFFTAG_MAXSAMPLEVALUE, mmmax );

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data32 + i*wX, wX * sizeof(float));
	}

	TIFFClose(image);
}


/*----------------------------------------------------------------------*/

int read_log()
{
	int			i, j;
	int		nnn;
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
	M = 0;
	nnn = 0;
	j   = 0;
	while(j>=0){
		j=j+1;
		lne[0]='\0';
		ss = fgets(lne, 100, f);
		if (i = strlen(lne) > 2){
			sscanf(lne, "%d %f %f %d", &nnn, &shottime[j], &shotangle[j], &flg_I0);
			if (flg_I0 == 0){
				II0[NI0] = nnn;
				NI0=NI0+1;				// number of I_0
			}
			else{
				Ishot[M] = nnn;
				M=M+1;				// number of I
			}
		}
		else{
			break;
		}
	}
	n_total = (short)(nnn + 2);
//	fprintf(stderr, "nshot = %d, NI0 = %d, total = %d \n", M, NI0, n_total);
	fclose(f);
	if(n_total<1000) iFlag = 0;
	if(n_total>=1000) iFlag = 1;
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
	int			i, j, k, jx, nshot, *ilp, iplc;
	double		t1, t2;
	double		a[MAXPIXL], b[MAXPIXL];
	char		ffi01[15], ffi02[15], fii[15];
	double		I01, I02;

	double		*poa;
	double		p_sum,p_ave;
	FILE		*fi;
	char		fname[15];

	TA1 =0.0;
	TA2 =0.0;

	ilp=(int *)malloc(M*sizeof(int));
	iplc=0;

	poa = (double *)malloc(N*sizeof(double));
	for(j=0;j<N;++j) *(poa+j)=0.0;

// counting for number of projections (M = nshot)
	nshot = 0;

// loop between I0_1st and I0_2nd
	for (j=0;j<NI0-1;++j){

//IIO[j] and IIO[j+1] are opened
		if(iFlag==0) sprintf(fname, "%s%03d.tif", flhead, II0[j]);
		if(iFlag==1) sprintf(fname, "%s%04d.tif", flhead, II0[j]);
		if ((i = Read16TiffFileLine(fname, 0)) != 0){
			printf("something wrong -- return value is %d(II01)", i);
			return -1;
		}
		for(i=0;i<N;++i) II01[i]=data[i];
//		printf("\n");
		if(iFlag==0) sprintf(fname, "%s%03d.tif", flhead, II0[j+1]);
		if(iFlag==1) sprintf(fname, "%s%04d.tif", flhead, II0[j+1]);
		if ((i = Read16TiffFileLine(fname, 0)) != 0){
			printf("something wrong -- return value is %d(II02)", i);
			return -1;
		}
		for(i=0;i<N;++i) II02[i]=data[i];
//		printf("\n");
// I0EV
		t1 = shottime[II0[j]];
		t2 = shottime[II0[j+1]];

// 1 layer(ln)
		for (jx=0;jx<N;++jx){
			I01   = (double)II01[jx]-dark[jx];
			I02   = (double)II02[jx]-dark[jx];
			a[jx] = (double)(((double)(I02   - I01))    / (t2 - t1));
			b[jx] = (double)(((double)I01*t2 - (double)I02*t1) / (t2 - t1));
		}
		for (k=II0[j]+1;k<II0[j+1];++k){
			// obtain p(x) from a[jx], b[jx] using shottime[k]
			if(iFlag==0) sprintf(fname, "%s%03d.tif", flhead, k);
			if(iFlag==1) sprintf(fname, "%s%04d.tif", flhead, k);
			if ((i = Read16TiffFileLine(fname, 0)) != 0){
				printf("something wrong -- return value is %d(I)", i);
				return -1;
			}
			for(i=0;i<N;++i) I[i]=data[i];
			*(ilp+nshot)=0;
			for (jx = 0; jx < N; ++jx){
				I0[jx]    = a[jx] * shottime[k] + b[jx];
				if ((I[jx]-dark[jx]) < 1){
					if(*(ilp+nshot)==0){
						printf("Warning \t");
//						printf("  jx = %d, I0 = %f, I = %d, dark = %d, ln =%d \n", jx, I0[jx], I[jx], dark[jx], ln);
//						printf("  II01 = %d, II02 = %d\n", II01[jx], II02[jx]);
//						printf("  t1   = %f, t2   = %f \n", t1, t2);
//						printf("  a = %f,  b = %f\n", a[jx], b[jx]);
						printf("  %s\t black\n", fii);
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
//			printf("%f\t%f\n",*(po+N*nshot+N/3),*(po+N*nshot+2*N/3));

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
				//	*(po+N*nshot+jx) = *(po+N*nshot+jx) - p_ave; /* comment out here */
					TA = TA + *(po+N*nshot+jx);
					*(poa+jx)=*(poa+jx)+*(po+N*nshot+jx);
				}
				TA1=TA1+TA;
				TA2=TA2+TA*TA;
				iplc=iplc+1;
			}
//
			nshot = nshot + 1;
		} // end of k loop
	} // end of j loop


// correcion for black projection
	for(k=0;k<M;++k){
		if(*(ilp+k)==1){
//			printf("%d\t%d\n",k,*(ilp+k) );
			for(j=0;j<N;++j){
				*(po+N*k+j)=*(poa+j)/(double)iplc;
			}
		}
	}

	return (0);
}

/*----------------------------------------------------------------------*/

float calc_center(){
int		i;
long	j, dx, dx0, N1;
double	mm1, mx1, mm2, mx2;
double	*p000, *p180, rsum, rd, rmsd0, rmsd;
float	center;
char		fname[15];

// file name and data for calculation of rotation center (data[y][x])
// odd = I_0, even = transmitted
char	c_file[13];
unsigned short	cent_1[MAXPIXL], cent_2[MAXPIXL];
unsigned short	cent_3[MAXPIXL], cent_4[MAXPIXL];

// read center file
//I_0 file
	if(iFlag==0) sprintf(fname, "%s%03d.tif", flhead, 1);
	if(iFlag==1) sprintf(fname, "%s%04d.tif", flhead, 1);
	if ((i = Read16TiffFileLine(fname, 0)) != 0){
		printf("something wrong -- return value is %d(II01)", i);
		return -1;
	}
	for(i=0;i<N;++i) cent_1[i]=data[i];

// zero degree file(transmitted)
	if(iFlag==0) sprintf(fname, "%s%03d.tif", flhead, 2);
	if(iFlag==1) sprintf(fname, "%s%04d.tif", flhead, 2);
	if ((i = Read16TiffFileLine(fname, 0)) != 0){
		printf("something wrong -- return value is %d(II01)", i);
		return -1;
	}
	for(i=0;i<N;++i) cent_2[i]=data[i];

//I_0 file
	if(iFlag==0) sprintf(fname, "%s%03d.tif", flhead, n_total-1);
	if(iFlag==1) sprintf(fname, "%s%04d.tif", flhead, n_total-1);
	if ((i = Read16TiffFileLine(fname, 0)) != 0){
		printf("something wrong -- return value is %d(II01)", i);
		return -1;
	}
	for(i=0;i<N;++i) cent_3[i]=data[i];

//  180 degree file(transmitted)
	if(iFlag==0) sprintf(fname, "%s%03d.tif", flhead, n_total);
	if(iFlag==1) sprintf(fname, "%s%04d.tif", flhead, n_total);
	if ((i = Read16TiffFileLine(fname, 0)) != 0){
		printf("something wrong -- return value is %d(II01)", i);
		return -1;
	}
	for(i=0;i<N;++i) cent_4[i]=data[i];

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
	double		Clock, t1,t2;					// timer setting
	int			i;
	FILE		*fo;
	Float		**P, **f;			// full size of reconstructed image
	double		delta, theta0;
	double		*rec_temp;
	double		data_max, data_min;
	float		base, f_center, f_size, *frec;
	long		lup;
	double		div;
	double		r0;
	int			vv, hh, jx;			// 
	char	*comm = NULL;
	char	fout[15];
	long	val, j;
	int	ifilter;

// parameter setting
	theta0 = 0.0;
	f_size = 1.0;
	if (argc==5){
		theta0 = M_PI/180.0*atof(argv[4]);
		f_center = atof(argv[2]);
		f_size = atof(argv[3]);
	}
	else if (argc==4){
		f_center = atof(argv[2]);
		f_size = atof(argv[3]);
	}
	else if (argc==3){
		f_center = atof(argv[2]);
	}
	else if (argc==2){
	}
	else{
		fprintf(stderr, "parameter was wrong!!!\n");
		fprintf(stderr, "usage : %s layer (center) (pixel size) (offsetangle)\n", argv[0]);
		fprintf(stderr, "default pixel size 1.0um\n");
		return(1);
	}
	s_layer = atoi(argv[1]);

// read shot log file
	if ((i = read_log()) != 0){
		fprintf(stderr, "canot read 'output.log' ! (%d)\n", i);
		return(1);
	}
// read dark image
	if ((i = Read16TiffFileLine(darkfile, 0)) != 0){
		fprintf(stderr, "something wrong in dark file (%d)\n", i);
		return(1);
	}
	for(i=0;i<N;++i) dark[i]=data[i];
	printf("\n");

// p initialization
	po = (double *)malloc(N*M*sizeof(double));

	cent_flag = 0;
	Clock=CLOCK();
	if((i=StoreProjection(s_layer)) !=0){
		printf("something wrong in StoreProjection (%d)\n",i);
		return(1);
	}
	TAM=TA1/M;
	TASD=sqrt(TA2/M-TAM*TAM);
//	fprintf(stderr,"TA mean %lf\tTASD %lf\tTASD/TA %lf\n", TAM, TASD, TASD/TAM);
	t1=CLOCK()-Clock;

	if(argc<3) {
		f_center= calc_center();
	}

// CBP
//	printf("init cbp ");
	if ((P=InitCBP(N,M))==NULL) Error("memory allocation error.");
//	printf("done \n");

	for (j=0; j<M; j++){
		for (i=0; i<N; i++){
			if(abs(*(po+N*j+i))<100){
				P[j][i]=*(po+N*j+i);
			}
		}
	}
/* ----------------  ring removal start ---------------- */
/*                                                       */
    int kernel_size = 5; // Default kernel size
    int num_threads = 40; // Default number of threads
    float		*image_data = NULL, *result_data = NULL;
    // Get kernel size from environment variable
    kernel_size = get_kernel_size_from_env();
	// Get number of threads from environment variable
    num_threads = get_num_threads_from_env();
    // Allocate memory
	image_data = (float *)malloc(N * M * sizeof(float));
	result_data = (float *)malloc(N * M * sizeof(float));

	for (j=0; j<M; j++){
		for (i=0; i<N; i++){
			*(image_data+N*j+i)=P[j][i];
		}
	}
	// Execute OpenMP image processing
	if (sort_filter_restore_omp(image_data, result_data, N, M, kernel_size, num_threads) != 0) {
		fprintf(stderr, "OpenMP image processing failed\n");
		return 5;
	}
	for (j=0; j<M; j++){
		for (i=0; i<N; i++){
				P[j][i]=*(result_data+N*j+i);
		}
	}
    if (image_data) free(image_data);
    if (result_data) free(result_data);
/* ----------------  ring removal finish --------------- */
/*                                                       */


//
	r0=-1.0*f_center;
	delta = 1.0; //normilized data

//	printf("cbp ");
// measure Convolution time
	Clock=CLOCK();
	f=CBP(delta,r0,theta0);
//	printf("done \n");

	rec_temp= (double *)malloc(N*N*sizeof(double));
	frec= (float *)malloc(N*N*sizeof(float));
	data_max =-32000.;
	data_min = 32000.;

	for (vv=0; vv<N; vv++){
		for (hh=0; hh<N; hh++){
			*(rec_temp+N*vv+hh) = f[vv][hh]*10000./f_size;	/* unit change  um -> cm */;
			if(data_max<*(rec_temp+N*vv+hh)) data_max=*(rec_temp+N*vv+hh);
			if(data_min>*(rec_temp+N*vv+hh)) data_min=*(rec_temp+N*vv+hh);
		}
	}
	t2=CLOCK()-Clock;

//	printf("MAX and min %f and %f\n", (float)data_max, (float)data_min);
//	printf("Image size:\t%d x %d\n", N, N);
//	printf("pixel size:\t%lf um\n",f_size);

// initialize fom
	jx=0;
	for(vv=0;vv<N;vv++){
		for(hh=0;hh<N;hh++){
			*(frec+jx)=(float)*(rec_temp+jx);
			jx=jx+1;
		}
	}

	if ((comm=(char *)malloc(150))==NULL)
		Error("comment memory allocation error.");
	sprintf(comm,"%f\t%f\t%d\t%f\t%f\t%f",f_size, f_center, M, theta0, (float)data_min, (float)data_max);
//	printf("%f\t%f\t%d\t%f\t%f\t%f",f_size, f_center, M, theta0, (float)data_max, (float)data_min);
	//sprintf(desc, "%lf\t%lf\n%d\t%lf\n%lf\t%lf\n%lf\t%lf", delta, r0, M, theta0, f1, f2, f0, df);

	(void)sprintf(fout, "rec%05d.tif", (int)s_layer);
	printf("%s\r",fout);
	Store32TiffFile(fout, N, N, 32, frec, comm);
	printf("%s\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n",fout,f_size, f_center, M, theta0, (float)data_min, (float)data_max, (float)t1, (float)t2);
//	printf("Stored to:\t%s\n",fout);
//	fprintf(stderr,"Store Sinogram:\t%lf / sec\n",t1);
//	fprintf(stderr,"reconstruction:\t%lf / sec\n",t2);

	free(po); free(rec_temp); free(frec);

	// append to log file
	FILE		*ff;
	if((ff = fopen("../cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(ff,"%s ",argv[i]);
	fprintf(ff,"\n");
	fclose(ff);

	return 0;
}
