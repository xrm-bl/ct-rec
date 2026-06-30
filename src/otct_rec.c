// program ofct_sino
// 
// Required files are q???.tif, dark.tif.
// output file "s????.sin"

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "cbp.h"
#include <math.h>
#include "tiffio.h"
#ifdef USE_GPU
  #include "sort_filter_g.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_gpu
#else
  #include "sort_filter_omp.h"
  #define SORT_FILTER_RESTORE sort_filter_restore_omp
#endif

/*----------------------------------------------------------------------*/
#ifndef M_PI
#define M_PI				3.1415926535897932385
#endif
#define INTEL
#define MAX_SHOT			20000
#define MAXPIXL				14001

/*----------------------------------------------------------------------*/

// main data for read transmitted images (data[y][x])
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

// file head character  q???.tif or q????.tif
char	flhead[]   = "q";
char	darkfile[16] = "dark.tif";

int useTiff = 0;   // 0 = HiPic .img, 1 = 16-bit TIFF .tif (auto-detected)

// output file name
char	fn[256];

// flag for q001.tif or q0001.tif
int		iFlag, nfq, BPS;

// projection
unsigned short	N, n_total, NN, NNST;

// calculation layer
long	s_layer;

int		Nx, Ny;

double	*pp;

// total absorption
double	TA, TA1, TA2, TAM, TASD;

// offset CT axis
float	dcnt;
int		idcnt;
long	ldcnt;

/*----------------------------------------------------------------------*/

static void Error(msg)
char        *msg;
{
	fputs(msg,stderr);
	fputc('\n',stderr);
	exit(1);
}


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
/*----------------------------------------------------------------------*/
/* HiPic (.img) header + single-line reader (added for img/tif switch).  */
#define HiPic_Header_Size  64
struct HiPic_Header{
    char            head[2];
    unsigned short  comment_length;
    short           width;
    short           height;
    short           x_offset;
    short           y_offset;
    short           type;
    char            reserved[63 - 13];
    char            *comment;
};
typedef struct HiPic_Header Header;

/* Read one layer line from a HiPic .img into the global data[]; mirrors
   Read16TiffFileLine (fills data[], sets N, iHead==1 -> size only).      */
int Read16ImgFileLine(char *rname, int iHead)
{
    FILE   *fi;
    Header h;
    int    j, i_res;

    if((fi = fopen(rname, "rb")) == NULL){
        printf("can not open %s for input\n", rname);
        return(-10);
    }
    if(fread(&h, sizeof(char), HiPic_Header_Size, fi) != HiPic_Header_Size){
        fclose(fi); fprintf(stderr, "EOF in %s (header)\n", rname); return(-1);
    }
    if(strncmp(h.head, "IM", 2)){
        fclose(fi); fprintf(stderr, "not a HiPic image: %s\n", rname); return(-11);
    }
    h.comment = (char *)malloc(h.comment_length + 1);
    if(h.comment == NULL){ fclose(fi); return(-12); }
    if(fread(h.comment, sizeof(char), h.comment_length, fi) != (unsigned short)h.comment_length){
        free(h.comment); fclose(fi); fprintf(stderr, "EOF (comment) in %s\n", rname); return(-13);
    }
    h.comment[h.comment_length] = '\0';
    N = h.width;
    if(iHead == 1){ free(h.comment); fclose(fi); return 1; }

    i_res = fseek(fi, (long)h.width * s_layer * sizeof(short), SEEK_CUR);
    if(i_res != 0){
        free(h.comment); fclose(fi);
        fprintf(stderr, " seek error in %s\n", rname); return((int)s_layer + 1);
    }
    if((j = fread(data, sizeof(short), h.width, fi)) != h.width){
        free(h.comment); fclose(fi);
        fprintf(stderr, " read error in %s (%d)\n", rname, j); return((int)s_layer + 1);
    }
    free(h.comment);
    fclose(fi);
    return 0;
}

/*----------------------------------------------------------------------*/
/* Dispatch to .tif or .img reader based on the detected format.
   In img mode a ".tif" request name is mapped to ".img".                 */
int ReadLine(char *rname, int iHead)
{
    char img[1024];
    int  L;

    if(useTiff) return Read16TiffFileLine(rname, iHead);

    snprintf(img, sizeof(img), "%s", rname);
    L = (int)strlen(img);
    if(L >= 4 && (strcmp(img + L - 4, ".tif") == 0 || strcmp(img + L - 4, ".TIF") == 0))
        strcpy(img + L - 4, ".img");
    return Read16ImgFileLine(img, iHead);
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
	short		nnn;
	FILE		*f;
	char		line[100], *ss;
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
		line[0]='\0';
		ss = fgets(line, 100, f);
		if (i = strlen(line) > 2){
			sscanf(line, "%hd %f %f %d", &nnn, &shottime[j], &shotangle[j], &flg_I0);
			if (flg_I0 == 0){
				II0[NI0] = nnn;
				NI0 = NI0 + 1;				// number of I_0 (direct beam)
			}
			else{
				Ishot[NST] = nnn;
				NST = NST + 1;				// number of I (sample image)
			}
		}
		else{
			break;
		}
	}
	n_total = (short)(nnn + 2);
	fclose(f);
	if(n_total<1000) iFlag = 0;
	if(n_total>=1000) iFlag = 1;
	if(n_total>MAX_SHOT){
		fprintf(stderr, "Too many projections!");
		return(1);
	}
	if(shotangle[n_total-2]<350.){
		fprintf(stderr, "It seems it's not a offset CT data.\n");
		fprintf(stderr, "%f\n", shotangle[n_total]);
		return 1;
	}
	NNST=NST/2;
	return(0);
}

/*----------------------------------------------------------------------*/

static int StoreProjection(ln)
long		ln;
{
	int			i, j, k, jx, nshot, ilp, idcN, l;
	double		t1, t2;
	double		a[MAXPIXL], b[MAXPIXL];
	char		ffi01[1024], ffi02[25], fii[25];
	double		I01, I02;

	double		*po;
	double		p_sum, p_ave, ddv, dva, dvb;
	double		black_thresh;
	double		id_sum;
	char		*env_bt;
	FILE		*fi;

	env_bt = getenv("CT_REC_BLACK_THRESH");
	if (env_bt != NULL){
		black_thresh = atof(env_bt);
	} else {
		black_thresh = 1.0;
	}

// initialize po and pp 
	po = (double *)malloc(N*NST*sizeof(double));
	pp = (double *)malloc(NN*NNST*sizeof(double));

	TA1 =0.0;
	TA2 =0.0;

//	idcN=N-idcnt;
	idcN=10;

// counting for number of projections (nshot = NST)
	nshot = 0;

// loop between I_0 and I_0 for reading sample images
	for (j=0;j<NI0-1;++j){

// open IIO[j] and IIO[j+1] 
		if(iFlag==0) snprintf(ffi01, sizeof(ffi01), "%s%03d.tif", flhead, II0[j]);
		if(iFlag==1) snprintf(ffi01, sizeof(ffi01), "%s%04d.tif", flhead, II0[j]);
		if ((i = ReadLine(ffi01, 0)) != 0){
			printf("something wrong -- return value is %d(II01)", i);
			return -1;
		}
		for(i=0;i<N;++i) II01[i]=data[i];
		if(iFlag==0) sprintf(ffi02, "%s%03d.tif", flhead, II0[j+1]);
		if(iFlag==1) sprintf(ffi02, "%s%04d.tif", flhead, II0[j+1]);
		if ((i = ReadLine(ffi02, 0)) != 0){
			printf("something wrong -- return value is %d(II02)", i);
			return -1;
		}
		for(i=0;i<N;++i) II02[i]=data[i];
// I0EV の係数を求めるための準備
		t1 = shottime[II0[j]];
		t2 = shottime[II0[j+1]];

// 1 layer(ln) だけI_0の値を内挿する
		for (jx = 0; jx < N; ++jx){
			I01   = (double)II01[jx]-dark[jx]; // 前後のI_0で一次補間
			I02   = (double)II02[jx]-dark[jx];
			a[jx] = (double)(((double)(I02   - I01))    / (t2 - t1));      //補間係数を求める
			b[jx] = (double)(((double)I01*t2 - (double)I02*t1) / (t2 - t1));
		}
		for ( k = II0[j] + 1; k < II0[j+1]; ++k){
// a[jx], b[jx] から shottime[k] を用いて p(x) を求める。
			if(iFlag==0) sprintf(fii, "%s%03d.tif", flhead, k);
			if(iFlag==1) sprintf(fii, "%s%04d.tif", flhead, k);
//			if ((i = read_tiff_line(fii, &I[0], ln)) != 0){
//				printf("something wrong -- return value is %d(II01)", i);
//				exit(-1);
//			}
			if ((i = ReadLine(fii, 0)) != 0){
				printf("something wrong -- return value is %d(I)", i);
				return -1;
			}
			for(i=0;i<N;++i) I[i]=data[i];
			ilp=0;
			for (jx = 0; jx < N; ++jx){
				I0[jx]    = a[jx] * shottime[k] + b[jx];
			}
			/* black check: average of (I[jx]-dark[jx]) over all pixels */
			id_sum = 0.0;
			for (jx = 0; jx < N; ++jx){
				id_sum += (double)(I[jx] - dark[jx]);
			}
			if (id_sum / (double)N < black_thresh){
				printf("Warning \t");
				printf("  %d\t black \t avg=%.2f (thresh=%.2f)\n",
					k, id_sum/(double)N, black_thresh);
				ilp=1;
			}
			for (jx = 0; jx < N; ++jx){
				*(po+nshot*N+jx)= log((double)(I0[jx])/(double)(I[jx]-dark[jx]));
			}
			if(ilp==1){
				for(jx=0;jx<N;++jx){
					*(po+nshot*N+jx)=0.0;
				}
			}

			p_sum=0.0;
			p_ave=0.0;
			for(jx=0;jx<10;++jx){
				p_sum=p_sum+*(po+nshot*N+jx);
			}
			p_ave=p_sum/10.;
			TA=0.0;
			for(jx=0;jx<N;++jx){
//				*(po+nshot*N+jx)=*(po+nshot*N+jx)-p_ave; /* comment out here */
				TA      =TA+*(po+nshot*N+jx);
			}
			TA1=TA1+TA;
			TA2=TA2+TA*TA;
//
			nshot = nshot + 1;
		} // end of k loop
	} // end of j loop
//	printf("1\n");

//
	ddv=(double)((2*idcN)-1);
	for(j=0;j<NNST;++j){
		k=0;
//		printf("%d\n",j);
		for(jx=0;jx<idcnt-idcN;++jx){
			*(pp+j*NN+k)=*(po+j*N+jx);
			k=k+1;
		}
// simple avarage
//		for(jx=-idcN;jx<idcN;++jx){
//			*(pp+k)=(*(po+j*N+jx+idcnt)+*(po+(j+NNST)*N+idcnt-jx-1))/2.0;
//			k=k+1;
//		}
		l=0;
// linear interpolation
		for(jx=-idcN;jx<idcN;++jx){
			dva=ddv-(double)l;
			dvb=(double)l;
			*(pp+j*NN+k)=(*(po+j*N+jx+idcnt)*dva+*(po+(j+NNST)*N+idcnt-jx-1)*dvb)/ddv;
			k=k+1;l=l+1;
		}
		for(jx=0;jx<idcnt-idcN;++jx){
			*(pp+j*NN+k)=*(po+(j+NNST)*N+idcnt-idcN-jx-1);
			k=k+1;
		}

//		for(jx=0;jx<idcnt+1;++jx){
//			*(pp+k)=*(po+j*N+jx);
//			k=k+1;
//		}
//		for(jx=1;jx<idcnt+1;++jx){
//			*(pp+k)=*(po+(j+NNST)*N+idcnt-jx);
//			k=k+1;
//		}
// roupe for correction
//		p_sum = 0.0;
//		p_ave = 0.0;
//		for (jx=0; jx<10; ++jx){
//			p_sum = p_sum + *(pp+jx) + *(po+NN-1-jx);
//		}
//		p_ave = p_sum / (10.-0.) / 2.;
//		for (jx=0; jx<NN; ++jx){
//			*(pp+jx) = *(pp+jx) - p_ave; /* comment out here */
//		}
		
//		if (i = fwrite(pp , sizeof(double), NN, fi) > NN*sizeof(double)){
//			fprintf(stderr, " Error writing data to %s at %d\n", fn, i);
//			fclose(fi);
//			return(i+1);
//		}
//		if(j==0){
//			for(jx=0;jx<idcnt+1;++jx){
//				printf("%lf\t%lf\n",*(pp+jx),*(po+j*N+jx));
//			}
//			for(jx=1;jx<idcnt+1;++jx){
//				printf("%lf\t%lf\n",*(pp+idcnt+jx),*(po+(j+NNST)*N+idcnt-jx));
//			}
//		}
	}



	
	free(po);
//	free(pp);
//	fclose(fi);
	return (0);
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
	int			i, j;
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
	char	fout[1024];
	long	val;
	int	ifilter;

// parameter setting
	theta0 = 0.0;
	f_size = 1.0;
	if (argc==5){
		theta0 = M_PI/180.0*atof(argv[4]);
		dcnt = atof(argv[2]);
		f_size = atof(argv[3]);
	}
	else if (argc==4){
		dcnt = atof(argv[2]);
		f_size = atof(argv[3]);
	}
	else if (argc==3){
		dcnt = atof(argv[2]);
	}
	else{
		fprintf(stderr, "parameter was wrong!!!\n");
		fprintf(stderr, "usage : %s layer offset  (pixel size) (offsetangle)\n", argv[0]);
		return(1);
	}
	s_layer = atoi(argv[1]);

    /* detect input format: dark.img -> HiPic .img, else dark.tif -> 16-bit TIFF */
    {
        FILE *ftest;
        if((ftest = fopen("dark.img", "rb")) != NULL){
            fclose(ftest); useTiff = 0; snprintf(darkfile, sizeof(darkfile), "dark.img");
        }
        else if((ftest = fopen("dark.tif", "rb")) != NULL){
            fclose(ftest); useTiff = 1; snprintf(darkfile, sizeof(darkfile), "dark.tif");
        }
        else{
            fprintf(stderr, "dark.img / dark.tif not found.\n");
            return(1);
        }
    }
	f_center=dcnt+0.5;
		
	NN=2*dcnt;
	idcnt=(short)dcnt;
	ldcnt=(long)(dcnt*10);

// read shot log file
	if ((i = read_log()) != 0){
		fprintf(stderr, "canot read 'output.log' ! (%d)\n", i);
		return(1);
	}
//	fprintf(stderr, " nshot = %d, NI0 = %d, total = %d \n", NST/2, NI0, n_total);

// read dark image
	if ((i = ReadLine(darkfile, 0)) != 0){
		fprintf(stderr, "something wrong in dark file (%d)\n", i);
		return(1);
	}
	for(i=0;i<N;++i) dark[i]=data[i];
	printf("\n");
//	fprintf(stderr, "width of projection data = %d \n",  h.width);


	cent_flag = 0;
	Clock=CLOCK();
	if((i=StoreProjection(s_layer)) !=0){
		printf("something wrong in StoreProjection (%d)\n",i);
		return(1);
	}
	TAM=TA1/NST;
	TASD=sqrt(TA2/NST-TAM*TAM);
//	fprintf(stderr,"TA mean %lf\tTASD %lf\tTASD/TA %lf\n", TAM, TASD, TASD/TAM);

	t1=CLOCK()-Clock;

// CBP
//	printf("init cbp ");
	if ((P=InitCBP(NN,NNST))==NULL) Error("memory allocation error.");
//	printf("done \n");

	for (j=0; j<NNST; j++){
		for (i=0; i<NN; i++){
			if(abs(*(pp+NN*j+i))<100){
				P[j][i]=*(pp+NN*j+i);
//				printf("%d\t%d\t%lf\t%lf\r",i,j,P[j][i],*(pp+NN*j+i));
			}
		}
	}
	if(P[0][100]==0.0){
		for(j=0;j<N;++j){
			P[0][i]=P[1][i];
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
	image_data = (float *)malloc(NN * NNST * sizeof(float));
	result_data = (float *)malloc(NN * NNST * sizeof(float));

	for (j=0; j<NNST; j++){
		for (i=0; i<NN; i++){
			*(image_data+NN*j+i)=P[j][i];
		}
	}
	// Execute OpenMP image processing
	if (SORT_FILTER_RESTORE(image_data, result_data, NN, NNST, kernel_size, num_threads) != 0) {
		fprintf(stderr, "OpenMP image processing failed\n");
		return 5;
	}
	for (j=0; j<NNST; j++){
		for (i=0; i<NN; i++){
				P[j][i]=*(result_data+NN*j+i);
		}
	}
    if (image_data) free(image_data);
    if (result_data) free(result_data);
/* ----------------  ring removal finish --------------- */
/*                                                       */

	r0=-1.0*f_center;
	delta = 1.0; //normilized data

	printf("\n");
//	printf("cbp ");
// measure Convolution time
	Clock=CLOCK();
	f=CBP(delta,r0,theta0);
//	printf("done \n");

	rec_temp= (double *)malloc(NN*NN*sizeof(double));
	frec= (float *)malloc(NN*NN*sizeof(float));
	data_max =-32000.;
	data_min = 32000.;

	for (vv=0; vv<NN; vv++){
		for (hh=0; hh<NN; hh++){
			*(rec_temp+NN*vv+hh) = f[vv][hh]*10000./f_size;	/* unit change  um -> cm */;
			if(data_max<*(rec_temp+NN*vv+hh)) data_max=*(rec_temp+NN*vv+hh);
			if(data_min>*(rec_temp+NN*vv+hh)) data_min=*(rec_temp+NN*vv+hh);
		}
	}
	t2=CLOCK()-Clock;
	
// initialize fom
	jx=0;
	for(vv=0;vv<NN;vv++){
		for(hh=0;hh<NN;hh++){
			*(frec+jx)=(float)*(rec_temp+jx);
			jx=jx+1;
		}
	}

	if ((comm=(char *)malloc(150))==NULL)
		Error("comment memory allocation error.");
	sprintf(comm,"%f\t%f\t%d\t%f\t%f\t%f",f_size, f_center, NNST, theta0, (float)data_min, (float)data_max);

	(void)snprintf(fout, sizeof(fout), "rec%05d.tif", (int)s_layer);
	printf("%s\r",fout);
	Store32TiffFile(fout, NN, NN, 32, frec, comm);
	printf("%s\t%f\t%f\t%d\t%f\t%f\t%f\t%f\t%f\n",fout, f_size, f_center, NNST, theta0, (float)data_min, (float)data_max, (float)t1, (float)t2);
	
	free(pp); free(rec_temp); free(frec);
	
	// append to log file
	FILE		*ff;
	if((ff = fopen("../cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(ff,"%s ",argv[i]);
	fprintf(ff, "   %% kernel_size %d", kernel_size);
	fprintf(ff,"\n");
	fclose(ff);

	return 0;
}
