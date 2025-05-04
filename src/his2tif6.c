/* split HIS format file to img files  */
/* usage is 'his2tif inputfile (head)' */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "tiffio.h"
//#include "cell.h"
//#include "sif.h"

#define INTEL
#define HIS_Header_Size 64

#define MA(cnt,ptr)	malloc((cnt)*sizeof(*(ptr)))

struct HIS_Header
{
	char			head[2];			/*0-1*/
	short			comment_length;		/*2-3*/
	short			width;				/*4-5*/
	short			height;				/*6-7*/
	short			x_offset;			/*8-9*/
	short			y_offset;			/*10-11*/
	short			type;				/*12-13*/
	unsigned short			n_image1;			/*14-17*/
	unsigned short			n_image2;			/*14-17*/
	short			reserve1;			/*18-19*/
	short			reserve2;			/*20-21*/
	double			time_stamp;			/*22-29*/
	long			maker;				/*30-33*/
	char			reserved[30];		/*34-63*/
	char			*comment;
};
typedef struct HIS_Header HISHeader;


unsigned short	*data;
unsigned char	*cdata;
long			NP;
char	flhead[]   = "a";

static FILE *fi1, *fo;
FILE	*fi;

static void	Error(msg)
char		*msg;
{
	fputs(msg,stderr); fputc('\n',stderr); exit(1);
}

/*----------------------------------------------------------------------*/
int read_hishead(fname,his)
HISHeader	*his;
char *fname;
{
	int		j;
	FILE	*fi;
	char	buffer[HIS_Header_Size];

//open input files
	if((fi = fopen(fname,"rb")) == NULL){
		printf("can not open %s for input\n", fname);
		free(his->comment); 
		return(-10);
	}
// reading header without comment and setup 
	if (fread(buffer, sizeof(char), HIS_Header_Size, fi) != HIS_Header_Size){
		fprintf(stderr, "Error [his header 64bytes]t\n"); 
		return(-21);
	}
	memcpy(&his->head,           buffer,       2);	// char x 2
	memcpy(&his->comment_length, buffer +  2,  2);	// short
	memcpy(&his->width,          buffer +  4,  2);	// short
	memcpy(&his->height,         buffer +  6,  2);	// short
	memcpy(&his->x_offset,       buffer +  8,  2);	// short
	memcpy(&his->y_offset,       buffer + 10,  2);	// short
	memcpy(&his->type,           buffer + 12,  2);	// short
	memcpy(&his->n_image1,       buffer + 14,  2);	// short
	memcpy(&his->n_image2,       buffer + 16,  2);	// short
	memcpy(&his->reserve1,       buffer + 18,  2);	// short
	memcpy(&his->reserve2,       buffer + 20,  2);	// short
	memcpy(&his->time_stamp,     buffer + 22,  8);	// double
	memcpy(&his->maker,          buffer + 30,  4);	// long
	memcpy(&his->reserved,       buffer + 34, 30);	// char x 30


	his->comment=(char *)malloc(his->comment_length+1);
	if (his->comment == NULL){
		fclose(fi); free(his->comment); 
		fprintf(stderr, "Connot allocate memory (HIS_Header.comment)\n"); 
		return(-12);
	}
	if (fread(his->comment, sizeof(char), his->comment_length, fi) != his->comment_length){ 
		fclose(fi); free(his->comment); 
		fprintf(stderr, "EOF found (comment)\n"); 
		return(-13);
	}
	his->comment[his->comment_length] = '\0';

// check for correct reading of HIS file
// printf("-------------------------------\n");
// printf("%s\n", his->head);
// printf("%d\t", his->comment_length);
// printf("%d\t", his->width);
// printf("%d\t", his->height);
// printf("%d\t", his->x_offset);
// printf("%d\t", his->y_offset);
// printf("%d\t", his->type);
// printf("%d\t", his->n_image1);			/*14-17*/
// printf("%ld\t", his->n_image2*65536);			/*14-17*/
// printf("%d\t", his->reserve1);			/*18-19*/
// printf("%d\t", his->reserve2);			/*20-21*/
//printf("%lf\t", his->time_stamp1);			/*22-29*/
//printf("%lf\t", his->time_stamp2);			/*22-29*/
// printf("%ld\t", his->maker);				/*30-33*/
// printf("%s\n", his->reserved);	/*34-64*/
// printf("%s\n", his->comment);
// printf("-------------------------------\n");
	
	
	fclose(fi);
	return(0);
}

/*----------------------------------------------------------------------*/
int read_his(his)
HISHeader	*his;
{
	int		j;

// read comment and image data from file1
	if (fread(his, sizeof(char), HIS_Header_Size, fi) != HIS_Header_Size){
		fprintf(stderr, "EOF in input file (header)\n"); free(his->comment); return(-1);
	}
	if (strncmp(his->head, "IM", 2)){
		fclose(fi); free(his->comment);
		fprintf(stderr, "This is not HIS image file.\n"); 
		return(-11);
	}

	his->comment=(char *)malloc(his->comment_length+1);
	if (his->comment == NULL){
		fclose(fi); free(his->comment); 
		fprintf(stderr, "Connot allocate memory (HIS_Header.comment)\n"); 
		return(-12);
	}
	if (fread(his->comment, sizeof(char), his->comment_length, fi) != his->comment_length){ 
		fclose(fi); free(his->comment); 
		fprintf(stderr, "EOF found (comment)\n"); 
		return(-13);
	}
	his->comment[his->comment_length] = '\0';

	if(his->type==2){
		if ((j = fread(data, sizeof(unsigned short), NP, fi)) != NP){
			fclose(fi); free(his->comment); 
			fprintf(stderr, " Error reading at %d\n",  j);
			return(-14);
		}
	}
	if(his->type==6){
		if ((j = fread(cdata, sizeof(unsigned char), NP, fi)) != NP){
			fclose(fi); free(his->comment); 
			fprintf(stderr, " Error reading at %d\n",  j);
			return(-14);
		}
	}

	return(0);
	}

/*----------------------------------------------------------------------*/

void Store16TiffFile(char *wname, int wX, int wY, int wBPS, unsigned short *data16, char *wdesc)
{
	TIFF *image;
	long i;

	image = TIFFOpen(wname, "w");

	TIFFSetField(image, TIFFTAG_IMAGEWIDTH, wX);
	TIFFSetField(image, TIFFTAG_IMAGELENGTH, wY);
	TIFFSetField(image, TIFFTAG_BITSPERSAMPLE, wBPS);
	TIFFSetField(image, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
	TIFFSetField(image, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
	TIFFSetField(image, TIFFTAG_SAMPLESPERPIXEL, 1);
	TIFFSetField(image, TIFFTAG_ROWSPERSTRIP, 1);
	TIFFSetField(image, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(image, TIFFTAG_IMAGEDESCRIPTION, wdesc);
	TIFFSetField(image, TIFFTAG_ARTIST, "his2tif6");

	for (i = 0; i<wY; i++) {
		TIFFWriteRawStrip(image, i, data16 + i*wX, wX * sizeof(unsigned short));
	}

	TIFFClose(image);
}


int	main(argc,argv)
int		argc;
char	*argv[];
{
	unsigned short	*image1, *outimg;
	long	i, j, jj;
	char	readfile[20], outfile[20];
	char	*cmt;
	short	clg;
	unsigned long	i_end, k;
	double ttt, fff;
	unsigned	bits;

//	Cell		**cell;
	short	Nx, Ny;

	int		y,x;
	int		dBPS;

	char	*comm = NULL;
	HISHeader		his;

	if(argc < 2){   // when bad number of parameters
		printf("Usage is 'his2tif hisfile (head)'\n");
	}
	if (argc==3){
		sprintf(flhead, "%s", argv[2]);
	}

	sprintf(readfile, ("%s"), argv[1]);

	if ((j = read_hishead(readfile, &his)) != 0){
		printf("something wrong -- return value is %d", j);
		free(his.comment); 
		return(-1);
	}
	i_end=his.n_image1+65536*his.n_image2;
	if(his.type==6) printf("12 bits %ld images\n", i_end);
	if(his.type==2) printf("16 bits %ld images\n", i_end);

	Nx=his.width;
	Ny=his.height;
	
	printf("reading %s.\n",readfile);

// initialize
	
	if(his.type==2){
		NP=Nx*Ny;
		data=(unsigned short *)malloc(NP*sizeof(unsigned short));
	}
	if(his.type==6){
		NP=(Nx*Ny)*3/2;
		cdata=(unsigned char *)malloc(NP*sizeof(unsigned char));
	}
//		if ((cell=(Cell **)MA(Ny,cell))==NULL)
//		    Error("cell memory allocation error (1).");
//		for (y=0; y<Ny; y++) {
//		    if ((cell[y]=(Cell *)MA(Nx,*cell))==NULL)
//			Error("cell memory allocation error (2).");
//		}

	bits=16;
	dBPS=16;
	

//open input files
	if((fi = fopen(argv[1],"rb")) == NULL){
		printf("can not open %s for input\n", argv[1]);
		return(-10);
	}
//	for(k=0;k<i_end;++k){
	for(k=0; !feof(fi) ;++k){
//		printf("%d\n",k);

	// read input images
		if ((j = read_his(&his)) != 0){
			printf("something wrong -- return value is %d", j);
			free(his.comment); 
			return(-1);
		}
		
		outimg = (unsigned short *) malloc(Nx*Ny*sizeof(unsigned short));
	
		if(his.type==2){
			for(jj=0;jj<NP;++jj){
				*(outimg+jj)=*(data+jj);
			}
		}
		if(his.type==6){
			j=0;
			for(jj=0;jj<NP;++jj){
				if ((jj%3)==0){
					*(outimg+j)=16**(cdata+jj)+*(cdata+jj+1)/16;
//					printf("%d\t%d\t%d\n",j,jj,*(outimg+j));
					j=j+1;
				}
				if ((jj%3)==2){
					*(outimg+j)=((*(cdata+jj-1)%16)*256)+*(cdata+jj);
//					printf("%d\t%d\t%d\n",j,jj,*(outimg+j));
					j=j+1;
				}
			}
		}

		if ((comm=(char *)malloc(150))==NULL)
			Error("comment memory allocation error.");
		sprintf(comm,"his2tif\t%ld\t/\t%ld",k,i_end);

		sprintf(outfile, ("a%06d.tif"), k+1);
		
		
//		j=0;
//		for (y=0; y<Ny; y++){
//			for (x=0; x<Nx; x++){
//				cell[y][x]=*(outimg+j);
//				j=j+1;
//			}
//		}
	

//		StoreTiffFile(outfile,his.width,his.height,16,outimg,comm);
//		StoreImageFile(outfile,(int)his.width,(int)his.height,(int)bits,cell,comm);
		Store16TiffFile(outfile, Nx, Ny, dBPS, outimg, comm);

//printf("%g\t", his.time_stamp);			/*22-29*/
//printf("%d\t", his.time_stamp2);			/*22-29*/
//printf("%d\t", his.time_stamp3);			/*22-29*/
//printf("%d\t", his.time_stamp4);			/*22-29*/
//printf("%d\t", his.time_stamp5);			/*22-29*/
//printf("%d\t", his.time_stamp6);			/*22-29*/
//printf("%d\t", his.time_stamp7);			/*22-29*/
//printf("%d\t", his.time_stamp8);			/*22-29*/

		free(outimg);
		free(comm);

		printf("a%06d.tif \r",k+1);
	}
	printf("\n");
	
//finish!
	free(data);
	free(his.comment);
//	for(y=0;y<Ny;y++) free(cell[y]);
//	free(cell);
	fclose(fi);
	return(0);
}
