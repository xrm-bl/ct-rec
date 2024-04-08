/* split HIS format file to img files  */
/* usage is 'his_spl_E inputfile partnum loopnum' */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define INTEL
#define HIS_Header_Size 64
#define IMG_Header_Size	64

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

struct IMG_Header{
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
typedef struct IMG_Header IMGHeader;


unsigned short	*data;
unsigned char	*cdata;
unsigned char	*odata;
long			NP;
char	flhead[]   = "a";

FILE	*fo;
FILE	*fi;

/*----------------------------------------------------------------------*/
int read_hishead(fname,his)
char *fname;
HISHeader	*his;
{
	int		j;
	FILE	*f;
	char	buffer[HIS_Header_Size];

//open input files
	if((f = fopen(fname,"rb")) == NULL){
		printf("can not open %s for input\n", fname);
		free(his->comment); 
		return(-20);
	}
// reading header without comment and setup 
	if (fread(buffer, sizeof(char), HIS_Header_Size, f) != HIS_Header_Size){
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
		fclose(f); free(his->comment); 
		fprintf(stderr, "Connot allocate memory (HIS_Header.comment)\n"); 
		return(-22);
	}
	if (fread(his->comment, sizeof(char), his->comment_length, f) != his->comment_length){
		fclose(f); free(his->comment); 
		fprintf(stderr, "EOF found (comment)\n"); 
		return(-23);
	}
	his->comment[his->comment_length] = '\0';

// check for correct reading of HIS file
//printf("-------------------------------\n");
//printf("%s\n", his->head);				/*0-1*/
//printf("%d\t", his->comment_length);		/*2-3*/
//printf("%d\t", his->width);				/*4-5*/
//printf("%d\t", his->height);				/*6-7*/
//printf("%d\t", his->x_offset);			/*8-9*/
//printf("%d\t", his->y_offset);			/*10-11*/
//printf("%d\t", his->type);				/*12-13*/
//printf("%d\t", his->n_image1);			/*14-15*/
//printf("%ld\t", his->n_image2*65536);		/*16-17*/
//printf("%d\t", his->reserve1);			/*18-19*/
//printf("%d\t", his->reserve2);			/*20-21*/
//printf("%lf\t", his->time_stamp);			/*22-29*/
//printf("%ld\t", his->maker);				/*30-33*/
//printf("%s\n", his->reserved);			/*34-64*/
//printf("%s\n", his->comment);
//printf("-------------------------------\n");
	
	
	fclose(f);
	return(0);
}

/*----------------------------------------------------------------------*/
int read_his(his,img)
HISHeader	*his;
IMGHeader	*img;
{
	int		j;

// read comment and image data from file1
	if (fread(his, sizeof(char), HIS_Header_Size, fi) != HIS_Header_Size){
		fprintf(stderr, "EOF in input file (header)\n"); free(his->comment); return(-1);
	}
//	if (strncmp(his->head, "IM", 2)){
//		fclose(fi); free(his->comment);
//		fprintf(stderr, "This is not HIS image file.\n"); 
//		return(-11);
//	}

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
//	his->comment[his->comment_length] = '\0';

	
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
	if(his->type==0){
		if ((j = fread(odata, sizeof(unsigned char), NP, fi)) != NP){
			fclose(fi); free(his->comment); 
			fprintf(stderr, " Error reading at %d\n",  j);
			return(-14);
		}
	}

	sprintf(img->head,"IM");
//	img->comment_length=his->comment_length;
	img->comment_length=0;
	img->width=his->width;
	img->height=his->height;
	img->x_offset=his->x_offset;
	img->y_offset=his->y_offset;
//	img->type=his->type;
	img->type=2;
//	sprintf(img->reserved,"");
//	printf("%s\n",his->comment);
//	sprintf(img->comment,"%s",his->comment);

		return(0);
	}

/*----------------------------------------------------------------------*/


int	main(argc,argv)
int		argc;
char	*argv[];
{
	unsigned short	*image1, *outimg;
	long	i, j, jj;
	unsigned long	PN, NL, ppp;
	char	readfile[20], outfile[20];
	char	*cmt;
	short	clg;
	unsigned long	i_end, k, kk;

	HISHeader		his;
	IMGHeader		img;

	if(argc != 4){   // when bad number of parameters
		printf("Usage is 'spl_E hisfile partnum loopnum'\n");
		return(-1);
	}
	sprintf(readfile, ("%s"), argv[1]);
	PN=atol(argv[2]);
	NL=atol(argv[3]);

	if ((j = read_hishead(readfile, &his)) != 0){
		printf("error, read his-head - %ld", j);
		free(his.comment); 
		return(-1);
	}
	i_end=his.n_image1+65536*his.n_image2;
	if(his.type==6) printf("12 bits %ld images\n", i_end); /* ŠÔˆá‚¦‚Ä‚È‚¢‚©H*/
	if(his.type==2) printf("16 bits %ld images\n", i_end);
	if(his.type==0) printf("8  bits %ld images\n", i_end);

	printf("total number of shot is %ld.\n", PN*NL);

	if ((PN*NL)>i_end){
		printf ("error: %ld is less than %ld. Looks strange.\n",i_end, PN*NL);
//		return 1;
	}
	
	printf("reading %s.\n",readfile);

// initialize
	
	if(his.type==2){
		NP=his.width*his.height;
		data=(unsigned short *)malloc(NP*sizeof(unsigned short));
	}
	if(his.type==6){
		NP=(his.width*his.height)*3/2;
		cdata=(unsigned char *)malloc(NP*sizeof(unsigned char));
	}
	if(his.type==0){
		NP=his.width*his.height;
		odata=(unsigned char *)malloc(NP*sizeof(unsigned char));
	}

//open input files
	if((fi = fopen(argv[1],"rb")) == NULL){
		printf("can not open %s for input\n", argv[1]);
		return(-10);
	}
	ppp=0;
	for(k=0;k<NL;++k){
		if (k!=0) ppp=100;
		for(kk=ppp;kk<PN;++kk){
// read input images
			if ((j = read_his(&his, &img)) != 0){
				printf("something wrong -- return value is %ld", i);
				free(his.comment);
				return(-1);
			}

			outimg = (unsigned short *) malloc(img.width * img.height * sizeof(unsigned short));
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
//						printf("%d\t%d\t%d\n",j,jj,*(outimg+j));
						j=j+1;
					}
					if ((jj%3)==2){
						*(outimg+j)=((*(cdata+jj-1)%16)*256)+*(cdata+jj);
//						printf("%d\t%d\t%d\n",j,jj,*(outimg+j));
						j=j+1;
					}
				}
			}
			if(his.type==0){
				for(jj=0;jj<NP;++jj){
					*(outimg+jj)=(unsigned short)*(odata+jj);
				}
			}
//		if(j!=img.width * img.height) printf("\nerror\n");
//		if(j==img.width * img.height) printf("\nok\n");

#ifdef WINDOWS
			if (PN<10)     {sprintf(outfile, ("%03d\\raw\\a%01d.img"), (int)(k+1), (int)(kk+1));}
			if (PN>9)      {sprintf(outfile, ("%03d\\raw\\a%02d.img"), (int)(k+1), (int)(kk+1));}
			if (PN>99)     {sprintf(outfile, ("%03d\\raw\\a%03d.img"), (int)(k+1), (int)(kk+1));}
			if (PN>999)    {sprintf(outfile, ("%03d\\raw\\a%04d.img"), (int)(k+1), (int)(kk+1));}
			if (PN>9999)   {sprintf(outfile, ("%03d\\raw\\a%05d.img"), (int)(k+1), (int)(kk+1));}
			if (PN>99999)  {sprintf(outfile, ("%03d\\raw\\a%06d.img"), (int)(k+1), (int)(kk+1));}
//			if (PN>999999) {sprintf(outfile, ("%03d\\raw\\a%07d.img"), (int)(k+1), (int)(kk+1));}
#else
			if (PN<10)     {sprintf(outfile, ("%03d/raw/a%01d.img"), (int)(k+1), (int)(kk+1));}
			if (PN>9)      {sprintf(outfile, ("%03d/raw/a%02d.img"), (int)(k+1), (int)(kk+1));}
			if (PN>99)     {sprintf(outfile, ("%03d/raw/a%03d.img"), (int)(k+1), (int)(kk+1));}
			if (PN>999)    {sprintf(outfile, ("%03d/raw/a%04d.img"), (int)(k+1), (int)(kk+1));}
			if (PN>9999)   {sprintf(outfile, ("%03d/raw/a%05d.img"), (int)(k+1), (int)(kk+1));}
			if (PN>99999)  {sprintf(outfile, ("%03d/raw/a%06d.img"), (int)(k+1), (int)(kk+1));}
//			if (PN>999999) {sprintf(outfile, ("%03d/raw/a%07d.img"), (int)(k+1), (int)(kk+1));}
#endif
		
// open output file
			if((fo  = fopen(outfile,"wb")) == NULL){
				printf("can not open %s for output\n", outfile);
				return(-1);
			}
// write dummy header
			if (fwrite(&img, sizeof(char), IMG_Header_Size, fo) != IMG_Header_Size){
				fprintf(stderr, "something wrong in %ld (header3)\n", k);
				return(-1);
			}
//			if (fwrite(img.comment, sizeof(char), img.comment_length, fo) != img.comment_length){
//				fprintf(stderr, "something wrong in comment of %d (header)\n", k);
//				return(-1);
//			}

// write image data
			if (i = fwrite(outimg, sizeof(unsigned short), img.width*img.height, fo) > img.width*img.height*sizeof(short)){
				fprintf(stderr, "something wrong in data3 area %ld  %d\n", i,  img.width*img.height);
				return(-1);
			}
// close files
			fclose(fo);
			free(outimg);

			printf("%s\r",outfile);
		}
	}
	printf("\ndone.\n");
	
//finish!
	if(his.type==2){
		free(data);
	}
	if(his.type==6){
		free(cdata);
	}
	if(his.type==0){
		free(odata);
	}
//	free(his.comment);
//	free(img.comment);
	fclose(fi);
	return(0);
}
