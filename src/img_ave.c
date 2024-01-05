/* average of two files                  */
/* usage is 'img_ave file1 file 2 file3' */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#define INTEL
#define HiPic_Header_Size 64

struct HiPic_Header
{
	char					head[2];
	unsigned short			comment_length;
	unsigned short			width;
	unsigned short			height;
	unsigned short			x_offset;
	unsigned short			y_offset;
	unsigned short			type;
	char					reserved[63 - 13];
	char					*comment;
};
typedef struct HiPic_Header Header;


unsigned short	*data;

static FILE *fi1, *fo;

/*----------------------------------------------------------------------*/
int read_hipic(fname, h)
char	*fname;
Header	*h;
{
		short		x_size, y_size;
		int		j;
		FILE	*fi;

//open input files
		if((fi = fopen(fname,"rb")) == NULL)
		{
			printf("can not open %s for input\n", fname);
			free(h->comment); 
			return(-10);
		}
// read comment and image data from file1
		if (fread(h, sizeof(char), HiPic_Header_Size, fi) != HiPic_Header_Size)
		{
			fprintf(stderr, "EOF in %s (header)\n", fname);
			free(h->comment); 
			return(-1);
		}
		if (strncmp(h->head, "IM", 2))
		{
			fclose(fi); 
			free(h->comment); 
			fprintf(stderr, "This is not HiPIc image file.\n"); 
			fprintf(stderr, "File name = %s\n", fname); 
			return(-11);
		}

		h->comment = (char *) malloc(h->comment_length + 1);
		if (h->comment == NULL)
		{
			fclose(fi); 
			free(h->comment); 
			fprintf(stderr, "Connot allocate memory (HiPic_Header.comment)\n"); 
			fprintf(stderr, "File name = %s\n", fname); 
			return(-12);
		}
		if (fread(h->comment, sizeof(char), h->comment_length, fi) != h->comment_length)
		{ 
			fclose(fi); 
			free(h->comment); 
			fprintf(stderr, "EOF found (comment)\n"); 
			fprintf(stderr, "File name = %s\n", fname); 
			return(-13);
		}

		h->comment[h->comment_length] = '\0';

		x_size = h->width;
		y_size = h->height;

		data = (unsigned short *) malloc(x_size * y_size * sizeof(unsigned short));
		if ((j = fread(data, sizeof(unsigned short), x_size*y_size, fi)) != x_size*y_size)
		{
			fclose(fi);
			free(h->comment); 
			fprintf(stderr, " Error reading %s at %d\n", fname, j);
			fprintf(stderr, "File name = %s\n", fname); 
			return(1);
		}

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

		fclose(fi);
		return(0);
	}

/*----------------------------------------------------------------------*/


int	main(argc,argv)
int		argc;
char	*argv[];
{
	unsigned short		*image1, *outimg;
	int			x_size1, y_size1, x_size2, y_size2;
	char		readfile[20];
	long		i, j, jj;
	unsigned long 		*sumdata;

	Header		h;

		if(argc < 3)   // when bad number of parameters
			{ printf("Usage is 'img_ave file1 file 2... output'\n"); exit(-1);}


// read input images
	for(i=1;i<argc-1;++i){
		sprintf(readfile, ("%s"), argv[i]);
		printf("read %s\n",readfile);
		if ((j = read_hipic(argv[i], &h)) != 0){
			printf("something wrong -- return value is %ld", i);
			exit(-1);
		}
		if (i==1){
			sumdata = (unsigned long *) malloc(h.width*h.height*sizeof(unsigned long));
			for(jj=0;jj<h.width*h.height;++jj) *(sumdata+jj)=0;
		}
		for(jj=0;jj<h.width*h.height;++jj){
			*(sumdata+jj)=*(sumdata+jj)+*(data+jj);
		}
		free(data);
	}

	outimg = (unsigned short *) malloc(h.width * h.height * sizeof(unsigned short));
	for(jj=0;jj<h.width * h.height;++jj) *(outimg+jj)=*(sumdata+jj)/(argc-2);

// open output file
		if((fo  = fopen(argv[argc-1],"wb")) == NULL){ // file3 open
			printf("can not open file3 for output\n"); exit(-1);
		}
// write dummy header to file3
		if (fwrite(&h, sizeof(char), HiPic_Header_Size, fo) != HiPic_Header_Size){
			fprintf(stderr, "something wrong in %s (header3)\n", argv[argc-1]);
			exit(-1);
		}
		if (fwrite(h.comment, sizeof(char), h.comment_length, fo) != h.comment_length){
			fprintf(stderr, "something wrong in comment3 of %s (header)\n", argv[argc-1]);
			exit(-1);
		}

// write comment and dark image to file3
		if (i = fwrite(outimg , sizeof(unsigned short), h.width*h.height, fo) > h.width*h.height*sizeof(unsigned short)){
			fprintf(stderr, "something wrong in data3 area %ld  %d\n", i,  h.width*h.height);
			exit(-1);
		}
// close files
		fclose(fo);

		printf("output %s \n",argv[argc-1]);

//finish!
	free(outimg);
	free(sumdata);
	free(h.comment);
	return(0);
}
//  printf("-------------------------------\n");
//  printf("%s\n", HiPic_Header.head);
//  printf("%d\n", HiPic_Header.comment_length);
//  printf("%d\n", HiPic_Header.width);
//  printf("%d\n", HiPic_Header.height);
//  printf("%d\n", HiPic_Header.x_offset);
//  printf("%d\n", HiPic_Header.y_offset);
//  printf("%d\n", HiPic_Header.type);
//  printf("%s\n", HiPic_Header.reserved);
//  printf("%s\n", HiPic_Header.comment);
//  printf("-------------------------------\n");
