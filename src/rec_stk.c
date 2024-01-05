//	program cthist
//	usage 'bmake > output_file'

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <malloc.h>

#define INTEL

#ifndef	M_PI
#define M_PI	3.14159265358979323846
#endif

char	line[100];

/*----------------------------------------------------------------------*/

int main(argc,argv)
int		argc;
char	**argv;
{

	long	i, i_num, i_sta, i_end, l_num;
	long	j, k, l, ix, iy, NN, nshot;


	if (argc!=4){
		printf("usage : rec_stk num_stack start end \n");
		return (-1);
	}

	i_num=atoi(argv[1]);
	i_sta=atoi(argv[2]);
	i_end=atoi(argv[3]);
	
	l_num=i_end-i_sta+1;

#ifdef WINDOWS
	sprintf(line, "mkdir whole\n");
	printf("%s", line);
	sprintf(line, "mkdir whole\\rec\n");
	printf("%s", line);

	for (j=0;j<i_num;j++){
		for(i=i_sta;i<i_end+1;i++){
			sprintf(line, "move %03d\\rec\\rec%05d.tif whole\\rec\\rec%05d.tif\n",(int)(j+1), (int)i, (int)(i+j*l_num));
			printf("%s", line);
		}
	}
#else
	sprintf(line, "mkdir whole\n");
	printf("%s", line);
	sprintf(line, "mkdir whole/rec\n");
	printf("%s", line);

	for (j=0;j<i_num;j++){
		for(i=i_sta;i<i_end+1;i++){
			sprintf(line, "mv %03d/rec/rec%05d.tif whole/rec/rec%05d.tif\n",(int)(j+1), (int)i, (int)(i+j*l_num));
			printf("%s", line);
		}
	}
#endif

	return(0);
}
