//	program cthist
//	usage 'sfa layer from to step > output_file'

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <malloc.h>

#define INTEL

long	i, j, n;
long	l_layer, l_from, l_to;
float	f_step, delta;
char	line[50];


/*----------------------------------------------------------------------*/

int main(argc,argv)
int		argc;
char	**argv;
{
	delta = (float)1.00;
	if (argc<4){
		printf("usage : ofsfa layer from to> output_file\n");
		exit(-1);
	}
	l_layer=atoi(argv[1]);
	l_from =atoi(argv[2]);
	l_to   =atoi(argv[3]);

	n=(long)abs(((float)l_from-(float)l_to));

//	printf("%d\t%d\t%d\t%f\t%d\n", l_layer, l_from, l_to, f_step, n);

	printf("mkdir check-rc\n");
	i=0;
	for (j=0; j<n+1; j++){
		sprintf(line, "ofct_srec_g_c raw %d 0 %d-%d 1.0 0.0 check-rc\n",l_from+j, l_layer, l_layer);
		printf("%s", line);
		sprintf(line, "move check-rc\\rec%05d.tif  check-rc\\%05d_%05d.tif\n",l_layer, l_layer, (long)((l_from+j)));
		printf("%s", line);
	}

	exit(0);
}
