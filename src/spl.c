// program spl.c

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*----------------------------------------------------------------------*/

int Error(msg)
char        *msg;
{
	fputs(msg,stderr);
	fputc('\n',stderr);
	return(1);
}

/*----------------------------------------------------------------------*/
int main(argc,argv)
int		argc;
char	**argv;
{
	long	i, j, k, l;
	long	L, M, N;
	

// parameter setting
	if (argc!=3){
		fprintf(stderr, "usage : %s N-shot N-split\n", argv[0]);
		return(1);
	}

	L=atoi(argv[1]);
	M=atoi(argv[2]);
	N=L*M;

	for(j=1;j<M+1;++j){
#ifdef WINDOWS
		printf("mkdir %03d\n", j);
		printf("mkdir %03d\\raw\n", j);
		printf("copy conv.bat %03d\\raw\n", j);
		printf("copy output.log %03d\\raw\n", j);
#else
		printf("mkdir %03d\n", j);
		printf("mkdir %03d/raw\n", j);
		printf("cp conv.bat %03d/raw\n", j);
		printf("cp output.log %03d/raw\n", j);
#endif
	}

	printf("his_spl_E a.his %d %d \n", L, M);

	// append to log file
	FILE		*f;
	if((f = fopen("cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(f,"%s ",argv[i]);
	fprintf(f,"\n");
	fclose(f);

	return 0;
}
