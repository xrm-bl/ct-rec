// program spl.c

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
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
	int	i, j, l;
	int	L, M, N;
	long	KK, k;

	char    command[50];

// parameter setting
	if (argc!=3){
		fprintf(stderr, "usage : %s N-shot N-split\n", argv[0]);
		return(1);
	}

	L=atoi(argv[1]);
	M=atoi(argv[2]);
	N=L*M;

#ifdef WINDOWS
//	sprintf(command, "sed -i -e 's/img/tif/g' conv.bat");
//	if (system(command) == -1) {printf("command error at conv.\n");}
#else
	sprintf(command, "sed -i -e 's/\r//g' conv.bat");				// CRLF -> LF
	if (system(command) == -1) {printf("command error at conv.\n");}
	sprintf(command, "sed -i -e 's/img/tif/g' conv.bat");			// img -> tif
	if (system(command) == -1) {printf("command error at conv.\n");}
	sprintf(command, "sed -i -e 's/ren/mv/g' conv.bat");			// ren -> mv
	if (system(command) == -1) {printf("command error at conv.\n");}
	sprintf(command, "sed -i -e 's/copy/cp/g' conv.bat");			// copy -> cp
	if (system(command) == -1) {printf("command error at conv.\n");}

	sprintf(command, "tail -n 1 output.log | cut -c28-35 > lastangle.dat");
	if (system(command) == -1) {printf("command error at lastangle.\n");}
	
#endif


	k=0;
	for(j=1;j<M+1;++j){
#ifdef WINDOWS
		sprintf(command, "mkdir %03d", j); printf("%s\n"command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
		sprintf(command, "mkdir %03d\\raw", j);printf("%s\n"command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
		sprintf(command, "copy conv.bat %03d\\raw", j);printf("%s\n"command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
		sprintf(command, "copy output.log %03d\\raw", j);printf("%s\n"command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
//		sprintf(command, "copy lastangle.dat %03d\\raw", j);printf("%s\n"command);
//		if (system(command) == -1) {printf("command error at %d\n",j); }

		for(l=1;l<L+1;++l){
			k=k+1;
			if (L>99)   {sprintf(command, "move a%06ld.tif %03d\\raw\\a%03d.tif", k,j,l);}
			if (L>999)  {sprintf(command, "move a%06ld.tif %03d\\raw\\a%04d.tif", k,j,l);}
			if (L>9999) {sprintf(command, "move a%06ld.tif %03d\\raw\\a%05d.tif", k,j,l);}
			if (system(command) == -1) {printf("command error at tif %d\n",k); }
//			printf("%s\n",command);
		}

#else
		sprintf(command, "mkdir %03d", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
		sprintf(command, "mkdir %03d/raw", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
		sprintf(command, "cp conv.bat %03d/raw", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
		sprintf(command, "cp output.log %03d/raw", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
		sprintf(command, "cp lastangle.dat %03d/raw", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}

		for(l=1;l<L+1;++l){
			k=k+1;
			if (L>99)   {sprintf(command, "mv a%06ld.tif %03d/raw/a%03d.tif", k,j,l);}
			if (L>999)  {sprintf(command, "mv a%06ld.tif %03d/raw/a%04d.tif", k,j,l);}
			if (L>9999) {sprintf(command, "mv a%06ld.tif %03d/raw/a%05d.tif", k,j,l);}
			if (system(command) == -1) {printf("command error at tif %d\n",k); }
//			printf("%s\n",command);
		}

#endif
	}

//	printf("his_spl_E a.his %d %d \n", L, M);
//

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
