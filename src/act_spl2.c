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

void replace_in_file(const char *filename, const char *search, const char *replace) {
	FILE *file = fopen(filename, "r");
	if (!file) {
		perror("ファイルを開けませんでした");
		return;
	}

	char temp_filename[] = "temp.txt";
	FILE *temp_file = fopen(temp_filename, "w");
	if (!temp_file) {
		perror("一時ファイルを作成できませんでした");
		fclose(file);
		return;
	}

	char buffer[1024];
	while (fgets(buffer, sizeof(buffer), file)) {
		char *pos;
		while ((pos = strstr(buffer, search)) != NULL) {
			*pos = '\0';
			fprintf(temp_file, "%s%s", buffer, replace);
			strcpy(buffer, pos + strlen(search));
		}
		fprintf(temp_file, "%s", buffer);
	}

	fclose(file);
	fclose(temp_file);

	remove(filename);
	rename(temp_filename, filename);
}

/*----------------------------------------------------------------------*/
int main(argc,argv)
int		argc;
char	**argv;
{
	int	i, j, l, mk;
	float	gk;
	int	L, M, N;
	long	KK, k;

	char    command[50];

// parameter setting
	gk=0;
	if (argc<3){
		fprintf(stderr, "usage : %s N-shot N-split (m_kernel_size) (g_kernel_size) \n", argv[0]);
		return(1);
	}else if (argc==3){
		mk=0;
		gk=0.0;
	}else if (argc==4){
		mk=atoi(argv[3]);
		gk=0.0;
	}else if (argc==5){
		mk=atoi(argv[3]);
		gk=atof(argv[4]);
	}

	L=atoi(argv[1]);
	M=atoi(argv[2]);
	N=L*M;

#ifdef WINDOWS
	replace_in_file("conv.bat", "img", "tif");
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

	for(j=1;j<M+1;++j){
		#ifdef WINDOWS
		sprintf(command, "mkdir %03d", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
		sprintf(command, "mkdir %03d\\raw", j);printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
		sprintf(command, "copy conv.bat %03d\\raw", j);printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
		sprintf(command, "copy output.log %03d\\raw", j);printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
//		sprintf(command, "copy lastangle.dat %03d\\raw", j);printf("%s\n"command);
//		if (system(command) == -1) {printf("command error at %d\n",j); }
		#else
		sprintf(command, "mkdir %03d", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
		sprintf(command, "mkdir %03d/raw", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
		sprintf(command, "cp conv.bat %03d/raw", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
		sprintf(command, "cp output.log %03d/raw", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
//		sprintf(command, "cp lastangle.dat %03d/raw", j); printf("%s\n",command);
//		if (system(command) == -1) {printf("command error at %d\n",j);}
		#endif
	}

	k=0;
	for(j=1;j<M+1;++j){
		for(l=1;l<L+1;++l){
			k=k+1;
			#ifdef WINDOWS
			if (L>99)   {sprintf(command, "start /b tif_mgf a%06ld.tif %03d\\raw\\a%03d.tif %d %f", k,j,l,mk,gk);}
			if (L>999)  {sprintf(command, "start /b tif_mgf a%06ld.tif %03d\\raw\\a%04d.tif %d %f", k,j,l,mk,gk);}
			if (L>9999) {sprintf(command, "start /b tif_mgf a%06ld.tif %03d\\raw\\a%05d.tif %d %f", k,j,l,mk,gk);}
			if (system(command) == -1) {printf("command error at tif %d\n",k); }
			#else
			if (L>99)   {sprintf(command, "tif_mgf a%06ld.tif %03d/raw/a%03d.tif %d %f", k,j,l,mk,gk);}
			if (L>999)  {sprintf(command, "tif_mgf a%06ld.tif %03d/raw/a%04d.tif %d %f", k,j,l,mk,gk);}
			if (L>9999) {sprintf(command, "tif_mgf a%06ld.tif %03d/raw/a%05d.tif %d %f", k,j,l,mk,gk);}
			if (system(command) == -1) {printf("command error at tif %d\n",k); }
			#endif
			printf("%s\r",command);
		}
	}
//printf("his_spl_E a.his %d %d \n", L, M);
printf("\n");
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
