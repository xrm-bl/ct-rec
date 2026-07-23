// program spl.c

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <string.h>

#ifndef WINDOWS
#include <unistd.h>
#include <sys/wait.h>

/*----------------------------------------------------------------------*/
/* 背景ジョブ実行(上限付き): Windows 側の "start /b" に相当。
   Linux では system() が完了までブロックするため tif_mgf が逐次実行に
   なっていた。fork+exec で背景実行し、同時実行数を ACT_SPL_JOBS
   (既定 8) で制限、終了前に wait_bg() で全ジョブを待ち合わせる。 */

static int	bg_njobs=0, bg_maxjobs=0;

static void	run_bg(const char *cmd)
{
	pid_t	pid;

	if (bg_maxjobs<=0) {
	    char *e=getenv("ACT_SPL_JOBS");
	    bg_maxjobs=(e!=NULL && atoi(e)>0)?atoi(e):8;
	}
	while (bg_njobs>=bg_maxjobs) {
	    if (wait(NULL)>0) bg_njobs--;
	    else { bg_njobs=0; break; }
	}
	if ((pid=fork())==0) {
	    execl("/bin/sh","sh","-c",cmd,(char *)NULL);
	    _exit(127);
	}
	if (pid>0) bg_njobs++;
	else {	/* fork 失敗時は従来どおり同期実行 */
	    printf("fork error: run synchronously\n");
	    if (system(cmd) == -1) printf("command error (sync)\n");
	}
}

static void	wait_bg(void)
{
	while (bg_njobs>0) {
	    if (wait(NULL)<0) break;
	    bg_njobs--;
	}
}
#endif

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
			memmove(buffer, pos + strlen(search), strlen(pos + strlen(search)) + 1);
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

	char    command[512];

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
	snprintf(command, sizeof(command), "sed -i -e 's/\r//g' conv.bat");				// CRLF -> LF
	if (system(command) == -1) {printf("command error at conv.\n");}
	snprintf(command, sizeof(command), "sed -i -e 's/img/tif/g' conv.bat");			// img -> tif
	if (system(command) == -1) {printf("command error at conv.\n");}
	snprintf(command, sizeof(command), "sed -i -e 's/ren/mv/g' conv.bat");			// ren -> mv
	if (system(command) == -1) {printf("command error at conv.\n");}
	snprintf(command, sizeof(command), "sed -i -e 's/copy/cp/g' conv.bat");			// copy -> cp
	if (system(command) == -1) {printf("command error at conv.\n");}

	snprintf(command, sizeof(command), "tail -n 1 output.log | cut -c28-35 > lastangle.dat");
	if (system(command) == -1) {printf("command error at lastangle.\n");}
	
#endif

	for(j=1;j<M+1;++j){
		#ifdef WINDOWS
		snprintf(command, sizeof(command), "mkdir %03d", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
		snprintf(command, sizeof(command), "mkdir %03d\\raw", j);printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
		snprintf(command, sizeof(command), "copy conv.bat %03d\\raw", j);printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
		snprintf(command, sizeof(command), "copy output.log %03d\\raw", j);printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j); }
//		snprintf(command, sizeof(command), "copy lastangle.dat %03d\\raw", j);printf("%s\n"command);
//		if (system(command) == -1) {printf("command error at %d\n",j); }
		#else
		snprintf(command, sizeof(command), "mkdir %03d", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
		snprintf(command, sizeof(command), "mkdir %03d/raw", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
		snprintf(command, sizeof(command), "cp conv.bat %03d/raw", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
		snprintf(command, sizeof(command), "cp output.log %03d/raw", j); printf("%s\n",command);
		if (system(command) == -1) {printf("command error at %d\n",j);}
//		snprintf(command, sizeof(command), "cp lastangle.dat %03d/raw", j); printf("%s\n",command);
//		if (system(command) == -1) {printf("command error at %d\n",j);}
		#endif
	}

	k=0;
	for(j=1;j<M+1;++j){
		for(l=1;l<L+1;++l){
			int wd=(L>9999)?5:(L>999)?4:3;	/* 出力名の桁数 (L<=99 も3桁: 従来は command が未設定で直前のコマンドを再実行していた) */

			k=k+1;
			#ifdef WINDOWS
			snprintf(command, sizeof(command), "start /b tif_mgf a%06ld.tif %03d\\raw\\a%0*d.tif %d %f", k,j,wd,l,mk,gk);
			if (system(command) == -1) {printf("command error at tif %d\n",k); }
			#else
			/* Windows の "start /b" と同様に背景実行 (上限 ACT_SPL_JOBS, 既定8) */
			snprintf(command, sizeof(command), "tif_mgf a%06ld.tif %03d/raw/a%0*d.tif %d %f", k,j,wd,l,mk,gk);
			run_bg(command);
			#endif
			printf("%s\r",command);
		}
	}
	#ifndef WINDOWS
	wait_bg();	/* 全 tif_mgf の完了を待ってから終了 (後段処理の前提を保証) */
	#endif
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
