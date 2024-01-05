//	program cthist
//	usage 'bmake > output_file'

/*----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/*----------------------------------------------------------------------*/

int main(argc,argv)
int		argc;
char	**argv;
{

	long	i, i_num, i_sta, i_end, l_num;
	short	nshot;

	float	f_ps, f_rc, f1, f2, f3, f4;
	FILE	*fo;

	if (argc!=3){
		printf("usage : chk-rc stk_num calc_pos > bat\n");
		return (-1);
	}

	i_num=atoi(argv[1]);
	i_sta=atoi(argv[2]);
	
	for(i=0;i<i_num;++i){
#ifdef WINDOWS
		printf("cd %03d\\raw\n",i+1);
		printf("call conv.bat\n");
		printf("ct_rec_g_c %d\n", i_sta);
		printf("move rec%05d.tif ..\\..\\%03d.tif\n", i_sta, i+1);
		printf("cd ..\\..\n");
#else
		printf("cd %03d/raw\n",i+1);
		printf("tr -d \"\\15\" < conv.bat |csh \n");
		printf("ct_rec_g_c %d\n", i_sta);
		printf("mv rec%05d.tif ../../%03d.tif\n", i_sta, i+1);
		printf("cd ../..\n");
#endif
		printf("pid %03d.tif >> image-description.txt\n", i+1);
	}

	printf("mkdir rc-check\n");
#ifdef WINDOWS
	printf("move *.tif rc-check\n");
#else
	printf("mv *.tif rc-check\n");
#endif
	

	
// append to log file
	FILE		*f;
	if((f = fopen("cmd-hst.log","a")) == NULL){
		return(-10);
	}
	for(i=0;i<argc;++i) fprintf(f,"%s ",argv[i]);
	fprintf(f,"\n");
	fclose(f);

	return(0);
}
