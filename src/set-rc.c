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

	if (argc!=2){
		printf("usage : set-rc pixel_size \n");
		return (-1);
	}

	fo=fopen("image-description.txt","r");
	f_ps=atof(argv[1]);
	
	i_num=0;
    while (fscanf(fo, "%f\t%f\t%hd\t%f\t%f\t%f", &f1, &f_rc, &nshot, &f2, &f3, &f4) != EOF) {
    	i_num=i_num+1;
    	printf("cd %03d\n", i_num);
    	printf("mkdir rec\n");
        printf("hp_tg_g_c raw %f %f 0.0 rec\n", f_ps, f_rc);
    	printf("cd ..\n");
    }

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
