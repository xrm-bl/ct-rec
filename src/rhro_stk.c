/*
 * rhro_stk - renumber the 16/8-bit slices of a Z scan into one stack
 *            (rec_stk for the rh / ro output of act_tg and ofact_srec)
 *
 *   rhro_stk num_stack start end {name ...} > ddd.sh      (Windows: > ddd.bat)
 *
 * act_tg / ofact_srec write every scan's slices as <scan>/rh/rh%05d.tif
 * (16 bit) and <scan>/ro/ro%05d.tif (8 bit), where <scan> is 001, 002, ...
 * This program prints, like rec_stk, the mkdir and mv (move) commands that
 * collect layers start..end of every scan into whole/rh and whole/ro with
 * one continuous numbering: layer i of scan j becomes i + j*(end-start+1).
 * Without names the two sets rh and ro are done; otherwise only the names
 * given (the directory name is also the file prefix, as act_tg writes it).
 * Nothing is moved by this program itself: run the generated script.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int	main(int argc,char **argv)
{
	const char	*dflt[2]={"rh","ro"},**names;
	int	nnames,i,j,k,i_num,i_sta,i_end,l_num;

	if (argc<4) {
	    (void)fprintf(stderr,
		"usage : rhro_stk num_stack start end {name ...} > ddd.sh\n"
		"        prints the commands that renumber <scan>/<name>/<name>%%05d.tif,\n"
		"        layers start..end of scans 001..num_stack, into whole/<name>/\n"
		"        (default names rh ro, the output of act_tg / ofact_srec)\n");
	    return 1;
	}
	i_num=atoi(argv[1]); i_sta=atoi(argv[2]); i_end=atoi(argv[3]);
	if (i_num<1 || i_sta<0 || i_end<i_sta) {
	    (void)fputs("rhro_stk : bad num_stack / start / end.\n",stderr);
	    return 1;
	}
	l_num=i_end-i_sta+1;
	if (argc>4) { names=(const char **)(argv+4); nnames=argc-4; }
	else	    { names=dflt; nnames=2; }

#ifdef WINDOWS
	(void)printf("mkdir whole\n");
	for (k=0; k<nnames; k++) (void)printf("mkdir whole\\%s\n",names[k]);
	for (k=0; k<nnames; k++)
	for (j=0; j<i_num; j++)
	for (i=i_sta; i<=i_end; i++)
	    (void)printf("move %03d\\%s\\%s%05d.tif whole\\%s\\%s%05d.tif\n",
			 j+1,names[k],names[k],i,names[k],names[k],i+j*l_num);
#else
	(void)printf("mkdir whole\n");
	for (k=0; k<nnames; k++) (void)printf("mkdir whole/%s\n",names[k]);
	for (k=0; k<nnames; k++)
	for (j=0; j<i_num; j++)
	for (i=i_sta; i<=i_end; i++)
	    (void)printf("mv %03d/%s/%s%05d.tif whole/%s/%s%05d.tif\n",
			 j+1,names[k],names[k],i,names[k],names[k],i+j*l_num);
#endif

	{   /* record command to cmd-hst.log, as rec_stk does */
	    FILE	*f;
	    if ((f=fopen("cmd-hst.log","a"))!=NULL) {
		for (i=0; i<argc; i++) (void)fprintf(f,"%s ",argv[i]);
		(void)fputc('\n',f); (void)fclose(f);
	    }
	}
	return 0;
}
