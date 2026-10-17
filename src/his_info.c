/*
 * his_info - what is in a.his + conv.bat: frame size, number of frames,
 *            frames per scan and number of scans (Z scan), and the
 *            centre layer.  For the rc-check-his scripts.
 *
 *   his_info {dir}
 *
 * Prints one line on stdout, six integers separated by blanks:
 *
 *   Nx Ny frames L M layer
 *
 *   Nx, Ny  frame size            frames  frames in the his file
 *   L       frames per scan: the highest frame number conv.bat names
 *   M       number of scans: frames / L when that is a whole number > 1
 *           (a Z scan, RHIS_SCAN = 1..M), otherwise 1
 *   layer   Ny / 2, the centre row of the frame (the layer rc-check-his
 *           reconstructs for the axis check unless one is given)
 *
 * A readable summary goes to stderr.  The his file is the one named on
 * the his2img line of conv.bat (or RHIS_H), the conv.bat name is RHIS_C;
 * only the frame headers are read, so this takes well under a second
 * even for a 35 GB file.  Exits with status 1 on any error.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rhis.h"

int	main(int argc,char **argv)
{
	const char	*dir=".";
	int	Nx,Ny,nframes,L,M;

	if (argc>2 || (argc==2 && (!strcmp(argv[1],"-h") || !strcmp(argv[1],"--help")))) {
	    (void)fprintf(stderr,"usage : %s {dir}\n"
			  "        prints \"Nx Ny frames L M layer\" of <dir>/a.his + conv.bat:\n"
			  "        L = frames per scan (highest frame conv.bat names), M = scans in\n"
			  "        the file (frames / L, 1 for a single scan), layer = Ny/2.\n",argv[0]);
	    return 1;
	}
	if (argc==2) dir=argv[1];

	HisLayout(dir,&Nx,&Ny,&nframes,&L,&M);

	if (L<=0)
	    (void)fprintf(stderr,"%s: warning: conv.bat names no frame; frames per scan unknown.\n",argv[0]);
	else if (nframes%L!=0)
	    (void)fprintf(stderr,"%s: warning: %d frames is not a whole number of scans of %d; "
			  "treated as one scan with %d extra frame%s.\n",
			  argv[0],nframes,L,nframes%L,(nframes%L==1)?"":"s");
	(void)fprintf(stderr,"%s: %d x %d, %d frames; %d frame%s per scan (conv.bat), %d scan%s; centre layer %d\n",
		      dir,Nx,Ny,nframes,L,(L==1)?"":"s",M,(M==1)?"":"s",Ny/2);
	(void)printf("%d %d %d %d %d %d\n",Nx,Ny,nframes,L,M,Ny/2);
	return 0;
}
