/*
 * his_info - what is in a.his + conv.bat: frame size, number of frames,
 *            frames per scan and number of scans (Z scan), and the
 *            centre layer.  For the rc-check-his scripts.
 *
 *   his_info {dir}
 *
 * Prints one line on stdout, seven integers separated by blanks:
 *
 *   Nx Ny frames L M layer kind
 *
 *   Nx, Ny  frame size            frames  frames in the his file
 *   L       frames per scan: the highest frame number conv.bat names
 *   M       number of scans: frames / L when that is a whole number > 1
 *           (a Z scan, RHIS_SCAN = 1..M), otherwise 1
 *   layer   Ny / 2, the centre row of the frame (the layer rc-check-his
 *           reconstructs for the axis check unless one is given)
 *   kind    180 or 360: the sample views of output.log reach 350 degrees
 *           (the test the his_* programs apply) -> 360-degree offset scan.
 *           rc-check-his needs it: for an offset scan the ImageDescription
 *           of his_rec carries the image centre dcnt+0.5 while his_srec
 *           wants the detector axis position Rc = dcnt, so 0.5 is taken
 *           off before the value goes into center.log.
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
	int	Nx,Ny,nframes,L,M,kind;

	if (argc>2 || (argc==2 && (!strcmp(argv[1],"-h") || !strcmp(argv[1],"--help")))) {
	    (void)fprintf(stderr,"usage : %s {dir}\n"
			  "        prints \"Nx Ny frames L M layer kind\" of <dir>/a.his + conv.bat +\n"
			  "        output.log: L = frames per scan (highest frame conv.bat names),\n"
			  "        M = scans in the file (frames / L, 1 for a single scan),\n"
			  "        layer = Ny/2, kind = 180 or 360 (offset scan).\n",argv[0]);
	    return 1;
	}
	if (argc==2) dir=argv[1];

	HisLayout(dir,&Nx,&Ny,&nframes,&L,&M);
	kind=(HisMaxSampleAngle(dir)>=350.0)?360:180;

	if (L<=0)
	    (void)fprintf(stderr,"%s: warning: conv.bat names no frame; frames per scan unknown.\n",argv[0]);
	else if (nframes%L!=0)
	    (void)fprintf(stderr,"%s: warning: %d frames is not a whole number of scans of %d; "
			  "treated as one scan with %d extra frame%s.\n",
			  argv[0],nframes,L,nframes%L,(nframes%L==1)?"":"s");
	(void)fprintf(stderr,"%s: %d x %d, %d frames; %d frame%s per scan (conv.bat), %d scan%s; centre layer %d; %d-degree scan\n",
		      dir,Nx,Ny,nframes,L,(L==1)?"":"s",M,(M==1)?"":"s",Ny/2,kind);
	(void)printf("%d %d %d %d %d %d %d\n",Nx,Ny,nframes,L,M,Ny/2,kind);
	return 0;
}
