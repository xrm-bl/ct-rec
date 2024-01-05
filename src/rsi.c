
#include <stdio.h>
#include <stdlib.h>
#include "cell.h"
#include "rif.h"

extern void	Error();

#define LEN	1024

void	ReadSliceImage(path,Nx,Ny,slice)
char	*path;
int	Nx,Ny;
Cell	**slice;
{
	int	X,Y,y,x;
	Cell	**cell;
	char	msg[LEN];

	ReadImageFile(path,&X,&Y,NULL,&cell,NULL);

	if (X>Nx || Y>Ny) {
	    (void)sprintf(msg,"%s : bad image size.",path); Error(msg);
	}
	for (y=0; y<Y; y++) {
	    for (x=0; x<X; x++) slice[y][x]=cell[y][x];

	    for (; x<Nx; x++) slice[y][x]=0;

	    free(cell[y]);
	}
	for (; y<Ny; y++) for (x=0; x<Nx; x++) slice[y][x]=0;

	free(cell);
}
