
#include <stdio.h>
#include "cell.h"
#include "rif.h"

static void	Error(msg)
char		*msg;
{
	fputs(msg,stderr); fputc('\n',stderr); exit(1);
}

static FILE	*hp;

#define Put1(data)	fputc(data,hp)
#define Put2(data)	Put1(data&0xff); Put1((data>>8)&0xff)

static void	PutCell(x,y,cell)
int		x,y;
Cell		cell;
{
	Put2((int)cell);
}

int	main(argc,argv)
int	argc;
char	**argv;
{
	int	Nx,Ny,len,idx;
	char	*desc;

	if (argc!=3) Error("usage : tiff2hp TIFF HiPic");

	ReadImageFile(argv[1],&Nx,&Ny,NULL,NULL,&desc);

	if (desc==NULL)
	    len=0;
	else
	    for (len=0; desc[len]!='\0'; len++) ;

	if (argv[2][0]=='-')
	    hp=stdout;
	else
	    if ((hp=fopen(argv[2],"wb"))==NULL) Error("file not stored.");

	Put1('I');	/* magic number */
	Put1('M');
	Put2(len);	/* comment length */
	Put2(Nx);	/* image size */
	Put2(Ny);
	Put2(0);	/* image offset */
	Put2(0);
	Put2(2);	/* image type */

	for (idx=14; idx<64; idx++) Put1(0);		/* reserved data */
	for (idx=0; idx<len; idx++) Put1(desc[idx]);	/* comment */

#ifdef	FUNCTION_CELL
	ReadImageFile(argv[1],NULL,NULL,NULL,PutCell,NULL);
#else
{
	Cell	**cell;
	int	y,x;

	ReadImageFile(argv[1],NULL,NULL,NULL,&cell,NULL);

	for (y=0; y<Ny; y++) for (x=0; x<Nx; x++) PutCell(x,y,cell[y][x]);
}
#endif
	if (hp!=stdout) fclose(hp); return 0;
}
