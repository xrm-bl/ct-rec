
#include <stdio.h>
#include <stdlib.h>
#include "cell.h"
#include "csi.h"
#include "sif.h"

extern void	Error(),ReadSliceImage();

static int	Unit[6][3]={{-1,0,0},{0,-1,0},{0,0,-1},{1,0,0},{0,1,0},{0,0,1}};

static int	Scan(axis,unit)
char		*axis;
int		**unit;
{
	int	sign;

	switch(axis[0]) {
	    case '-': sign=0; break;
	    case '+': sign=3; break;
	    default : return -1;
	}
	switch(axis[1]) {
	    case 'x': case 'X': *unit=Unit[sign  ]; return 0;
	    case 'y': case 'Y': *unit=Unit[sign+1]; return 1;
	    case 'z': case 'Z': *unit=Unit[sign+2]; return 2;
	}
	return -1;
}

static int	*Uh,*Uv,*Ud,h0,v0,d0,Nh,Nv;

static int	Check(Ah,Av,Ad,line,Nx,Ny,Nz)
char		*Ah,*Av,*Ad,*line;
int		Nx,Ny,Nz;
{
	int	h,v,d,Nd;

	if ((h=Scan(Ah,&Uh))<0 ||
	    (v=Scan(Av,&Uv))<0 || h==v ||
	    (d=Scan(Ad,&Ud))<0 || v==d || d==h)
	    if (line!=NULL) {
		(void)fputs("bad axial direction :\n",stderr); Error(line);
	    }
	    else
		Error("bad axial direction.");

	if ((Nh=Uh[0]*Nx+Uh[1]*Ny+Uh[2]*Nz)>0) h0=0; else h0=(Nh=(-Nh))-1;
	if ((Nv=Uv[0]*Nx+Uv[1]*Ny+Uv[2]*Nz)>0) v0=0; else v0=(Nv=(-Nv))-1;
	if ((Nd=Ud[0]*Nx+Ud[1]*Ny+Ud[2]*Nz)>0) d0=0; else d0=(Nd=(-Nd))-1;

	return Nd;
}

static void	Store(org,d,new,BPS,path)
Cell		***org,**new;
int		d,BPS;
char		*path;
{
	int	v,h;

	for (v=0; v<Nv; v++)
	for (h=0; h<Nh; h++)
	    new[v][h]=org[Uh[2]*(h-h0)+Uv[2]*(v-v0)+Ud[2]*(d-d0)]
			 [Uh[1]*(h-h0)+Uv[1]*(v-v0)+Ud[1]*(d-d0)]
			 [Uh[0]*(h-h0)+Uv[0]*(v-v0)+Ud[0]*(d-d0)];

	StoreImageFile(path,Nh,Nv,BPS,new,NULL);
}

#define LEN	1024

int	main(argc,argv)
int	argc;
char	**argv;
{
	char	**paths,line[LEN],Ah[LEN],Av[LEN],Ad[LEN],path[LEN],form[LEN];
	int	Nx,Ny,Nz,BPS,Nhv,yhv,z,len,d,Nd;
	Cell	***org,**new;
	
	int		i;

	if (argc!=3 && argc!=7 && argc!=8)
    Error("usage : si_rar orgDir nameFile\n"
	  "        si_rar orgDir nameFile hAxis vAxis dAxis newDir\n"
	  "        si_rar orgDir nameFile hAxis vAxis dAxis sliceNo newTIFF");

	paths=CheckSliceImages(argv[1],argv[2],&Nx,&Ny,&Nz,&BPS);

	Nhv=(Nx>Ny)?(Nx>Nz)?Nx:Nz:(Ny>Nz)?Ny:Nz;

#define ALLOC(type,noe)	(type *)malloc(sizeof(type)*(noe))

	if ((  org=ALLOC(Cell **,(size_t)Nz))==NULL ||
	    ( *org=ALLOC(Cell * ,(size_t)Nz*(size_t)Ny))==NULL ||
	    (**org=ALLOC(Cell   ,(size_t)Nz*(size_t)Ny*(size_t)Nx))==NULL ||
	    (  new=ALLOC(Cell * ,(size_t)Nhv))==NULL ||
	    ( *new=ALLOC(Cell   ,(size_t)Nhv*(size_t)Nhv))==NULL)
	    Error("no allocatable memory for cell.");

	for (yhv=1; yhv<Ny; yhv++) org[0][yhv]=org[0][yhv-1]+Nx;

	for (z=1; z<Nz; z++) {
	    org[z]=org[z-1]+Ny;
	    for (yhv=0; yhv<Ny; yhv++) org[z][yhv]=org[z][yhv-1]+Nx;
	}
	for (yhv=1; yhv<Nhv; yhv++) new[yhv]=new[yhv-1]+Nhv;

#define READ()	for (z=0; z<Nz; z++) ReadSliceImage(paths[z],Nx,Ny,org[z])

switch(argc) {
    case 3:
	READ();

	while (fgets(line,LEN,stdin)!=NULL) {
	    for (len=0; line[len]!='\0'; len++) ; /* if (len==0) continue; */

	    if (line[--len]!='\n') {
		(void)fputs("missing end of line :\n",stderr); Error(line);
            }
	    line[len]='\0';
	    switch(sscanf(line,"%s %s %s %d %s",Ah,Av,Ad,&d,path)) {
		case 0 : continue;
		case 5 : break;
		default: (void)fputs("bad parameters :\n",stderr); Error(line);
	    }
	    Nd=Check(Ah,Av,Ad,line,Nx,Ny,Nz);
	    if (d<0 || d>=Nd) {
		(void)fputs("bad slice number :\n",stderr); Error(line);
	    }
	    Store(org,d,new,BPS,path);
	}
	break;
    case 7:
	Nd=Check(argv[3],argv[4],argv[5],(char *)NULL,Nx,Ny,Nz); READ();

	len=1; for (d=10; d<Nd; d*=10) ++len;
	(void)sprintf(form,"%s/%%0%dd.tif",argv[6],len);

	for (d=0; d<Nd; d++) {
	    (void)sprintf(path,form,d); Store(org,d,new,BPS,path);
	}
	break;
    case 8:
	Nd=Check(argv[3],argv[4],argv[5],(char *)NULL,Nx,Ny,Nz);
	if (sscanf(argv[6],"%d",&d)!=1 ||
	    d<0 || d>=Nd) Error("bad slice number.");

	READ(); Store(org,d,new,BPS,argv[7]); break;
}
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
