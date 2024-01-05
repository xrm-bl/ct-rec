
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include "cell.h"

#ifdef	CSI_RGB
#include "rif_rgb.h"
#else
#include "rif.h"
#endif

#include "csi.h"

#define PATH_LEN	1023

static void	Error(head,tail)
char		*head,*tail;
{
	(void)fputs(head,stderr); (void)fputs(" : ",stderr);
	(void)fputs(tail,stderr); (void)fputs(".\n",stderr); exit(1);
}

#define ALLOC(type,size)	(type *)malloc(sizeof(type)*(size_t)(size))

static char	Directory[PATH_LEN];

#ifdef	CSI_BS
static PBS	*Alloc(Nz)
#else
static char	**Alloc(Nz)
#endif
int		Nz;
{
#ifdef	CSI_BS
	PBS	*paths=ALLOC(PBS,Nz);
#else
	char	**paths=ALLOC(char *,Nz);
#endif
	if (paths==NULL) Error(Directory,"no memory for paths");

	return paths;
}

static int	DirLen,nx,ny,bps;

static char	*Peep(name)
char		*name;
{
	char	*path=ALLOC(char,DirLen+strlen(name)+1);
	int	Nx,Ny,BPS;

	if (path==NULL) Error(Directory,"no memory for path");

	(void)strcat(strcpy(path,Directory),name);
#ifdef	CSI_RGB
	ReadImageFile_RGB(path,&Nx,&Ny,&BPS,NULL,NULL,NULL,NULL);
#else
	ReadImageFile(path,&Nx,&Ny,&BPS,NULL,NULL);
#endif
	if (nx<Nx) nx=Nx;
	if (ny<Ny) ny=Ny; if (bps<BPS) bps=BPS; return path;
}

static int	Normal(path1,path2)
#ifdef	CSI_BS
PBS		*path1,*path2;
#else
char		**path1,**path2;
#endif
{
#ifdef	CSI_BS
	return strcmp(path1->path+DirLen,path2->path+DirLen);
#else
	return strcmp(*path1+DirLen,*path2+DirLen);
#endif
}

static int	Reverse(path1,path2)
#ifdef	CSI_BS
PBS		*path1,*path2;
#else
char		**path1,**path2;
#endif
{
	return Normal(path2,path1);
}

static int	Numeric(path,d)
char		*path;
double		*d;
{
	for (path+=DirLen; *path!='\0'; path++)
	    if (*path=='-' || *path=='+' ||
	       (*path>='0' && *path<='9')) return sscanf(path,"%lf",d);

	return 0;
}

static int	Ascending(path1,path2)
#ifdef	CSI_BS
PBS		*path1,*path2;
#else
char		**path1,**path2;
#endif
{
	double	d1,d2;
#ifdef	CSI_BS
	int	i1=Numeric(path1->path,&d1),
		i2=Numeric(path2->path,&d2);
#else
	int	i1=Numeric(*path1,&d1),
		i2=Numeric(*path2,&d2);
#endif
	return (i1==1 && i2==1)?(d1>d2)?1:(d1<d2)?-1:0:
	       (i1==1)?1:(i2==1)?-1:Normal(path1,path2);
}

static int	Descending(path1,path2)
#ifdef	CSI_BS
PBS		*path1,*path2;
#else
char		**path1,**path2;
#endif
{
	return Ascending(path2,path1);
}

#ifdef	CSI_BS
PBS	*CheckSliceImages(directory,nameFile,Nx,Ny,Nz,BPS)
#else
char	**CheckSliceImages(directory,nameFile,Nx,Ny,Nz,BPS)
#endif
char	*directory,*nameFile;
int	*Nx,*Ny,*Nz,*BPS;
{
	int	nz,z;
	char	name[PATH_LEN],
		*nifn="no image file name",
		*mifn="missing image file name";
#ifdef	CSI_BS
	PBS	*paths;
#else
	char	**paths;
#endif

	DirLen=strlen(strcat(strcpy(Directory,directory),"/")); nx=ny=bps=0;

	if (nameFile[0]!='-')
{
	FILE	*file;
	char	line[PATH_LEN],
		*tll="too long line";

	if ((file=fopen(nameFile,"r"))==NULL)
	    Error(nameFile,"name file not open");

	for (nz=0; fgets(line,PATH_LEN,file)!=NULL; )
	    if (*line!='\0') {
		if (line[strlen(line)-1]!='\n') Error(nameFile,tll);

		if (sscanf(line,"%s",name)==1) ++nz;
	    }

	if (nz==0) Error(nameFile,nifn);

	paths=Alloc(nz);

	rewind(file);

	for (z=0; z<nz && fgets(line,PATH_LEN,file)!=NULL; )
	    if (*line!='\0') {
		if (line[strlen(line)-1]!='\n') Error(nameFile,tll);
#ifdef	CSI_BS
		switch(sscanf(line,"%s %lf %lf",name,&(paths[z].base),
						     &(paths[z].step))) {
		    case 1 :
			paths[z].base=0.0; paths[z].step=1.0; break;
		    case 2 :
			Error(nameFile,"missing base or step");
		    case 3 :
			break;
		    default :
			continue;
		}
		paths[z++].path=Peep(name);
#else
		if (sscanf(line,"%s",name)==1) paths[z++]=Peep(name);
#endif
	    }

	if (z!=nz) Error(nameFile,mifn);

	(void)fclose(file);
}
	else if (nameFile[1]!='\0' && nameFile[1]!='-')
{
	char	*ptr;

	ptr=nameFile+1;
	for (nz=0; sscanf(ptr,"%s",name)==1; nz++)
	    ptr=strstr(ptr,name)+strlen(name);

	if (nz==0) Error(nameFile+1,nifn);

	paths=Alloc(nz);

	ptr=nameFile+1;
	for (z=0; z<nz && sscanf(ptr,"%s",name)==1; z++) {
#ifdef	CSI_BS
	    paths[z].path=Peep(name); paths[z].base=0.0; paths[z].step=1.0;
#else
	    paths[z]=Peep(name);
#endif
	    ptr=strstr(ptr,name)+strlen(name);
	}
	if (z!=nz) Error(nameFile+1,mifn);
}
	else
{
	DIR		*dir;
	struct dirent	*ent;

	if ((dir=opendir(directory))==NULL)
	    Error(directory,"directory not open");

	for (nz=0; (ent=readdir(dir))!=NULL; )
	    if (strcmp(ent->d_name,".") && strcmp(ent->d_name,"..")) ++nz;

	if (nz==0) Error(directory,nifn);

	paths=Alloc(nz);

	rewinddir(dir);

	for (z=0; z<nz && (ent=readdir(dir))!=NULL; )
	    if (strcmp(ent->d_name,".") && strcmp(ent->d_name,".."))
#ifdef	CSI_BS
	    {
		paths[z].path=Peep(ent->d_name);
		paths[z].base=0.0;
		paths[z].step=1.0;
		++z;
	    }
#else
		paths[z++]=Peep(ent->d_name);
#endif
	if (z!=nz) Error(directory,mifn);

	(void)closedir(dir);

	if (nameFile[1]=='\0')
	    qsort(paths,(size_t)nz,sizeof(*paths),Normal);
	else if (nameFile[2]=='r' && nameFile[3]=='\0')
	    qsort(paths,(size_t)nz,sizeof(*paths),Reverse);
	else if (nameFile[2]=='n' && nameFile[3]=='\0')
	    qsort(paths,(size_t)nz,sizeof(*paths),Ascending);
	else if ((nameFile[2]=='r' && nameFile[3]=='n') ||
		 (nameFile[2]=='n' && nameFile[3]=='r'))
	    qsort(paths,(size_t)nz,sizeof(*paths),Descending);
	else
	    Error(nameFile+1,"unknown sorting option");
}
#ifdef	CSI_SLICE
	for (z=0; z<nz; z++)
#ifdef	CSI_BS
	    paths[z].path+=DirLen;
#else
	    paths[z]+=DirLen;
#endif
#endif
	if (Nx!=NULL) *Nx=nx;
	if (Ny!=NULL) *Ny=ny;
	if (Nz!=NULL) *Nz=nz; if (BPS!=NULL) *BPS=bps; return paths;
}
