
#include <stddef.h>

extern void	Error();

void	RangeList(string,limit,on_off)
char	*string,*on_off;
size_t	limit;
{
	size_t	n,n1,n2;
	int	c;
	char	*brl="bad range list.";

	for (n=0; n<=limit; n++) on_off[n]=0;

#define GETC()	(*(string++))
#define STON(n)	for (n=0; c>='0' && c<='9'; c=GETC()) \
		    if ((n=n*10+c-'0')>limit) Error(brl)

	do {
	    while ((c=GETC())==',') ;

	    if (c!='\0') {
		if (c=='-') {
		    n1=0;
		    if ((c=GETC())==',' || c=='\0')
			n2=limit;
		    else {
			STON(n2);
			if (c!=',' && c!='\0') Error(brl);
		    }
		}
		else {
		    STON(n1);
		    if (c==',' || c=='\0')
			n2=n1;
		    else if (c=='-')
			if ((c=GETC())==',' || c=='\0')
			    n2=limit;
			else {
			    STON(n2);
			    if ((c!=',' && c!='\0') || n2<n1) Error(brl);
			}
		    else
			Error(brl);
		}
		for (n=n1; n<=n2; n++) on_off[n]=1;
	    }
	} while (c==',');
}
