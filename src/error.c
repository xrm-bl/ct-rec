
#include <stdio.h>
#include <stdlib.h>

void	Error(msg)
char	*msg;
{
	fputs(msg,stderr); fputc('\n',stderr); exit(1);
}
