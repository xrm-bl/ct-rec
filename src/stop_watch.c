
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if	defined(_WIN32)
#include <windows.h>

static double	Time()
{
	union {
		FILETIME	ft;
		ULONGLONG	ull;
	} t;

	GetSystemTimeAsFileTime(&(t.ft));
	return (double)(t.ull-116444736000000000LL)*1.0e-7;
}

#define PRECISION	"7"
#elif	defined(__unix__)
#include <sys/time.h>

static double	Time()
{
	struct timeval	t;

	(void)gettimeofday(&t,NULL);
	return (double)t.tv_sec+(double)t.tv_usec*1.0e-6;
}

#define PRECISION	"6"
#else
#include <time.h>

#define Time()		(double)time(NULL)
#define PRECISION	"0"
#endif

int	main(int argc,char **argv)
{
	size_t	l;
	int	i;
	char	*c;
	double	t;

	if (argc<2) {
	    (void)fputs("usage : stop_watch command\n",stderr); exit(1);
	}
	l=0; for (i=1; i<argc; i++) l+=(strlen(argv[i])+1);

	if ((c=(char *)malloc(l))==NULL) {
	    (void)fputs("memory allocation erorr.\n",stderr); exit(1);
	}
	(void)strcpy(c,argv[1]);
	for (i=2; i<argc; i++) (void)strcat(strcat(c," "),argv[i]);

	t=Time(); i=system(c);
	(void)fprintf(stderr,"%0." PRECISION "lf\n",Time()-t); return i;
}
