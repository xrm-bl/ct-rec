
#ifndef	Float
#define Float	float
#endif

#ifdef	__cplusplus
#define EXTERN	extern "C"
#else
#define EXTERN	extern
#endif

#ifndef	Filter
#define Filter	Chesler
#else
EXTERN	double	Filter(int i);
#endif

EXTERN Float	**InitCBP(int N,int M);
EXTERN Float	**EndCBP();
EXTERN Float	**CBP(double dr,double r0,double t0);
EXTERN void	PrepareCBP();
EXTERN void	ExecuteCBP(double dr,double r0,double t0);
EXTERN void	BeginCBP(double dr,double r0,double t0);
EXTERN void	TermCBP();

#ifndef	M_PI
#define M_PI	3.14159265358979323846
#endif
