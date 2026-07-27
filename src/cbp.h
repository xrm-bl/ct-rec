
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

/* 打ち切り(カッピング)補正パッドの既定しきい値(比率)。環境変数 PAD_THRESH
   未設定のときに使う。cbp_thread*.c / cbp.cu が参照。0 以下で既定 OFF、
   PAD_THRESH を明示的に 0 (以下) にすると個別に強制 OFF。
   コンパイル時に -DPAD_THRESH_DEFAULT=... で上書き可。 */
#ifndef	PAD_THRESH_DEFAULT
#define PAD_THRESH_DEFAULT	0.3
#endif
