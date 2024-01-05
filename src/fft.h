
typedef struct {
		double	r,i;
	} Complex;

extern void	FFT(int S,int L,Complex *E,Complex *C);
extern int	SetUpFFT(int N,Complex **E);
