/*----------------------------------------------------------------------*/
/* paganin.h - single-distance phase retrieval for a homogeneous object */
/*                                                                      */
/*   D. Paganin, S. C. Mayo, T. E. Gureyev, P. R. Miller, S. W. Wilkins, */
/*   J. Microsc. 206, 33-40 (2002), Eq. (10), parallel beam (M = 1):    */
/*                                                                      */
/*     T(r) = -(1/mu) ln( F^-1[ mu F{ I(r,R2)/I_in } / (R2 d |k|^2 + mu) ] ) */
/*                                                                      */
/*   |k| is the angular wavenumber (2 pi / length); on a grid of pitch  */
/*   W and padded size N, k = 2 pi m / (W N), m = -N/2 .. N/2-1.        */
/*                                                                      */
/* What this module returns is the filtered transmission                */
/*     T'(r) = F^-1[ F{I/I_in} / (1 + (R2 d/mu) |k|^2) ]                 */
/* so that the caller can form the projection value mu*T = -ln T' with  */
/* the same log and floor it uses for absorption data.  mu*T has the    */
/* units of -ln(I/I0), so the reconstruction still yields mu.           */
/*                                                                      */
/* I/I_in is supplied by the caller as (I - dark)/I0, i.e. the usual    */
/* dark- and flat-field-corrected transmission.                         */
/*                                                                      */
/* Boundary handling: the kernel's impulse response decays as           */
/* K0(r/rc) with rc = sqrt(R2 d/mu), so the frame is mirror-padded by   */
/* PAGANIN_PAD_RC * rc (at least PAGANIN_PAD_MIN pixels) on every side  */
/* before the FFT.  Mirror (symmetric) padding keeps the field          */
/* continuous at the frame edge whether or not the sample extends       */
/* beyond the field of view.  The padded dimensions are rounded up to   */
/* the next power of two (the built-in FFT is radix-2).                 */
/*                                                                      */
/* Parameter file (pr.par), five numbers, one per line:                 */
/*   MU      linear attenuation coefficient     [cm^-1]   > 0           */
/*   DELTA   refractive index decrement         [-]       > 0           */
/*   R1      source - sample distance           [cm]      > 0 (unused)  */
/*   R2      sample - detector distance         [cm]      >= 0          */
/*   P_SIZE  effective pixel size               [um]      > 0           */
/* A missing line, a non-number or an out-of-range value is an error.   */
/* R1 is accepted for compatibility; with parallel illumination M = 1.  */
/*                                                                      */
/* FFT backend: a self-contained radix-2 complex FFT - no external      */
/* library.  The transform runs in SINGLE precision and exploits the    */
/* real-valued input: two real rows are packed into one complex row     */
/* FFT, and only the Mx/2+1 non-redundant (Hermitian) columns are       */
/* transformed.  Single precision keeps ~1e-6 relative accuracy, far    */
/* below the photon noise of the transmission data.  The loops are      */
/* parallelised with OpenMP when compiled with /openmp (MSVC) or        */
/* -fopenmp (GCC); without it the code runs single-threaded.  Thread    */
/* count: PAGANIN_THREADS environment variable if set and > 0,          */
/* otherwise the number of online CPUs (independent of OMP_NUM_THREADS  */
/* so the ring-removal setting used elsewhere does not leak in here).   */
/* Memory: about 12 bytes per padded pixel (spectrum + kernel), e.g. a  */
/* 2048 x 2048 frame padded to 4096 x 4096 needs ~100 MB.               */
/*----------------------------------------------------------------------*/

#ifndef PAGANIN_H
#define PAGANIN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#ifdef _OPENMP
#  include <omp.h>
#  ifdef _WIN32
#    define WIN32_LEAN_AND_MEAN
#    define NOMINMAX
#    include <windows.h>
#  else
#    include <unistd.h>
#  endif
#endif

#ifndef M_PI
#define M_PI	3.14159265358979323846
#endif

#define PAGANIN_PAD_RC	10.0		/* pad = PAD_RC * rc  (exp(-10) ~ 5e-5) */
#define PAGANIN_PAD_MIN	32		/* but never less than this many pixels */
#define PAGANIN_MAX_N	32768		/* largest padded dimension accepted */

typedef struct {
	/* pr.par */
	double		MU, DELTA, R1, R2;	/* cm^-1, -, cm, cm */
	double		P_SIZE_um;		/* as given */
	double		W;			/* pixel pitch in cm */
	double		c;			/* R2 DELTA / MU, cm^2 */
	/* geometry */
	int		nthreads;		/* threads for the loops below */
	int		Nx, Ny;			/* frame */
	int		px, py;			/* pad on each side */
	int		Mx, My;			/* padded FFT size (powers of two) */
	int		nk;			/* Mx/2 + 1 stored columns */
	int		sc1;			/* scratch length, max(Mx, My) */
	double		rc_px;			/* sqrt(R2 d/mu) / W */
	/* work areas (single precision) */
	float		*sr, *si;		/* My * nk half spectrum, split re/im */
	float		*filt;			/* nk * My kernel, column-major */
	float		*twxr, *twxi;		/* exp(-2 pi i k / Mx), k < Mx/2 */
	float		*twyr, *twyi;		/* exp(-2 pi i k / My), k < My/2 */
	float		*scr;			/* nthreads * 2 * sc1 scratch */
} Paganin;

/*----------------------------------------------------------------------*/
/* pr.par: read, validate, convert units.  0 on success, -1 on error     */
/* (message on stderr).  *pg is untouched on error.                      */
static int PaganinReadPar(const char *path, Paganin *pg)
{
	static const char	*name[5] = {"MU", "DELTA", "R1", "R2", "P_SIZE"};
	FILE	*f;
	char	line[256], *end;
	double	v[5];
	int	i;

	if ((f = fopen(path, "r")) == NULL) {
		(void)fprintf(stderr, "paganin: cannot open %s\n", path);
		return -1;
	}
	for (i = 0; i < 5; ++i) {
		if (fgets(line, sizeof(line), f) == NULL) {
			(void)fprintf(stderr, "paganin: %s: line %d (%s) is missing\n",
			    path, i + 1, name[i]);
			fclose(f);
			return -1;
		}
		errno = 0;
		v[i] = strtod(line, &end);
		if (end == line || errno != 0 || !isfinite(v[i])) {
			line[strcspn(line, "\r\n")] = '\0';
			(void)fprintf(stderr, "paganin: %s: line %d (%s) is not a number: '%s'\n",
			    path, i + 1, name[i], line);
			fclose(f);
			return -1;
		}
		while (*end == ' ' || *end == '\t') ++end;
		if (*end != '\0' && *end != '\r' && *end != '\n' &&
		    *end != '#' && !(end[0] == '/' && end[1] == '/')) {
			line[strcspn(line, "\r\n")] = '\0';
			(void)fprintf(stderr, "paganin: %s: line %d (%s) has trailing text: '%s'\n",
			    path, i + 1, name[i], line);
			fclose(f);
			return -1;
		}
	}
	fclose(f);

	if (!(v[0] > 0.0)) { (void)fprintf(stderr, "paganin: MU must be > 0 (got %g)\n", v[0]); return -1; }
	if (!(v[1] > 0.0)) { (void)fprintf(stderr, "paganin: DELTA must be > 0 (got %g)\n", v[1]); return -1; }
	if (!(v[2] > 0.0)) { (void)fprintf(stderr, "paganin: R1 must be > 0 (got %g)\n", v[2]); return -1; }
	if (!(v[3] >= 0.0)) { (void)fprintf(stderr, "paganin: R2 must be >= 0 (got %g)\n", v[3]); return -1; }
	if (!(v[4] > 0.0)) { (void)fprintf(stderr, "paganin: P_SIZE [um] must be > 0 (got %g)\n", v[4]); return -1; }
	if (v[4] < 1.0e-3) {		/* below 1 nm: almost certainly an old cm-unit file */
		(void)fprintf(stderr,
		    "paganin: P_SIZE = %g um is below 1 nm; pr.par takes the pixel size "
		    "in micrometres (old files used cm: %g cm = %g um)\n",
		    v[4], v[4], v[4] * 1.0e4);
		return -1;
	}

	pg->MU = v[0]; pg->DELTA = v[1]; pg->R1 = v[2]; pg->R2 = v[3];
	pg->P_SIZE_um = v[4];
	pg->W = v[4] * 1.0e-4;			/* um -> cm */
	return 0;
}

/*----------------------------------------------------------------------*/
/* smallest power of two >= v (radix-2 FFT) */
static int pg_pow2_size(int v)
{
	int	n;

	for (n = 2; n < v; n <<= 1) ;
	return n;
}

/* mirror (symmetric) index: ..., 2, 1, 0 | 0, 1, ..., n-1 | n-1, n-2, ... */
static int pg_mirror(int i, int n)
{
	int	p = 2 * n;

	i %= p;
	if (i < 0) i += p;
	return (i < n) ? i : p - 1 - i;
}

/*----------------------------------------------------------------------*/
/* threads: PAGANIN_THREADS if set and > 0, else number of online CPUs.  */
/* Only meaningful when compiled with OpenMP; 1 otherwise.               */
static int pg_thread_count(void)
{
	int	n = 1;
#ifdef _OPENMP
	const char	*e = getenv("PAGANIN_THREADS");

	if (e != NULL && atoi(e) > 0) {
		n = atoi(e);
	} else {
#  ifdef _WIN32
		SYSTEM_INFO	si;
		GetSystemInfo(&si);
		n = (int)si.dwNumberOfProcessors;
#  else
		long	c = sysconf(_SC_NPROCESSORS_ONLN);
		if (c > 0) n = (int)c;
#  endif
	}
	if (n < 1) n = 1;
#endif
	return n;
}

/*----------------------------------------------------------------------*/
/* forward twiddle table: tw[k] = exp(-2 pi i k / n), k = 0 .. n/2-1     */
/* (computed in double, stored in float)                                 */
static int pg_twiddle(int n, float **twr, float **twi)
{
	int	k;
	double	a;

	*twr = (float *)malloc((size_t)(n / 2) * sizeof(float));
	*twi = (float *)malloc((size_t)(n / 2) * sizeof(float));
	if (*twr == NULL || *twi == NULL) return -1;
	for (k = 0; k < n / 2; ++k) {
		a = 2.0 * M_PI * (double)k / (double)n;
		(*twr)[k] = (float)cos(a);
		(*twi)[k] = (float)(-sin(a));
	}
	return 0;
}

/*----------------------------------------------------------------------*/
/* in-place radix-2 DIT complex FFT, single precision.  n is a power of  */
/* two >= 2.  sign -1: forward exp(-2 pi i kn/N); sign +1: inverse       */
/* WITHOUT the 1/n factor (the caller normalises once).                  */
static void pg_fftf(int n, int sign, const float *twr, const float *twi,
		    float *xr, float *xi)
{
	int	i, j, k, t, len, half, step, p;
	float	ur, ui, vr, vi, wr, wi, s;

	/* bit-reversal permutation */
	for (i = 1, j = 0; i < n; ++i) {
		for (k = n >> 1; j >= k; k >>= 1) j -= k;
		j += k;
		if (i < j) {
			s = xr[i]; xr[i] = xr[j]; xr[j] = s;
			s = xi[i]; xi[i] = xi[j]; xi[j] = s;
		}
	}
	/* butterflies */
	for (len = 2; len <= n; len <<= 1) {
		half = len >> 1;
		step = n / len;
		for (i = 0; i < n; i += len) {
			for (j = 0, p = 0; j < half; ++j, p += step) {
				wr = twr[p];
				wi = (sign < 0) ? twi[p] : -twi[p];
				t  = i + j;
				vr = xr[t + half] * wr - xi[t + half] * wi;
				vi = xr[t + half] * wi + xi[t + half] * wr;
				ur = xr[t];
				ui = xi[t];
				xr[t]        = ur + vr;
				xi[t]        = ui + vi;
				xr[t + half] = ur - vr;
				xi[t + half] = ui - vi;
			}
		}
	}
}

/*----------------------------------------------------------------------*/
/* Set up padding, twiddles, kernel and buffers for an Nx x Ny frame.    */
/* 0 or -1.                                                              */
static int PaganinInit(Paganin *pg, int Nx, int Ny)
{
	int	j, k, mx, my, pad;
	double	rc, dkx, dky, k1, k2x;
	size_t	nn;

	pg->Nx = Nx;
	pg->Ny = Ny;
	pg->c  = pg->R2 * pg->DELTA / pg->MU;		/* cm^2 */

	/* real-space decay length of the kernel, in pixels */
	rc = sqrt(pg->c);				/* cm; 0 when R2 = 0 */
	pg->rc_px = rc / pg->W;
	pad = (int)ceil(PAGANIN_PAD_RC * pg->rc_px);
	if (pad < PAGANIN_PAD_MIN) pad = PAGANIN_PAD_MIN;
	pg->px = pg->py = pad;
	pg->Mx = pg_pow2_size(Nx + 2 * pad);
	pg->My = pg_pow2_size(Ny + 2 * pad);
	if (pg->Mx > PAGANIN_MAX_N || pg->My > PAGANIN_MAX_N) {
		(void)fprintf(stderr,
		    "paganin: padded field %d x %d exceeds %d "
		    "(frame %d x %d, pad %d = %.0f x rc, rc = %.3g px)\n",
		    pg->Mx, pg->My, PAGANIN_MAX_N, Nx, Ny, pad, PAGANIN_PAD_RC,
		    pg->rc_px);
		(void)fprintf(stderr,
		    "paganin: a pad this large usually means P_SIZE has the wrong "
		    "units: pr.par takes it in MICROMETRES here (%g um given); the "
		    "old mkpms/mkpmbg chain used cm - a cm value makes rc 10^4 "
		    "pixels too large\n", pg->P_SIZE_um);
		return -1;
	}
	pg->nk  = pg->Mx / 2 + 1;
	pg->sc1 = (pg->Mx > pg->My) ? pg->Mx : pg->My;
	nn = (size_t)pg->My * pg->nk;

	pg->nthreads = pg_thread_count();

	pg->sr   = (float *)malloc(nn * sizeof(float));
	pg->si   = (float *)malloc(nn * sizeof(float));
	pg->filt = (float *)malloc(nn * sizeof(float));
	pg->scr  = (float *)malloc((size_t)pg->nthreads * 2 * pg->sc1 * sizeof(float));
	if (pg->sr == NULL || pg->si == NULL || pg->filt == NULL || pg->scr == NULL ||
	    pg_twiddle(pg->Mx, &pg->twxr, &pg->twxi) != 0 ||
	    pg_twiddle(pg->My, &pg->twyr, &pg->twyi) != 0) {
		(void)fprintf(stderr, "paganin: no memory for FFT buffers (%d x %d)\n",
		    pg->Mx, pg->My);
		return -1;
	}

	/* kernel 1 / (1 + (R2 d/mu) |k|^2), column-major for the column pass;
	   frequency index m = i (i <= M/2), i - M (i > M/2) */
	dkx = 2.0 * M_PI / (pg->W * (double)pg->Mx);	/* cm^-1 per index */
	dky = 2.0 * M_PI / (pg->W * (double)pg->My);
	for (k = 0; k < pg->nk; ++k) {
		mx  = k;				/* 0 .. Mx/2 */
		k1  = dkx * mx;
		k2x = k1 * k1;
		for (j = 0; j < pg->My; ++j) {
			my = (j <= pg->My / 2) ? j : j - pg->My;
			k1 = dky * my;
			pg->filt[(size_t)k * pg->My + j] =
			    (float)(1.0 / (1.0 + pg->c * (k2x + k1 * k1)));
		}
	}

	(void)fprintf(stderr,
	    "paganin: MU=%g /cm DELTA=%g R2=%g cm P_SIZE=%g um; "
	    "rc=%.1f px, pad=%d, FFT %d x %d (built-in radix-2, real-packed, "
	    "single precision, %.0f MB), %d thread%s\n",
	    pg->MU, pg->DELTA, pg->R2, pg->P_SIZE_um,
	    pg->rc_px, pad, pg->Mx, pg->My,
	    (double)nn * 3.0 * sizeof(float) / (1024.0 * 1024.0),
	    pg->nthreads, pg->nthreads == 1 ? "" : "s");
	return 0;
}

/*----------------------------------------------------------------------*/
/* In place: T[Nx*Ny] = (I-dark)/I0  ->  filtered transmission T'.       */
/*                                                                      */
/* Row pass: two mirror-padded real rows are packed into one complex    */
/* FFT of length Mx and untangled into the half spectra (k = 0..Mx/2)   */
/* of both rows.  Column pass: each of the nk stored columns gets a     */
/* forward FFT, the kernel, and the inverse FFT in one sweep.  Inverse  */
/* row pass: only the row pairs that intersect the frame are rebuilt    */
/* (Hermitian extension) and transformed back.                          */
static void PaganinFilter(Paganin *pg, double *T)
{
	int	Nx = pg->Nx, Ny = pg->Ny, Mx = pg->Mx, My = pg->My;
	int	px = pg->px, py = pg->py, nk = pg->nk;
	int	jp0 = py / 2, jp1 = (py + Ny - 1) / 2;
	int	jp, k;
	float	inv = 1.0f / ((float)Mx * (float)My);

	/* rows, forward: pack two real rows as one complex FFT */
#ifdef _OPENMP
#pragma omp parallel num_threads(pg->nthreads)
#endif
	{
		int	tid = 0, i, kk;
		float	*cr, *ci;
#ifdef _OPENMP
		tid = omp_get_thread_num();
#endif
		cr = pg->scr + (size_t)tid * 2 * pg->sc1;
		ci = cr + pg->sc1;
#ifdef _OPENMP
#pragma omp for
#endif
		for (jp = 0; jp < My / 2; ++jp) {
			const double	*r0 = T + (size_t)pg_mirror(2*jp   - py, Ny) * Nx;
			const double	*r1 = T + (size_t)pg_mirror(2*jp+1 - py, Ny) * Nx;
			float	*s0r = pg->sr + (size_t)(2*jp) * nk;
			float	*s0i = pg->si + (size_t)(2*jp) * nk;
			float	*s1r = s0r + nk, *s1i = s0i + nk;

			for (i = 0; i < Mx; ++i) {
				int	m = pg_mirror(i - px, Nx);
				cr[i] = (float)r0[m];
				ci[i] = (float)r1[m];
			}
			pg_fftf(Mx, -1, pg->twxr, pg->twxi, cr, ci);
			for (kk = 0; kk < nk; ++kk) {
				int	mk = (Mx - kk) & (Mx - 1);
				s0r[kk] = 0.5f * (cr[kk] + cr[mk]);
				s0i[kk] = 0.5f * (ci[kk] - ci[mk]);
				s1r[kk] = 0.5f * (ci[kk] + ci[mk]);
				s1i[kk] = 0.5f * (cr[mk] - cr[kk]);
			}
		}
	}

	/* columns: forward FFT, kernel, inverse FFT in one sweep */
#ifdef _OPENMP
#pragma omp parallel num_threads(pg->nthreads)
#endif
	{
		int	tid = 0, j;
		float	*cr, *ci;
#ifdef _OPENMP
		tid = omp_get_thread_num();
#endif
		cr = pg->scr + (size_t)tid * 2 * pg->sc1;
		ci = cr + pg->sc1;
#ifdef _OPENMP
#pragma omp for
#endif
		for (k = 0; k < nk; ++k) {
			const float	*ft = pg->filt + (size_t)k * My;

			for (j = 0; j < My; ++j) {
				cr[j] = pg->sr[(size_t)j * nk + k];
				ci[j] = pg->si[(size_t)j * nk + k];
			}
			pg_fftf(My, -1, pg->twyr, pg->twyi, cr, ci);
			for (j = 0; j < My; ++j) {
				cr[j] *= ft[j];
				ci[j] *= ft[j];
			}
			pg_fftf(My, +1, pg->twyr, pg->twyi, cr, ci);
			for (j = 0; j < My; ++j) {
				pg->sr[(size_t)j * nk + k] = cr[j];
				pg->si[(size_t)j * nk + k] = ci[j];
			}
		}
	}

	/* rows, inverse: only the pairs that intersect the frame */
#ifdef _OPENMP
#pragma omp parallel num_threads(pg->nthreads)
#endif
	{
		int	tid = 0, i, kk;
		float	*cr, *ci;
#ifdef _OPENMP
		tid = omp_get_thread_num();
#endif
		cr = pg->scr + (size_t)tid * 2 * pg->sc1;
		ci = cr + pg->sc1;
#ifdef _OPENMP
#pragma omp for
#endif
		for (jp = jp0; jp <= jp1; ++jp) {
			const float	*s0r = pg->sr + (size_t)(2*jp) * nk;
			const float	*s0i = pg->si + (size_t)(2*jp) * nk;
			const float	*s1r = s0r + nk, *s1i = s0i + nk;
			int	j0 = 2 * jp, j1 = 2 * jp + 1;

			/* C = S0 + i S1; k > Mx/2 by Hermitian symmetry */
			for (kk = 0; kk < nk; ++kk) {
				cr[kk] = s0r[kk] - s1i[kk];
				ci[kk] = s0i[kk] + s1r[kk];
			}
			for (kk = nk; kk < Mx; ++kk) {
				int	mk = Mx - kk;
				cr[kk] = s0r[mk] + s1i[mk];
				ci[kk] = s1r[mk] - s0i[mk];
			}
			pg_fftf(Mx, +1, pg->twxr, pg->twxi, cr, ci);

			if (j0 >= py && j0 < py + Ny)
				for (i = 0; i < Nx; ++i)
					T[(size_t)(j0 - py) * Nx + i] = cr[i + px] * inv;
			if (j1 >= py && j1 < py + Ny)
				for (i = 0; i < Nx; ++i)
					T[(size_t)(j1 - py) * Nx + i] = ci[i + px] * inv;
		}
	}
}

/*----------------------------------------------------------------------*/
static void PaganinFree(Paganin *pg)
{
	free(pg->sr);   free(pg->si);   free(pg->filt);
	free(pg->twxr); free(pg->twxi);
	free(pg->twyr); free(pg->twyi);
	free(pg->scr);
	memset(pg, 0, sizeof(*pg));
}

#endif	/* PAGANIN_H */
