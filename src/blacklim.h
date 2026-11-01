/*----------------------------------------------------------------------*/
/* blacklim.h - low-transmission guard for the sinogram builders.        */
/*                                                                      */
/* CT_REC_BLACK_THRESH is a PER-PIXEL floor on the measured signal, in   */
/* counts above dark.  A pixel at or below it is treated as opaque and   */
/* the floor is used in its place, so                                    */
/*                                                                      */
/*     p = log(I0 / (I - dark))                                          */
/*                                                                      */
/* can no longer become +Inf (I-dark == 0) or NaN (I-dark < 0).  Those   */
/* values used to reach the ring removal, where the GPU per-column sort  */
/* pads the column with FLT_MAX / index -1 and relies on the padding     */
/* sorting to the top.  +Inf beats FLT_MAX and NaN compares false        */
/* against everything, so the padding sank into the real ranks and       */
/* perm[] came out holding (unsigned)-1; the scatter in median_scatter() */
/* then wrote about 35 TB past the buffer and CUDA reported              */
/*                                                                      */
/*     CUDA error sort_filter_g.cu:266: an illegal memory access ...     */
/*     ring removal image processing failed                              */
/*                                                                      */
/* on the first synchronising call after the kernel.  The CPU/OpenMP     */
/* backend did not crash but sorted such columns incorrectly, so the     */
/* ring removal was silently wrong there instead.                        */
/*                                                                      */
/* CT_REC_BLACK_FRAC is the fraction of clipped pixels above which a     */
/* whole projection is judged black (the "missing angle" case) and       */
/* set to 0 (in every builder since 2026-10-27; ct_rec and the hp_tg     */
/* family used to substitute the mean).  It replaces the older           */
/* test on the line average, which could not work here: with open-beam   */
/* margins in the field of view - which the air reference at both ends   */
/* of the line requires anyway - the average has a hard floor of         */
/* (1 - coverage) * signal and never approaches the threshold, however   */
/* opaque the centre of the sample becomes.  For a uniformly attenuated  */
/* line (a plate seen edge-on) the two criteria fire at the same         */
/* transmission, so plate-like samples keep behaving as before.          */
/*                                                                      */
/* CT_REC_BLACK_MEAN (2026-10-27) brings the old average test back as an  */
/* additional criterion: a projection whose mean (I - dark) over the line */
/* (ct_rec, ofct_rec, sinog, ...) or over the whole frame (the HiPic      */
/* readers used by hp_tg, ofct_srec, ...) is below it is black as well.   */
/* It is off unless set; it does not touch the per-pixel floor, and a     */
/* black projection is set to 0 whichever test fired.  The value is in    */
/* counts above dark, exactly like CT_REC_BLACK_THRESH before 2026-08-06, */
/* so an old "CT_REC_BLACK_THRESH=2000" becomes "CT_REC_BLACK_MEAN=2000". */
/*                                                                      */
/* All variables are read once per process and validated.  A normal run  */
/* prints nothing: the values in effect are echoed only when BOTH        */
/* thresholds have been exceeded (a projection was judged black), and    */
/* the end-of-run summary appears only when at least one pixel actually  */
/* hit the floor.                                                        */
/*----------------------------------------------------------------------*/

#ifndef BLACKLIM_H
#define BLACKLIM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* counts above dark; 2 caps the absorbance at log(I0/2) */
#define BLACK_THRESH_DEFAULT	2.0
/* clipped fraction that makes a projection "black" */
#define BLACK_FRAC_DEFAULT	0.5
/* line/frame average below which a projection is black; 0 = test off */
#define BLACK_MEAN_DEFAULT	0.0

static double	blackThresh	= -1.0;	/* <0 = not read yet */
static double	blackFrac	= -1.0;
static double	blackMean	= -1.0;
static double	blackMinSignal	= 1.0e300;	/* smallest (I-dark) seen */
static double	blackMaxP	= 0.0;	/* largest absorbance emitted */
static long	blackClipped	= 0;	/* pixels that hit the floor */
static long	blackSeen	= 0;	/* pixels converted */
static long	blackProj	= 0;	/* projections judged black */

/*----------------------------------------------------------------------*/
/* Per-pixel signal floor, counts above dark.  A value <= 0 would let    */
/* Inf/NaN back in, so it is rejected rather than honoured.              */
static double BlackThresh(void)
{
	char	*e;
	double	v;

	if (blackThresh >= 0.0) return blackThresh;

	blackThresh = BLACK_THRESH_DEFAULT;
	if ((e = getenv("CT_REC_BLACK_THRESH")) != NULL) {
		v = atof(e);
		if (v > 0.0) {
			blackThresh = v;
		} else {
			(void)fprintf(stderr,
			    "CT_REC_BLACK_THRESH: '%s' is not a positive number, "
			    "using %.4g\n", e, blackThresh);
		}
	}
	return blackThresh;
}

/*----------------------------------------------------------------------*/
/* Fraction of clipped pixels above which the projection is black.       */
static double BlackFrac(void)
{
	char	*e;
	double	v;

	if (blackFrac >= 0.0) return blackFrac;

	blackFrac = BLACK_FRAC_DEFAULT;
	if ((e = getenv("CT_REC_BLACK_FRAC")) != NULL) {
		v = atof(e);
		if (v > 0.0 && v <= 1.0) {
			blackFrac = v;
		} else {
			(void)fprintf(stderr,
			    "CT_REC_BLACK_FRAC: '%s' is not in (0,1], using %.4g\n",
			    e, blackFrac);
		}
	}
	return blackFrac;
}

/*----------------------------------------------------------------------*/
/* Average (I - dark) below which the projection is black; 0 = off.      */
/* Unlike the two guards above this one may be switched off, because it  */
/* only discards projections and never lets Inf/NaN through.             */
static double BlackMean(void)
{
	char	*e, *end;
	double	v;

	if (blackMean >= 0.0) return blackMean;

	blackMean = BLACK_MEAN_DEFAULT;
	if ((e = getenv("CT_REC_BLACK_MEAN")) != NULL && *e != '\0') {
		v = strtod(e, &end);
		if (end == e) {
			(void)fprintf(stderr,
			    "CT_REC_BLACK_MEAN: '%s' is not a number, "
			    "the average test stays off\n", e);
		} else if (v > 0.0) {
			blackMean = v;
		}
	}
	return blackMean;
}

/* True when the average test is on and sum/n (the mean (I - dark) of the */
/* line or frame, n pixels) is below it.                                 */
static int BlackMeanHit(double sum, double n)
{
	double	m = BlackMean();

	return m > 0.0 && sum/n < m;
}

/* Suffix for the "black" warning: the average and its limit when the    */
/* test is on, an empty string otherwise.                                */
static const char *BlackMeanNote(double sum, double n, char *buf, size_t len)
{
	double	m = BlackMean();

	if (m > 0.0)
		(void)snprintf(buf, len, " avg=%.1f (mean=%.4g)", sum/n, m);
	else
		buf[0] = '\0';
	return buf;
}

/*----------------------------------------------------------------------*/
/* Absorbance of one pixel with the floor applied.                       */
/*   num  = I0                     (counts above dark)                   */
/*   den  = I - dark               (counts above dark, may be <= 0)      */
/*   nclip, if not NULL, is incremented when the floor had to be used.   */
/* The return value is always finite.                                    */
static double BlackLog(double num, double den, double flr, int *nclip)
{
	double	p;

	++blackSeen;
	if (den < blackMinSignal) blackMinSignal = den;

	/* !(x >= flr) instead of (x < flr): every comparison with NaN is
	   false, so a plain (x < flr) would let NaN straight through.  The
	   I0 interpolation divides by (t2 - t1) and produces NaN/Inf when
	   two I0 timestamps coincide, so the inputs are not guaranteed to
	   be ordinary numbers.  The upper clamp bounds +Inf the same way. */
	if (!(den >= flr)) {
		den = flr;
		++blackClipped;
		if (nclip != NULL) ++*nclip;
	} else if (den > 1.0e300) {
		den = 1.0e300;
	}
	if (!(num >= flr)) {
		num = flr;
	} else if (num > 1.0e300) {
		num = 1.0e300;
	}

	p = log(num/den);
	if (p > blackMaxP) blackMaxP = p;
	return p;
}

/*----------------------------------------------------------------------*/
/* Counts one projection that was judged black.  The first one triggers  */
/* the threshold echo: both limits were exceeded, so the values that did */
/* it become worth showing.  Quiet runs stay quiet.                      */
static void BlackCountProjection(void)
{
	if (blackProj == 0L) {
		(void)fprintf(stderr,
		    "CT_REC_BLACK_THRESH = %.4g counts above dark (per-pixel floor, "
		    "max absorbance = log(I0/%.4g))\n", BlackThresh(), BlackThresh());
		(void)fprintf(stderr,
		    "CT_REC_BLACK_FRAC   = %.4g of the line\n", BlackFrac());
		if (BlackMean() > 0.0)
			(void)fprintf(stderr,
			    "CT_REC_BLACK_MEAN   = %.4g counts above dark "
			    "(line/frame average)\n", BlackMean());
	}
	++blackProj;
}

/*----------------------------------------------------------------------*/
/* One-line summary; call once the sinogram is complete.  Turns what used */
/* to be a CUDA crash into a readable data-quality statement.  Prints    */
/* nothing when no pixel hit the floor.                                  */
static void BlackReport(void)
{
	if (blackSeen == 0L || blackClipped == 0L) return;

	(void)fprintf(stderr,
	    "low-transmission guard: clipped %ld / %ld pixel(s) (%.3f%%), "
	    "minimum signal %.1f counts, max absorbance %.3f, "
	    "%ld projection(s) judged black\n",
	    blackClipped, blackSeen,
	    100.0*(double)blackClipped/(double)blackSeen,
	    blackMinSignal, blackMaxP, blackProj);

	if (blackClipped > 0L) {
		(void)fprintf(stderr,
		    "  the sample is too thick or too dense for the current "
		    "exposure at those pixels; raise the exposure or lower the "
		    "energy if the clipped fraction matters\n");
	}
}

#endif	/* BLACKLIM_H */
