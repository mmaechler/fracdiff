/*-*- mode: C; kept-old-versions: 12;  kept-new-versions: 20; -*-
 *
 * fdsim.f -- translated by f2c (version 20031025).
 *
 * and produced by  f2c-clean,v 1.10 2002/03/28 16:37:27 maechler
 */
#include <Rmath.h>

#include "fracdiff.h"

extern double dgamr_(double *);
extern double dgamma_(double *);

/* Common Block Declarations --- included as "extern" */
#define FD_EXTERNAL extern

#include "mach_comm.h"
#include "gamm_comm.h"

/* Subroutine */
void fdsim(int *n, int *ip, int *iq, double *ar, double *ma,
	   double *d__, double *mu, double *y, double *s,
	   double *flmin, double *flmax, double *epmin, double *epmax)
{
/*  generates a random time series for use with fracdf

  Input :

  n	 int  length of the time series
  ip	 int  number of autoregressive parameters
  iq	 int  number of moving average parameters
  ar	 float	  (ip) autoregressive parameters
  ma	 float	  (iq) moving average parameters
  d	 float	   fractional differencing parameter
  rmu	 float	   time series mean
  y	 float	  (n+iq) 1st n : normalized random numbers
  s	 float	  (n+iq) workspace

  Output :

  s	 float	 (n) the generated time series
 -----------------------------------------------------------------------------

	Simulates a series of length n from an ARIMA (p,d,q) model
	with fractional d (0 < d < 0.5).

 -----------------------------------------------------------------------------
     float		 ar(ip), ma(iq), d, rmu
     float		 y(n+iq), s(n+iq)
 --------------------------------------------------------------------------
 */

    /* System generated locals */
    double d__1;

    /* Local variables */
    int i, j, k;
    double dj, vk, dk1, amk, sum, dk1d, temp;

    /*	   Parameter adjustments */
    --y;
    --s;

    /* Common Block -- Initializations: Input & Output for gamma() functions */
    gammfd_.igamma = 0;
    gammfd_.jgamma = 0;
    machfd_.fltmin = *flmin;
    machfd_.fltmax = *flmax;
    machfd_.epsmin = *epmin;
    machfd_.epsmax = *epmax;

    /* Calculate vk[0] = 'g0' */
    d__1 = 1. - *d__;
    temp = dgamr_(&d__1);
    if (gammfd_.igamma != 0) {
	for (i = 1; i <= *n; ++i)
	    s[i] = 0.;
	return;
    }
    /* else : */
    d__1 = 1. - *d__ * 2.;
    vk = dgamma_(&d__1) * (temp * temp);
    if (gammfd_.igamma != 0) {
	for (i = 1; i <= *n; ++i)
	    s[i] = 0.;
	return;
    }
    /* else -- Gamma values ok, compute	 : */

    /*	 Generate y(1) */

    y[1] *= sqrt(vk);

/*	 Generate y(2) and initialise vk,phi(j) */

    temp = *d__ / (1. - *d__);
    vk *= 1. - temp * temp;
    amk = temp * y[1];
    s[1] = temp;
    y[2] = amk + y[2] * sqrt(vk);

/*	 Generate y(3),...,y(n+iq) */

    for (k = 3; k <= (*n + *iq); ++k) {
	dk1 = (double) k - 1.;
	dk1d = dk1 - *d__;

	/*	 Update the phi(j) using the recursion formula on W498 */

	for (j = 1; j <= (k - 2); ++j) {
	    dj = dk1 - (double) j;
	    s[j] *= dk1 * (dj - *d__) / (dk1d * dj);
	}
	temp = *d__ / dk1d;
	s[k - 1] = temp;

	/*	 Update vk */

	vk *= 1. - temp * temp;

	/*	 Form amk */
	amk = 0.;
	for (j = 1; j <= (k - 1); ++j)
	    amk += s[j] * y[k - j];

	/*	 Generate y(k) */

	y[k] = amk + y[k] * sqrt(vk);
    }

/* We now have an ARIMA (0,d,0) realisation of length n+iq in
  y(k), k=1,n+iq. We now run this through an inverse ARMA(p,q)
	 filter to get the final output in s(k), k=1,n. */

    for (k = 1; k <= *n; ++k) {
	sum = 0.;
	j = imin2(*ip,k);
	for (i = 0; i < j; ++i)
	    sum += ar[i] * s[k - i - 1];
	for (j = 0; j < *iq; ++j)
	    sum -= ma[j] * y[k + *iq - j - 1];
	s[k] = sum + y[k + *iq];
    }
    /* now add the global mean */
    if (*mu != 0.) {
	for (i = 1; i <= *n; ++i)
	    s[i] += *mu;
    }
    return;
} /* fdsim */

