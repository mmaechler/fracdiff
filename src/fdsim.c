/* fdsim.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal fltmin, fltmax, epsmin, epsmax;
} machfd_;

#define machfd_1 machfd_

struct {
    integer igamma, jgamma;
} gammfd_;

#define gammfd_1 gammfd_

/* Subroutine */ int fdsim_(integer *n, integer *ip, integer *iq, doublereal *
	ar, doublereal *ma, doublereal *d__, doublereal *rmu, doublereal *y, 
	doublereal *s, doublereal *flmin, doublereal *flmax, doublereal *
	epmin, doublereal *epmax)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal dj, vk, dk1, amk, sum, dk1d, temp;
    extern doublereal dgamr_(doublereal *), dgamma_(doublereal *);

/*  generates a random time series for use with fracdf */

/*  Input : */

/*  n      integer  length of the time series */
/*  ip     integer  number of autoregressive parameters */
/*  iq     integer  number of moving average parameters */
/*  ar     real    (ip) autoregressive parameters */
/*  ma     real    (iq) moving average parameters */
/*  d      real     fractional differencing parameter */
/*  rmu    real     time series mean */
/*  y      real    (n+iq) 1st n : normalized random numbers */
/*  s      real    (n+iq) workspace */

/*  Output : */

/*  s      real   (n) the generated time series */
/* ----------------------------------------------------------------------------- */

/*        Simulates a series of length n from an ARIMA (p,d,q) model */
/*        with fractional d (0 < d < 0.5). */

/* ----------------------------------------------------------------------------- */
/*     real               ar(ip), ma(iq), d, rmu */
/*     real               y(n+iq), s(n+iq) */
/* -------------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ar;
    --ma;
    --y;
    --s;

    /* Function Body */
    gammfd_1.igamma = 0;
    gammfd_1.jgamma = 0;
    machfd_1.fltmin = *flmin;
    machfd_1.fltmax = *flmax;
    machfd_1.epsmin = *epmin;
    machfd_1.epsmax = *epmax;

/* 	 Calculate vk[0] = 'g0' */
    d__1 = 1. - *d__;
    temp = dgamr_(&d__1);
    if (gammfd_1.igamma != 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[i__] = 0.;
	}
	return 0;
    }
    d__1 = 1. - *d__ * 2.;
    vk = dgamma_(&d__1) * (temp * temp);
    if (gammfd_1.igamma != 0) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[i__] = 0.;
	}
	return 0;
    }

/* 	 Generate y(1) */

    y[1] *= sqrt(vk);

/* 	 Generate y(2) and initialise vk,phi(j) */

    temp = *d__ / (1. - *d__);
    vk *= 1. - temp * temp;
    amk = temp * y[1];
    s[1] = temp;
    y[2] = amk + y[2] * sqrt(vk);

/* 	 Generate y(3),...,y(n+iq) */

    i__1 = *n + *iq;
    for (k = 3; k <= i__1; ++k) {
	dk1 = (doublereal) k - 1.;
	dk1d = dk1 - *d__;

/* 	 Update the phi(j) using the recursion formula on W498 */

	i__2 = k - 2;
	for (j = 1; j <= i__2; ++j) {
	    dj = dk1 - (doublereal) j;
	    s[j] *= dk1 * (dj - *d__) / (dk1d * dj);
	}
	temp = *d__ / dk1d;
	s[k - 1] = temp;

/* 	 Update vk */

	vk *= 1. - temp * temp;

/* 	 Form amk */

	amk = 0.;
	i__2 = k - 1;
	for (j = 1; j <= i__2; ++j) {
	    amk += s[j] * y[k - j];
	}

/* 	 Generate y(k) */

	y[k] = amk + y[k] * sqrt(vk);
    }

/* 	 We now have an ARIMA (0,d,0) realisation of length n+iq in */
/* 	 y(k),k=1,n+iq. We now run this through an inverse ARMA(p,q) */
/* 	 filter to get the final output in s(k), k=1,n. */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	sum = 0.;
	j = min(*ip,k);
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    sum += ar[i__] * s[k - i__];
	}
	i__2 = *iq;
	for (j = 1; j <= i__2; ++j) {
	    sum -= ma[j] * y[k + *iq - j];
	}
	s[k] = sum + y[k + *iq];
    }
/* now add the global mean */
    if (*rmu != 0.) {
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    s[i__] += *rmu;
	}
    }
    return 0;
} /* fdsim_ */

