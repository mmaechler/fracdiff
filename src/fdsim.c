/* Produced by
 * $Id: f2c-clean,v 1.10 2002/03/28 16:37:27 maechler Exp $
 */
/* fdsim.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip


-- be brave, try without >>>
 * #include "f2c.h" /* <<<<< ------*/
#include <math.h>
#ifndef max
# define	max(a, b) 		((a) < (b) ? (b) : (a))
#endif
#ifndef min
# define	min(a, b)		((a) > (b) ? (b) : (a))
#endif
#ifndef abs
# define	abs(x)			((x) >= 0 ? (x) : -(x))
#endif


/* Common Block Declarations */

struct {
    double fltmin, fltmax, epsmin, epsmax;
} machfd_;

#define machfd_1 machfd_

struct {
    int igamma, jgamma;
} gammfd_;

#define gammfd_1 gammfd_

/* Subroutine */ int fdsim_(int *n, int *ip, int *iq, double *
	ar, double *ma, double *d__, double *rmu, double *y,

	double *s, double *flmin, double *flmax, double *
	epmin, double *epmax)
{
    /* System generated locals */
    int i__1, i__2;
    double d__1;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int i__, j, k;
    static double dj, vk, dk1, amk, sum, dk1d, temp;
    extern double dgamr_(double *), dgamma_(double *);

/*  generates a random time series for use with fracdf 

  Input : 

  n      int  length of the time series 
  ip     int  number of autoregressive parameters 
  iq     int  number of moving average parameters 
  ar     float    (ip) autoregressive parameters 
  ma     float    (iq) moving average parameters 
  d      float     fractional differencing parameter 
  rmu    float     time series mean 
  y      float    (n+iq) 1st n : normalized random numbers 
  s      float    (n+iq) workspace 

  Output : 

  s      float   (n) the generated time series 
 ----------------------------------------------------------------------------- 

        Simulates a series of length n from an ARIMA (p,d,q) model 
        with fractional d (0 < d < 0.5). 

 ----------------------------------------------------------------------------- 
     float               ar(ip), ma(iq), d, rmu 
     float               y(n+iq), s(n+iq) 
 -------------------------------------------------------------------------- 
     Parameter adjustments */
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
	for (i__ = 1; i__ <= *n; ++i__) { /* f2c-clean: s {i__1} {*n} */
	    s[i__] = 0.;
	}
	return 0;
    }
    d__1 = 1. - *d__ * 2.;
    vk = dgamma_(&d__1) * (temp * temp);
    if (gammfd_1.igamma != 0) {
	for (i__ = 1; i__ <= *n; ++i__) { /* f2c-clean: s {i__1} {*n} */
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

    for (k = 3; k <= (*n + *iq); ++k) { /* f2c-clean: s {i__1} {*n + *iq} */
	dk1 = (double) k - 1.;
	dk1d = dk1 - *d__;

/* 	 Update the phi(j) using the recursion formula on W498 */

	for (j = 1; j <= (k - 2); ++j) { /* f2c-clean: s {i__2} {k - 2} */
	    dj = dk1 - (double) j;
	    s[j] *= dk1 * (dj - *d__) / (dk1d * dj);
	}
	temp = *d__ / dk1d;
	s[k - 1] = temp;

/* 	 Update vk */

	vk *= 1. - temp * temp;

/* 	 Form amk */

	amk = 0.;
	for (j = 1; j <= (k - 1); ++j) { /* f2c-clean: s {i__2} {k - 1} */
	    amk += s[j] * y[k - j];
	}

/* 	 Generate y(k) */

	y[k] = amk + y[k] * sqrt(vk);
    }

/* 	 We now have an ARIMA (0,d,0) realisation of length n+iq in 
 	 y(k),k=1,n+iq. We now run this through an inverse ARMA(p,q) 
 	 filter to get the final output in s(k), k=1,n. */

    for (k = 1; k <= *n; ++k) { /* f2c-clean: s {i__1} {*n} */
	sum = 0.;
	j = min(*ip,k);
	for (i__ = 1; i__ <= j; ++i__) { /* f2c-clean: s {i__2} {j} */
	    sum += ar[i__] * s[k - i__];
	}
	for (j = 1; j <= *iq; ++j) { /* f2c-clean: s {i__2} {*iq} */
	    sum -= ma[j] * y[k + *iq - j];
	}
	s[k] = sum + y[k + *iq];
    }
/* now add the global mean */
    if (*rmu != 0.) {
	for (i__ = 1; i__ <= *n; ++i__) { /* f2c-clean: s {i__1} {*n} */
	    s[i__] += *rmu;
	}
    }
    return 0;
} /* fdsim_ */

