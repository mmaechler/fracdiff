/* fdgam.f -- translated by f2c (version 20031025).

 * and produced by
 * $Id: f2c-clean,v 1.10 2002/03/28 16:37:27 maechler Exp $
 *
 * and manually pretty edited by Martin Maechler, 2004-10-01
*/

#include <Rmath.h>


#ifndef max
# define	max(a, b) 		((a) < (b) ? (b) : (a))
#endif
#ifndef min
# define	min(a, b)		((a) > (b) ? (b) : (a))
#endif
#ifndef abs
# define	abs(x)			((x) >= 0 ? (x) : -(x))
#endif

/* EXPORTS */
double dgamma_(double *x);
double dgamr_ (double *x);

static int    dlgams_(double *, double *, double *);
static double dlngam_(double *);

static void d9gaml_(double *xmin, double *xmax);
static double d9lgmc_(double *);

static double dcsevl_(double *x, double *a, int *n);

static int initds_(double *, int *, float *);


/* Common Block Declarations */

struct {
    double fltmin, fltmax, epsmin, epsmax;
} machfd_;

#define machfd_1 machfd_

struct {
    int igamma, jgamma;
} gammfd_;

#define gammfd_1 gammfd_

/* Table of constant values */

static int c__42 = 42;
static double c_b9 = 2.;
static int c__15 = 15;

double dgamma_(double *x)
{
    /* Initialized data */

    static double gamcs[42] = { .008571195590989331421920062399942,
	    .004415381324841006757191315771652,
	    .05685043681599363378632664588789,
	    -.004219835396418560501012500186624,
	    .001326808181212460220584006796352,
	    -1.893024529798880432523947023886e-4,
	    3.606925327441245256578082217225e-5,
	    -6.056761904460864218485548290365e-6,
	    1.055829546302283344731823509093e-6,
	    -1.811967365542384048291855891166e-7,
	    3.117724964715322277790254593169e-8,
	    -5.354219639019687140874081024347e-9,
	    9.19327551985958894688778682594e-10,
	    -1.577941280288339761767423273953e-10,
	    2.707980622934954543266540433089e-11,
	    -4.646818653825730144081661058933e-12,
	    7.973350192007419656460767175359e-13,
	    -1.368078209830916025799499172309e-13,
	    2.347319486563800657233471771688e-14,
	    -4.027432614949066932766570534699e-15,
	    6.910051747372100912138336975257e-16,
	    -1.185584500221992907052387126192e-16,
	    2.034148542496373955201026051932e-17,
	    -3.490054341717405849274012949108e-18,
	    5.987993856485305567135051066026e-19,
	    -1.027378057872228074490069778431e-19,
	    1.762702816060529824942759660748e-20,
	    -3.024320653735306260958772112042e-21,
	    5.188914660218397839717833550506e-22,
	    -8.902770842456576692449251601066e-23,
	    1.527474068493342602274596891306e-23,
	    -2.620731256187362900257328332799e-24,
	    4.496464047830538670331046570666e-25,
	    -7.714712731336877911703901525333e-26,
	    1.323635453126044036486572714666e-26,
	    -2.270999412942928816702313813333e-27,
	    3.896418998003991449320816639999e-28,
	    -6.685198115125953327792127999999e-29,
	    1.146998663140024384347613866666e-29,
	    -1.967938586345134677295103999999e-30,
	    3.376448816585338090334890666666e-31,
	    -5.793070335782135784625493333333e-32 };
    static double pi = 3.1415926535897932384626433832795;
    static double sq2pil = .91893853320467274178032973640562;
    static int ngam = 0;
    static double xmin = 0.;
    static double xmax = 0.;
    static double xsml = 0.;
    static double dxrel = 0.;

    /* System generated locals */
    int i__1;
    float r__1;
    double ret_val, d__1, d__2;

    /* Local variables */
    static int i__, n;
    static double y, temp, sinpiy;

/*     jan 1984 edition.  w. fullerton, c3, los alamos scientific lab. */
/*     double precision x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax, */
/*     1  xmin, y, d9lgmc, dcsevl, d1mach, dexp, dint, dlog, */
/*     2  dsin, dsqrt */
/*     external d1mach, d9lgmc, dcsevl, dexp, dint, dlog, dsin, dsqrt, */
/*     1  initds */

/*     series for gam        on the interval  0.          to  1.00000e+00 */
/*     with weighted error   5.79e-32 */
/*     log weighted error  31.24 */
/*     significant figures required  30.00 */
/*     decimal places required  32.05 */


/*     sq2pil is 0.5*alog(2*pi) = alog(sqrt(2*pi)) */
    ret_val = -999.;

    if (ngam == 0) {
/*        ngam = initds (gamcs, 42, 0.1*sngl(  d1mach) ) */
	r__1 = (float) machfd_1.epsmin * .1f;
	ngam = initds_(gamcs, &c__42, &r__1);

	d9gaml_(&xmin, &xmax);
	if (gammfd_1.igamma != 0) {
	    return ret_val;
	}
/*        xsml = dexp (dmax1 (dlog(d1mach(1)), -dlog(d1mach(2)))+0.01d0) */
/* Computing MAX */
	d__1 = log(machfd_1.fltmin), d__2 = -log(machfd_1.fltmax);
	xsml = exp(max(d__1,d__2) + .01);
/*        dxrel = dsqrt (d1mach(4)) */
	dxrel = sqrt(machfd_1.epsmax);

    }
/*     y = fabs(x) */
    y = abs(*x);
    if (y > 10.) {
	goto L50;
    }

/*     compute gamma(x) for -xbnd .le. x .le. xbnd.  reduce interval and find */
/*     gamma(1+y) for 0.0 .le. y .lt. 1.0 first of all. */

    n = (int) (*x);
    if (*x < 0.) {
	--n;
    }
    y = *x - (double) ((float) n);
    --n;
/*     dgamma = 0.9375d0 + dcsevl (2.d0*y-1.d0, gamcs, ngam) */
    d__1 = y * 2. - 1.;
    temp = dcsevl_(&d__1, gamcs, &ngam);
    if (gammfd_1.igamma != 0) {
	return ret_val;
    }
    ret_val = temp + .9375;
    if (n == 0) {
	return ret_val;
    }

    if (n > 0) {
	goto L30;
    }

/*     compute gamma(x) for x .lt. 1.0 */

    n = -n;
/*     if (x.eq.0.d0) call seteru (14hdgamma  x is 0, 14, 4, 2) */
/*     if (x.lt.0d0 .and. x+dble(float(n-2)).eq.0.d0) call seteru ( */
/*     1  31hdgamma  x is a negative integer, 31, 4, 2) */
/*     if (x.lt.(-0.5d0) .and. fabs((x-dint(x-0.5d0))/x).lt.dxrel) call */
/*     1  seteru (68hdgamma  answer lt half precision because x too near n */
/*     2egative integer, 68, 1, 1) */
/*     if (y.lt.xsml) call seteru ( */
/*     1  54hdgamma  x is so close to 0.0 that the result overflows, */
/*     2  54, 5, 2) */
    if (*x == 0.) {
/*     write(6,*) 'dgamma : x is 0' */
	gammfd_1.igamma = 11;
	return ret_val;
    }
    if (*x < 0. && *x + (double) ((float) (n - 2)) == 0.) {
/*     write( 6, *) 'dgamma : x is a negative integer' */
	gammfd_1.igamma = 12;
	return ret_val;
    }
    if (*x < -.5 && (d__1 = (*x - (double) ((int) (*x - .5))) / *x,
	    abs(d__1)) < dxrel) {
	gammfd_1.jgamma = 11;
    }
/*     1  write(6,*) 'dgamma : answer lt half precision because */
/*     2                       x too near a negative integer' */
    if (y < xsml) {
/*     write(6,*)  'dgamma :, */
/*     1               x is so close to 0.0 that the result overflows' */
	gammfd_1.igamma = 13;
	return ret_val;
    }

    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val /= *x + (double) ((float) (i__ - 1));
/* L20: */
    }
    return ret_val;

/*     gamma(x) for x .ge. 2.0 and x .le. 10.0 */

L30:
    i__1 = n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ret_val = (y + (double) ((float) i__)) * ret_val;
/* L40: */
    }
    return ret_val;

/*     gamma(x) for fabs(x) .gt. 10.0.  recall y = fabs(x). */

L50:
    if (*x > xmax) {
/*     write(6,*) 'dgamma : x so big gamma overflows' */
	gammfd_1.igamma = 14;
	return ret_val;
    }

    ret_val = 0.;
    if (*x < xmin) {
/*     write(6,*) 'dgamma : x so small gamma underflows' */
	gammfd_1.jgamma = 12;
	return ret_val;
    }

/*     dgamma = dexp ((y-0.5d0)*dlog(y) - y + sq2pil + d9lgmc(y) ) */
    temp = d9lgmc_(&y);
    if (gammfd_1.igamma != 0) {
	return ret_val;
    }
    ret_val = exp((y - .5) * log(y) - y + sq2pil + temp);
    if (*x > 0.) {
	return ret_val;
    }

/*     if (fabs((x-dint(x-0.5d0))/x).lt.dxrel) call seteru ( */
/*     1  61hdgamma  answer lt half precision, x too near negative integer */
/*     2  , 61, 1, 1) */
    if ((d__1 = (*x - (double) ((int) (*x - .5))) / *x, abs(d__1)) <
	    dxrel) {
	gammfd_1.jgamma = 11;
    }

/*     sinpiy = dsin (pi*y) */
    sinpiy = sin(pi * y);
    if (sinpiy == 0.) {
/*     write(6,*) 'dgamma : x is a negative integer' */
	gammfd_1.igamma = 12;
	return ret_val;
    }

    ret_val = -pi / (y * sinpiy * ret_val);
    return ret_val;
} /* dgamma_ */

double dgamr_(double *x)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    static double temp, alngx, sgngx;

/*     july 1977 edition.  w. fullerton, c3, los alamos scientific lab. */
/*     this routine, not dgamma(x), should be the fundamental one. */
/*     ============  ============= */
/* Calls dgamma(), only if |x| < 10; otherwise dlgams() -> dlngam() -> d9lgmc() */
/*     external dexp, dgamma, dint, d1mach */

    ret_val = 0.;
    if (*x <= 0. && (double) ((int) (*x)) == *x) {
	return ret_val;
    }

    if (abs(*x) <= 10.) {
/*     dgamr = 1.0d0/dgamma(x) */
	temp = dgamma_(x);
	if (gammfd_1.igamma != 0) {
	    ret_val = machfd_1.fltmax;
	    return ret_val;
	}
	ret_val = 1. / temp;
    } else {
/*     x > 10. : */
	dlgams_(x, &alngx, &sgngx);
	if (gammfd_1.igamma != 0) {
	    return ret_val;
	}
	ret_val = sgngx * exp(-alngx);
    }
    return ret_val;
} /* dgamr_ */

/* Subroutine */
int dlgams_(double *x, double *dlgam, double *sgngam)
{
    /* System generated locals */
    double d__1;

    /* Local variables */
    int intx;

/* july 1977 edition.  w. fullerton, c3, los alamos scientific lab. */

/* evaluate log abs (gamma(x)) and return the sign of gamma(x) in sgngam. */
/* sgngam is either +1.0 or -1.0. */


    *dlgam = dlngam_(x);
    if (gammfd_1.igamma != 0) {
	return 0;
    }
    *sgngam = 1.;
    if (*x > 0.) {
	return 0;
    }

    d__1 = -((double) ((int) (*x)));
    intx = (int) (fmod(d__1, c_b9) + .1);
    if (intx == 0) {
	*sgngam = -1.;
    }

    return 0;
} /* dlgams_ */

int initds_(double *dos, int *nos, float *eta)
{
    /* System generated locals */
    int ret_val, i__1;
    float r__1;

    /* Local variables */
    static int i__, ii;
    static double err;

/* june 1977 edition.   w. fullerton, c3, los alamos scientific lab. */

/* initialize the double precision orthogonal series dos so that initds */
/* is the number of terms needed to insure the error is no larger than */
/* eta.  ordinarily eta will be chosen to be one-tenth machine precision. */

/*             input arguments -- */
/* dos    dble prec array of nos coefficients in an orthogonal series. */
/* nos    number of coefficients in dos. */
/* eta    requested accuracy of series. */


/*     if (nos.lt.1) call seteru ( */
/*    1  35hinitds  number of coefficients lt 1, 35, 2, 2) */
    /* Parameter adjustments */
    --dos;

    /* Function Body */
    if (*nos < 1) {
	gammfd_1.jgamma = 31;
    }

    i__ = -1;
    err = 0.f;
    i__1 = *nos;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *nos + 1 - ii;
	err += (r__1 = (float) dos[i__], fabs(r__1));
	if (err > *eta) {
	    goto L20;
	}
/* L10: */
    }

/* 20   if (i.eq.nos) call seteru (28hinitds  eta may be too small, 28, */
/*    1  1, 2) */
L20:
/*     if (i.eq.nos) write(6,*) 'initds : eta may be too small' */
    if (i__ == *nos) {
	gammfd_1.jgamma = 32;
    }
    ret_val = i__;

    return ret_val;
} /* initds_ */

/* Subroutine */
static void d9gaml_(double *xmin, double *xmax)
{
    /* System generated locals */
    double d__1, d__2;

    /* Local variables */
    static int i__;
    static double xln, xold, alnbig, alnsml;

/* june 1977 edition.   w. fullerton, c3, los alamos scientific lab. */

/* calculate the minimum and maximum legal bounds for x in gamma(x). */
/* xmin and xmax are not the only bounds, but they are the only non- */
/* trivial ones to calculate. */

/*             output arguments -- */
/* xmin   dble prec minimum legal value of x in gamma(x).  any smaller */
/*        value of x might result in underflow. */
/* xmax   dble prec maximum legal value of x in gamma(x).  any larger */
/*        value of x might cause overflow. */

/*     double precision xmin, xmax, alnbig, alnsml, xln, xold, d1mach, */
/*    1  dlog */
/*     external d1mach, dlog */

/*     alnsml = dlog(d1mach(1)) */
    alnsml = log(machfd_1.fltmin);
    *xmin = -alnsml;
    for (i__ = 1; i__ <= 10; ++i__) {
	xold = *xmin;
/*       xln = dlog(xmin) */
	xln = log(*xmin);
	*xmin -= *xmin * ((*xmin + .5) * xln - *xmin - .2258 + alnsml) / (*
		xmin * xln + .5);
/*       if (fabs(xmin-xold).lt.0.005d0) go to 20 */
	if ((d__1 = *xmin - xold, abs(d__1)) < .005) {
	    goto L20;
	}
/* L10: */
    }
/*     call seteru (27hd9gaml  unable to find xmin, 27, 1, 2) */
/*     write(6,*) 'd9gaml : unable to find xmin' */
    gammfd_1.igamma = 21;
    return;

L20:
    *xmin = -(*xmin) + .01;

/*     alnbig = dlog (d1mach(2)) */
    alnbig = log(machfd_1.fltmax);
    *xmax = alnbig;
    for (i__ = 1; i__ <= 10; ++i__) {
	xold = *xmax;
/*       xln = dlog(xmax) */
	xln = log(*xmax);
	*xmax -= *xmax * ((*xmax - .5) * xln - *xmax + .9189 - alnbig) / (*
		xmax * xln - .5);
/*       if (fabs(xmax-xold).lt.0.005d0) go to 40 */
	if ((d__1 = *xmax - xold, abs(d__1)) < .005) {
	    goto L40;
	}
/* L30: */
    }
/*     call seteru (27hd9gaml  unable to find xmax, 27, 2, 2) */
/*     write(6,*) 'd9gaml : unable to find xmax' */
    gammfd_1.igamma = 22;
    return;

L40:
    *xmax += -.01;
/* Computing MAX */
    d__1 = *xmin, d__2 = -(*xmax) + 1.;
    *xmin = max(d__1,d__2);

    return;

} /* d9gaml_ */

double d9lgmc_(double *x)
{
    /* Initialized data */

    static double algmcs[15] = { .1666389480451863247205729650822,
	    -1.384948176067563840732986059135e-5,
	    9.810825646924729426157171547487e-9,
	    -1.809129475572494194263306266719e-11,
	    6.221098041892605227126015543416e-14,
	    -3.399615005417721944303330599666e-16,
	    2.683181998482698748957538846666e-18,
	    -2.868042435334643284144622399999e-20,
	    3.962837061046434803679306666666e-22,
	    -6.831888753985766870111999999999e-24,
	    1.429227355942498147573333333333e-25,
	    -3.547598158101070547199999999999e-27,1.025680058010470912e-28,
	    -3.401102254316748799999999999999e-30,
	    1.276642195630062933333333333333e-31 };
    static int nalgm = 0;
    static double xbig = 0.;
    static double xmax = 0.;

    /* System generated locals */
    float r__1;
    double ret_val, d__1, d__2;

    /* Local variables */
    static double temp;

/* august 1977 edition.  w. fullerton, c3, los alamos scientific lab. */

/* compute the log gamma correction factor for x .ge. 10. so that */
/* dlog (dgamma(x)) = dlog(dsqrt(2*pi)) + (x-.5)*dlog(x) - x + d9lgmc(x) */

/*     double precision x, algmcs(15), xbig, xmax, dcsevl, d1mach, */
/*    1  dexp, dlog, dsqrt */
/*     external d1mach, dcsevl, dexp, dlog, dsqrt, initds */

/* series for algm       on the interval  0.          to  1.00000e-02 */
/*                                        with weighted error   1.28e-31 */
/*                                         log weighted error  30.89 */
/*                               significant figures required  29.81 */
/*                                    decimal places required  31.48 */



    if (nalgm != 0) {
	goto L10;
    }
/*     nalgm = initds (algmcs, 15, sngl(d1mach(3)) ) */
    r__1 = (float) machfd_1.epsmin;
    nalgm = initds_(algmcs, &c__15, &r__1);
/*     xbig = 1.0d0/dsqrt(d1mach(3)) */
    xbig = 1. / sqrt(machfd_1.epsmin);
/*     xmax = dexp (dmin1(dlog(d1mach(2)/12.d0), -dlog(12.d0*d1mach(1)))) */
/* Computing MIN */
    d__1 = log(machfd_1.fltmax / 12.), d__2 = -log(machfd_1.fltmin * 12.);
    xmax = exp((min(d__1,d__2)));

/* 10   if (x.lt.10.d0) call seteru (23hd9lgmc  x must be ge 10, 23, 1, 2) */

L10:
    if (*x < 10.) {
/*       write(6,*) 'd9lgmc : x must be ge 10' */
	gammfd_1.igamma = 51;
/*       d9lgmc = d1mach(2) */
	ret_val = machfd_1.fltmax;
	return ret_val;
    }
    if (*x >= xmax) {
	goto L20;
    }

    ret_val = 1. / (*x * 12.);
/*     if (x.lt.xbig) d9lgmc = dcsevl (2.0d0*(10.d0/x)**2-1.d0, algmcs, */
/*    1  nalgm) / x */
    if (*x < xbig) {
/* Computing 2nd power */
	d__2 = 10. / *x;
	d__1 = d__2 * d__2 * 2. - 1.;
	temp = dcsevl_(&d__1, algmcs, &nalgm);
	if (gammfd_1.igamma != 0) {
/*         d9lgmc = d1mach(2) */
	    ret_val = machfd_1.fltmax;
	} else {
	    ret_val = temp / *x;
	}
    }
    return ret_val;

L20:
    ret_val = 0.;
/*     call seteru (34hd9lgmc  x so big d9lgmc underflows, 34, 2, 0) */
/*     write(6,*) 'd9lgmc : x so big d9lgmc underflows' */
    gammfd_1.jgamma = 51;
    return ret_val;

} /* d9lgmc_ */

double dcsevl_(double *x, double *a, int *n)
{
    /* System generated locals */
    int i__1;
    double ret_val;

    /* Local variables */
    static int i__;
    static double b0, b1, b2;
    static int ni;
    static double twox;


/* evaluate the n-term chebyshev series a at x.  adapted from */
/* r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973). */

/*             input arguments -- */
/* x      dble prec value at which the series is to be evaluated. */
/* a      dble prec array of n terms of a chebyshev series.  in eval- */
/*        uating a, only half the first coef is summed. */
/* n      number of terms in array a. */

/*     double precision d1mach */
/*     external         d1mach */

    /* Parameter adjustments */
    --a;

    /* Function Body */
    b2 = 0.f;
/*     if (n.lt.1) call seteru (28hdcsevl  number of terms le 0, 28, 2,2) */
/*     if (n.gt.1000) call seteru (31hdcsevl  number of terms gt 1000, */
/*    1  31, 3, 2) */
/*     if (x.lt.(-1.1d0) .or. x.gt.1.1d0) call seteru ( */
/*    1  25hdcsevl  x outside (-1,+1), 25, 1, 1) */

    if (*n < 1) {
/*       write(6,*) 'dcsevl : number of terms le 0' */
	gammfd_1.igamma = 41;
/*       dcsevl = d1mach(2) */
	ret_val = machfd_1.fltmax;
	return ret_val;
    }
    if (*n > 1000) {
/*       write(6,*) 'dcsevl : number of terms gt 1000' */
	gammfd_1.igamma = 42;
/*       dcsevl = d1mach(2) */
	ret_val = machfd_1.fltmax;
	return ret_val;
    }
    if (*x < -1.1 || *x > 1.1) {
/*       write(6,*) 'dcsevl : x outside (-1,+1)' */
	gammfd_1.igamma = 43;
/*       dcsevl = d1mach(2) */
	ret_val = machfd_1.fltmax;
	return ret_val;
    }

    twox = *x * 2.;
    b1 = 0.;
    b0 = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	b2 = b1;
	b1 = b0;
	ni = *n - i__ + 1;
	b0 = twox * b1 - b2 + a[ni];
/* L10: */
    }

    ret_val = (b0 - b2) * .5;

    return ret_val;
} /* dcsevl_ */

double dlngam_(double *x)
{
    /* Initialized data */

    static double sq2pil = .91893853320467274178032973640562;
    static double sqpi2l = .225791352644727432363097614947441;
    static double pi = 3.1415926535897932384626433832795;
    static double xmax = 0.;
    static double dxrel = 0.;

    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    static double y, temp, sinpiy;

/*     august 1980 edition.   w. fullerton, c3, los alamos scientific lab. */
/*     double precision x, dxrel, pi, sinpiy, sqpi2l, sq2pil, */
/*     1  y, xmax, dint, dgamma, d9lgmc, d1mach, dlog, dsin, dsqrt */
/*     external d1mach, d9lgmc, dgamma, dint, dlog, dsin, dsqrt */

/*     sq2pil = alog (sqrt(2*pi)),  sqpi2l = alog(sqrt(pi/2)) */


    ret_val = 0.;
    if (xmax == 0.) {
/*     xmax = d1mach(2)/dlog(d1mach(2)) */
	xmax = machfd_1.fltmax / log(machfd_1.fltmax);
/*     dxrel = dsqrt (d1mach(4)) */
	dxrel = sqrt(machfd_1.fltmax);
    }
    y = abs(*x);
    if (y <= 10.) {

/*     |x| <= 10 :  Compute  dlngam := dlog (fabs (dgamma(x)) ) */

	temp = dgamma_(x);
	if (gammfd_1.igamma != 0) {
	    ret_val = machfd_1.fltmax;
	    return ret_val;
	}
	ret_val = log((abs(temp)));
	return ret_val;
    }
/*     ELSE  |x| > 10 :  Compute dlog ( fabs (dgamma(x)) ) */

    if (y > xmax) {
/*     write(6,*) 'dlngam : abs(x) so big dlngam overflows' */
	gammfd_1.igamma = 61;
	ret_val = machfd_1.fltmax;
	return ret_val;
    }

/*     if (x.gt.0.d0) dlngam = sq2pil + (x-0.5d0)*dlog(x) - x + d9lgmc(y) */
    temp = d9lgmc_(&y);
    if (gammfd_1.igamma != 0) {
	ret_val = machfd_1.fltmax;
	return ret_val;
    }
    if (*x > 0.) {
	ret_val = sq2pil + (*x - .5) * log(*x) - *x + temp;
    }
    if (*x > 0.) {
	return ret_val;
    }

    sinpiy = (d__1 = sin(pi * y), abs(d__1));
    if (sinpiy == 0.) {
/*     write(6,*) 'dlngam : x is a negative integer' */
	gammfd_1.igamma = 62;
	ret_val = machfd_1.fltmax;
	return ret_val;
    }

/*     dlngam = sqpi2l + (x-0.5d0)*dlog(y) - x - dlog(sinpiy) - d9lgmc(y) */
    temp = d9lgmc_(&y);
    if (gammfd_1.igamma != 0) {
	ret_val = machfd_1.fltmax;
	return ret_val;
    }
    ret_val = sqpi2l + (*x - .5) * log(y) - *x - log(sinpiy) - temp;

/*     if (fabs((x-dint(x-0.5d0))*dlngam/x).lt.dxrel) call seteru ( */
/*     1  68hdlngam  answer lt half precision because x too near negative */
/*     2integer, 68, 1, 1) */
    if ((d__1 = (*x - (double) ((int) (*x - .5))) * ret_val / *x, abs(
	    d__1)) < dxrel) {
	gammfd_1.jgamma = 61;
    }
    return ret_val;

} /* dlngam_ */

