/* fdcore.f -- translated by f2c (version 20031025).
 *
 * and produced by
 * $Id: f2c-clean,v 1.10 2002/03/28 16:37:27 maechler Exp $
 *
 * and manually pretty edited by Martin Maechler, 2004-09-18
 */

#include <Rmath.h>

/* dcopy() and ddot() only:*/
#include <R_ext/BLAS.h>

#include "fracdiff.h"

typedef int /* Unknown procedure type */ (*U_fp)();

extern double dgamr_(double *);
extern double dgamma_(double *);

/* Subroutine */
void fracdf_(double *x, int *n, int *m, int *nar, int *nma,
	     double *dtol, double *drange, double *hood,
	     double *d__, double *ar, double *ma, double *w,
	     int *lenw, int *iw, int *inform__,
	     double *flmin, double *flmax, double *epmin, double *epmax);

static
double dopt(double *x, double *dinit, double *drange,
	    double *hood, double *delta, double *w, int *iw);

static
double pqopt_(double *x, double *d__, double *w, int *iw);


/* These + ajqp_(..)  are passed to LMDER1() to be optimized: */
static
int ajp_(double *p, double *a, double *ajac,
	 int *lajac, int *iflag, double *y);
static
int ajq_(double *qp, double *a, double *ajac,
	 int *lajac, int *iflag, double *y);

/* Common Block Declarations */

struct {
    double fltmin, fltmax, epsmin, epsmax;
} machfd_;

#define machfd_1 machfd_

struct {
    double epsp25, epspt3, epspt5, epsp75, bignum;
} mauxfd_;

#define mauxfd_1 mauxfd_

union {
    struct { int nn, mm, np, nq, npq, npq1, maxpq, maxpq1, minpq, nm; } _1;
    struct { int n,  m,  np, nq, npq, npq1, maxpq, maxpq1, minpq, nm; } _2;
} dimsfd_;

#define dimsfd_1 (dimsfd_._1)
#define dimsfd_2 (dimsfd_._2)

struct {
    int maxopt, maxfun, nopt, nfun, ngrd, ifun, igrd, info;
} cntrfd_;

#define cntrfd_1 cntrfd_

union {
    struct {
	double told, tolf, tolx, tolg, anorm, deltax, gnorm;
    } _1;
    struct {
	double dtol, ftol, xtol, gtol, anorm, deltax, gnorm;
    } _2;
} tolsfd_;

#define tolsfd_1 (tolsfd_._1)
#define tolsfd_2 (tolsfd_._2)

struct {
    int ly, lamk, lak, lvk, lphi, lpi;
} wfilfd_;

#define wfilfd_1 wfilfd_

struct {
    int lqp, la, lajac, ipvt, ldiag, lqtf, lwa1, lwa2, lwa3, lwa4;
} woptfd_;

#define woptfd_1 woptfd_

struct {
    int ilimit, jlimit;
} limsfd_;

#define limsfd_1 limsfd_

struct {
    int igamma, jgamma;
} gammfd_;

#define gammfd_1 gammfd_

struct {
    int iminpk, jminpk;
} mnpkfd_;

#define mnpkfd_1 mnpkfd_

struct {
    int ksvd, kcov, kcor;
} hessfd_;

#define hessfd_1 hessfd_

union {
    struct {
	double hatmu, wnv, cllf;
    } _1;
    struct {
	double hatmu, wnv, hood;
    } _2;
} filtfd_;

#define filtfd_1 (filtfd_._1)
#define filtfd_2 (filtfd_._2)

/* Table of constant values */

static double c_b2 = -99.f;
static int ic__1 = 1;
static int ic__0 = 0;
static double c__1 = 1.;

/*****************************************************************************
 ******************************************************************************/
void fracdf_(double *x, int *n, int *m, int *nar, int *nma,
	     double *dtol, double *drange, double *hood,
	     double *d__, double *ar, double *ma, double *w,
	     int *lenw, int *iw, int *inform__,
	     double *flmin, double *flmax, double *epmin, double *epmax)
{
/*   float              xa(n)
     double precision   ar(*), ma(*), drange(2)
     double precision   w(*)
*/

    /* System generated locals */
    double d__1;

    /* Local variables */
    double delta;
    int lfree, lwfree, lenthw;

/* ------------------------------------------------------------------------------

   Input :

  x       double   time series for the ARIMA model
  n       int  length of the time series
  M       int  number of terms in the likelihood approximation
                   suggested value 100 (see Haslett and Raftery 1989)
  nar     int  number of autoregressive parameters
  nma     int  number of moving average parameters
  dtol    double   desired length of final interval of uncertainty for d
                   suggested value : 4th root of machine precision
                   if dtol < 0 it is automatically set to this value
                   dtol will be altered if necessary by the program
  drange  double   array of length 2 giving minimum and maximum values f
                   for the fractional differencing parameter
  d       double   initial guess for optimal fractional differencing parameter
  w       double   work array
  lenw    int  length of double precision workspace w, must be at least
 		max( p+q+2*(n+M), 3*n+(n+6.5)*(p+q) +1, (3+2*(p+q+1))*(p+q+1)+1)
   MM:		max( p+q+2*(n+M), 3*n+(n+6.5)*(p+q) +1,     31 * 12)
 	 is what the code below rather checks

  Output :

  dtol    double   value of dtol ultimately used by the algorithm
  d       double   final value optimal fractional differencing parameter
  hood    double   logarithm of the maximum likelihood
  ar      double   optimal autoregressive parameters
  ma      double   optimal moving average parameters

 ------------------------------------------------------------------------------
  copyright 1991 Department of Statistics, University of Washington
  written by Chris Fraley
 -----------------------------------------------------------------------------
     Parameter adjustments */
    --ar;
    --ma;
    --drange;
    --w;

    /* Function Body */
    if (*m <= 0) {
	*m = 100;
    }
/*     MM: Using 'fdcom' instead of 'code copy' -- FIXME: use #include in C
     initialize several of the above common blocks: */
    fdcom(n, m, nar, nma, &c_b2, flmin, flmax, epmin, epmax);
    lfree = woptfd_1.lwa4 + *n - dimsfd_1.minpq;
/* 	= 1+ ipvt + 5.5*npq + n - minpq
 	= 2+ 6.5*npq + 3*n - 2*minpq + (n-maxpq)*npq
 and               lvk+M = 1 + npq + 2(n + M)
 */

    lwfree = imax2(372, imax2(wfilfd_1.lvk + *m, lfree));
/*                                 ^^^^^^^ MM: where is this needed? */
    if (lwfree > *lenw + 1) {
	limsfd_1.ilimit = lwfree - *lenw;
/*       write( 6, *) 'insufficient storage : ',
    *               'increase length of w by at least', ILIMIT */
	*inform__ = 1;
/* 	return the *desired* workspace storage: */
	*lenw = lwfree;
	return;
    }
    lenthw = *lenw;
    cntrfd_1.maxopt = 100;
    cntrfd_1.maxfun = 100;
/* set error and warning flags */
    *inform__ = 0;
    gammfd_1.igamma = 0;
    mnpkfd_1.iminpk = 0;
    limsfd_1.ilimit = 0;
    gammfd_1.jgamma = 0;
    mnpkfd_1.jminpk = 0;
    limsfd_1.jlimit = 0;
    if (*dtol > .1) {
	*dtol = .1;
    }
    if (*dtol <= 0.) {
	tolsfd_1.told = mauxfd_1.epsp25;
	tolsfd_1.tolf = mauxfd_1.epspt3;
    } else {
	tolsfd_1.told = fmax2(*dtol,mauxfd_1.epspt5);
/* Computing MAX */
	d__1 = *dtol / 10.;
	tolsfd_1.tolf = fmax2(d__1,mauxfd_1.epsp75);
    }
    tolsfd_1.tolg = tolsfd_1.tolf;
    tolsfd_1.tolx = tolsfd_1.told;
    *dtol = tolsfd_1.told;
/*     if (npq != 0) call dcopy( npq, zero, 0, w(lqp), 1) */
    if (dimsfd_1.npq != 0) {
	F77_CALL(dcopy)(&dimsfd_1.np, &ar[1], &ic__1, &w[woptfd_1.lqp + dimsfd_1.nq], &
		ic__1);
	F77_CALL(dcopy)(&dimsfd_1.nq, &ma[1], &ic__1, &w[woptfd_1.lqp], &ic__1);
    }
    cntrfd_1.nopt = 0;
    cntrfd_1.nfun = 0;
    cntrfd_1.ngrd = 0;
/* 	   ==== */
    *d__ = dopt(x, d__, &drange[1], hood, &delta, &w[1], iw);
/* 	   ==== */
    if (cntrfd_1.nopt >= cntrfd_1.maxopt) {
	limsfd_1.jlimit = 1;
    }
/*       write( 6, *)
       write( 6, *) 'WARNING : optimization limit reached'
     end if */
    if (gammfd_1.igamma != 0 || mnpkfd_1.iminpk != 0) {
	*d__ = machfd_1.fltmax;
	*hood = machfd_1.fltmax;
	F77_CALL(dcopy)(&dimsfd_1.np, &machfd_1.fltmax, &ic__0, &ar[1], &ic__1);
	F77_CALL(dcopy)(&dimsfd_1.nq, &machfd_1.fltmax, &ic__0, &ma[1], &ic__1);
	if (gammfd_1.igamma != 0) {
	    *inform__ = 2;
	}
	if (mnpkfd_1.iminpk != 0) {
	    *inform__ = 3;
	}
	return;
    }
    F77_CALL(dcopy)(&dimsfd_1.np, &w[woptfd_1.lqp + dimsfd_1.nq], &ic__1, &ar[1], &ic__1)
	    ;
    F77_CALL(dcopy)(&dimsfd_1.nq, &w[woptfd_1.lqp], &ic__1, &ma[1], &ic__1);
    if (gammfd_1.jgamma != 0) {
	*inform__ = 4;
    }
    if (mnpkfd_1.jminpk != 0) {
	*inform__ = 5;
    }
    if (limsfd_1.jlimit != 0) {
	*inform__ = 6;
    }
    return;
/* 900  format( 4h itr, 14h     d          ,   14h    est mean  ,
     *                16h     white noise,  17h     log likelihd,
     *                 4h  nf, 3h ng) */

} /* fracdf() {main} */


/******************************************************************************
 *****************************************************************************

 optimization with respect to d based on Brent's fmin algorithm */

double dopt(double *x, double *dinit, double *drange,
	    double *hood, double *delta, double *w, int *iw)
{
/*     float              x(n) */

    /* Initialized data */

/*  cc is the squared inverse of the golden ratio:
    cc = half*(three-sqrt(5.0d0))
*/
    static double cc = .38196601125011;

    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    static double d__, aa, bb, fa, dd, fb, ee, hh, fu, fv, fw, fx, rr, ss,
	     tt, uu, vv, ww, xx, eps, tol, tol1, tol2, tol3;

/*  copyright 1991 Department of Statistics, University of Washington
  written by Chris Fraley
 ------------------------------------------------------------------------------
*/


/* eps is approximately the square root of the relative machine
  precision. */
    eps = machfd_1.epsmax;
    tol1 = eps + 1.;
    eps = sqrt(eps);
/* -Wall: */
    ret_val = -1.;
    dd = 0.;

    aa = drange[0];
    bb = drange[1];
    if (*dinit > aa + tolsfd_2.dtol &&
	*dinit < bb - tolsfd_2.dtol) {
	vv = *dinit;
    } else {
	vv = aa + cc * (bb - aa);
    }
    ww = vv;
    xx = vv;
    uu = xx;
    ee = 0.;
    cntrfd_1.nopt = 1;
    fx = pqopt_(x, &xx, w, iw);
/*       ===== */
    fv = fx;
    fw = fx;
    tol = fmax2(tolsfd_2.dtol,0.);
    tol3 = tol / 3.;

/*  main loop starts here */

L10:
    if (gammfd_1.igamma != 0 || mnpkfd_1.iminpk != 0) {
	d__ = uu;
	*hood = machfd_1.fltmax;
	return ret_val;
    }
    hh = (aa + bb) * .5;
    tol1 = eps * (fabs(xx) + 1.) + tol3;
    tol2 = tol1 * 2.;

/*  check stopping criterion */

    *delta = (d__1 = xx - hh, fabs(d__1)) + (bb - aa) * .5;
/*     if (abs(xx-hh) .le. (tol2-half*(bb-aa))) goto 100 */
    if (*delta <= tol2) {
	goto L100;
    }
    if (cntrfd_1.nopt >= cntrfd_1.maxopt) {
	goto L100;
    }
/*     if (delpq <= EPSMAX*(one+pqnorm)) goto 100 */
    rr = 0.;
    ss = 0.;
    tt = 0.;
    if (fabs(ee) > tol1) {

/*  fit parabola */

	rr = (xx - ww) * (fx - fv);
	ss = (xx - vv) * (fx - fw);
	tt = (xx - vv) * ss - (xx - ww) * rr;
	ss = (ss - rr) * 2.;
	if (ss <= 0.) {
	    ss = -ss;
	} else {
	    tt = -tt;
	}
	rr = ee;
	ee = dd;
    }
    if (fabs(tt) >= (d__1 = ss * .5 * rr, fabs(d__1)) || tt <= ss * (aa - xx) ||
	     tt >= ss * (bb - xx)) {

/*  a golden-section step */

	if (xx >= hh) {
	    ee = aa - xx;
	} else {
	    ee = bb - xx;
	}
	dd = cc * ee;
    } else {

/*  a parabolic-interpolation step */

	dd = tt / ss;
	uu = xx + dd;

/*  f must not be evaluated too close to aa or bb */

	if (uu - aa < tol2 || bb - uu < tol2) {
	    dd = tol1;
	    if (xx >= hh) {
		dd = -dd;
	    }
	}
    }

/*  f must not be evaluated too close to xx */

    if (fabs(dd) >= tol1) {
	uu = xx + dd;
    } else {
	if (dd <= 0.) {
	    uu = xx - tol1;
	} else {
	    uu = xx + tol1;
	}
    }
    ++cntrfd_1.nopt;
    fu = pqopt_(x, &uu, w, iw);

/*  update  aa, bb, vv, ww, and xx */

    if (fx >= fu) {
	if (uu >= xx) {
	    aa = xx;
	    fa = fx;
	} else {
	    bb = xx;
	    fb = fx;
	}
	vv = ww;
	fv = fw;
	ww = xx;
	fw = fx;
	xx = uu;
	fx = fu;
    } else {
	if (uu >= xx) {
	    bb = uu;
	    fb = fu;
	} else {
	    aa = uu;
	    fa = fu;
	}
	if (fu > fw && ww != xx) {
	    if (fu <= fv || vv == xx || vv == ww) {
		vv = uu;
		fv = fu;
	    }
	} else {
	    vv = ww;
	    fv = fw;
	    ww = uu;
	    fw = fu;
	}
    }
    goto L10;

/*  end of main loop */

L100:
    ret_val = xx;
    *hood = -fx;
    filtfd_1.cllf = *hood;
    return ret_val;
/* 900  format( i4, 2(1pe14.6), 1pe16.7, 1pe17.8, 1x, 2(i3))
 901  format( i4, 3(1pe10.2), 1pe11.2, 2(i3), 3(1pe8.1), i2) */
} /* dopt */


/**************************************************************************
 ************************************************************************** */
static
double pqopt_(double *x, double *d__, double *w, int *iw)
{
/* x: double x(n) */
/* w: work array exactly as in main  fracdf() */

    /* Initialized data */

    static int modelm = 1;
    static double factlm = 100.;

    /* System generated locals */
    int i__1, i__2;
    double ret_val;

    /* Local variables */
    static double t, u, bic, slogvk;

    /* Parameter adjustments */
    --w;

    /* copyright 1991 Department of Statistics, University of Washington
     * written by Chris Fraley
 ---------------------------------------------------------------------------- */
    fdfilt_(x, d__,
	    &w[(0 + (0 + (wfilfd_1.ly   << 3))) / 8], &slogvk,
	    &w[(0 + (0 + (wfilfd_1.lamk << 3))) / 8],
	    &w[(0 + (0 + (wfilfd_1.lak  << 3))) / 8],
	    &w[(0 + (0 + (wfilfd_1.lvk  << 3))) / 8],
	    &w[(0 + (0 + (wfilfd_1.lphi << 3))) / 8],
	    &w[(0 + (0 + (wfilfd_1.lpi  << 3))) / 8]);
    if (gammfd_1.igamma != 0) {
	ret_val = machfd_1.fltmax;
	filtfd_2.wnv = machfd_1.fltmax;
	filtfd_2.hood = -machfd_1.fltmax;
	return ret_val;
    }
    t = (double) dimsfd_2.n;
    if (dimsfd_2.npq == 0) {
/* 	trivial case  p = q = 0 : */
	filtfd_2.wnv = F77_CALL(ddot)(&dimsfd_2.n,
				      &w[wfilfd_1.ly], &ic__1,
				      &w[wfilfd_1.ly], &ic__1) / t;
	cntrfd_1.ifun = 0;
	cntrfd_1.igrd = 0;
	cntrfd_1.info = -1;
    } else {

/*     optimize as an unconstrained optimization problem */

	if (modelm == 2) {
	    F77_CALL(dcopy)(&dimsfd_2.npq, &c__1, &ic__0,
			    &w[woptfd_1.ldiag], &ic__1);
	}
	if (cntrfd_1.nopt < 0) {
	    if (dimsfd_2.np != 0) {
		i__1 = dimsfd_2.n - dimsfd_2.np;
		i__2 = dimsfd_2.n - dimsfd_2.np;
		lmder1_((U_fp)ajp_, &i__1, &dimsfd_2.np,
			&w[woptfd_1.lqp + dimsfd_2.nq], &w[woptfd_1.la],
			&w[woptfd_1.lajac], &i__2, &tolsfd_2.ftol, &tolsfd_2.xtol,
			&tolsfd_2.gtol, &cntrfd_1.maxfun, &w[woptfd_1.ldiag],
			&modelm, &factlm, &cntrfd_1.info, &cntrfd_1.ifun,
			&cntrfd_1.igrd, iw /* was &w[woptfd_1.ipvt] */,
			&w[woptfd_1.lqtf],
			&w[woptfd_1.lwa1], &w[woptfd_1.lwa2], &w[woptfd_1.lwa3],
			&w[woptfd_1.lwa4], &w[wfilfd_1.ly]);
	    }
	    if (dimsfd_2.nq != 0) {
		i__1 = dimsfd_2.n - dimsfd_2.nq;
		i__2 = dimsfd_2.n - dimsfd_2.nq;
		lmder1_((U_fp)ajq_, &i__1, &dimsfd_2.nq, &w[woptfd_1.lqp],
			&w[woptfd_1.la], &w[woptfd_1.lajac], &i__2,
			&tolsfd_2.ftol, &tolsfd_2.xtol, &tolsfd_2.gtol,
			&cntrfd_1.maxfun, &w[woptfd_1.ldiag], &modelm, &factlm,
			&cntrfd_1.info, &cntrfd_1.ifun, &cntrfd_1.igrd,
			iw /* was &w[woptfd_1.ipvt] */, &w[woptfd_1.lqtf],
			&w[woptfd_1.lwa1], &w[woptfd_1.lwa2],
			&w[woptfd_1.lwa3], &w[woptfd_1.lwa4],
			&w[wfilfd_1.ly]);
	    }
	}
	lmder1_((U_fp)ajqp_, &dimsfd_2.nm, &dimsfd_2.npq, &w[woptfd_1.lqp],
		&w[woptfd_1.la], &w[woptfd_1.lajac], &dimsfd_2.nm,
		&tolsfd_2.ftol, &tolsfd_2.xtol, &tolsfd_2.gtol,
		&cntrfd_1.maxfun, &w[woptfd_1.ldiag], &modelm, &factlm,
		&cntrfd_1.info, &cntrfd_1.ifun, &cntrfd_1.igrd,
		iw /* was &w[woptfd_1.ipvt] */, &w[woptfd_1.lqtf],
		&w[woptfd_1.lwa1], &w[woptfd_1.lwa2],
		&w[woptfd_1.lwa3], &w[woptfd_1.lwa4],
		&w[wfilfd_1.ly]);
	if (cntrfd_1.info == 0) {
/*     write( 6, *) 'MINPACK : improper input parameters */
	    mnpkfd_1.iminpk = 10;
	    ret_val = machfd_1.fltmax;
	    filtfd_2.wnv = machfd_1.fltmax;
	    filtfd_2.hood = -machfd_1.fltmax;
	    return ret_val;
	}
	if (cntrfd_1.info == 5) {
/*     write( 6, *) 'MINPACK : function evaluation limit reached' */
	    mnpkfd_1.jminpk = 5;
	}
	if (cntrfd_1.info == 6) {
/*     write( 6, *) 'MINPACK : ftol is too small' */
	    mnpkfd_1.jminpk = 6;
	}
	if (cntrfd_1.info == 7) {
/*     write( 6, *) 'MINPACK : xtol is too small' */
	    mnpkfd_1.jminpk = 7;
	}
	if (cntrfd_1.info == 8) {
/*     write( 6, *) 'MINPACK : gtol is too small' */
	    mnpkfd_1.jminpk = 8;
	}
/*     call daxpy( npq, (-one), w(lpq), 1, w(lqp), 1
     delpq  = sqrt(ddot( npq, w(lqp), 1, w(lqp), 1))
     pqnorm = sqrt(ddot( npq, w(lpq), 1, w(lpq), 1)) */
	filtfd_2.wnv = tolsfd_2.anorm * tolsfd_2.anorm / (double) (
		dimsfd_2.nm - 1);
    }
    u = t * (log(filtfd_2.wnv) + 2.8378) + slogvk;
    ret_val = u / 2.;
    bic = u + (double) (dimsfd_2.np + dimsfd_2.nq + 1) * log(t);
    filtfd_2.hood = -ret_val;
    return ret_val;
} /* End pqopt_() */

/*************************************************************************** */

/* Subroutine */ int
fdfilt_(double *x, double *d__, double *y,
	double *slogvk, double *amk, double *ak, double *vk,
	double *phi, double *pi)
{
    /* System generated locals */
    double d__1;

    /* Local variables */
    static int j, k;
    static double r__, s, t, u, v, z__, g0;
    static int km, mcap, mcap1;

/* called as	 fdfilt( x, d, w(ly), slogvk,
 				     w(lamk), w(lak), w(lvk), w(lphi), w(lpi))
     float              x(n)
     double precision  y(n), amk(n), ak(n)
     double precision  vk(M), phi(M), pi(M)
 **************************************************************************
 input  :
          x       float    original time series
          d       double  estimated value of d
 output :
          y       double  filtered series
          slogvk  double  the sum of the logarithms of the vk
 notes  :
          y can use the same storage as either ak or amk
          phi and pi can use the same storage
          can be arranged so that phi, pi and vk share the same storage

 MM:  Which filtering exactly ????
 --   --> look at ./fdsim.f  which is similar (but simpler)

 **************************************************************************
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 -----------------------------------------------------------------------
     Parameter adjustments */
    --pi;
    --phi;
    --vk;
    --ak;
    --amk;
    --y;
    --x;

    /* Function Body */
    mcap = imin2(dimsfd_2.m,dimsfd_2.n);
    mcap1 = mcap + 1;

/* calculate amk(k), vk(k), and ak(k) for k=1,n (see W522-4 for notation). */


/*  k = 1 */

    amk[1] = 0.;
    ak[1] = 1.;

/*  k = 2 ;  initialize phi(1) */

    z__ = *d__ / (1. - *d__);
    amk[2] = z__ * x[1];
    ak[2] = 1. - z__;
    phi[1] = z__;
    d__1 = 1. - *d__;
    t = dgamr_(&d__1);
    if (gammfd_1.igamma != 0) {
	return 0;
    }
    d__1 = 1. - *d__ * 2.;
    g0 = dgamma_(&d__1) * (t * t);
    if (gammfd_1.igamma != 0) {
	return 0;
    }
    vk[1] = g0;
    vk[2] = g0 * (1. - z__ * z__);

/*  k = 3, mcap */

    for (k = 3; k <= mcap; ++k) {
	km = k - 1;
	t = (double) km;
	u = t - *d__;

	/*  calculate phi() and vk() using the recursion formula on W498 */

	for (j = 1; j <= (km - 1); ++j) { /* f2c-clean: s {i__2} {km - 1} */
	    s = t - (double) j;
	    phi[j] *= t * (s - *d__) / (u * s);
	}
	v = *d__ / u;
	phi[km] = v;
	vk[k] = vk[km] * (1. - v * v);

	/*  form amk(k) and ak(k) */

	u = 0.;
	v = 1.;
	for (j = 1; j <= km; ++j) { /* f2c-clean: s {i__2} {km} */
	    t = phi[j];
	    u += t * x[k - j];
	    v -= t;
	}
	amk[k] = u;
	ak[k] = v;
    }

/*     k = mcap+1, n */

    if (dimsfd_2.m < dimsfd_2.n) { /* i.e. mcap = min(M,n) != n */

      /* calculate pi(j), j = 1,mcap */

	pi[1] = *d__;
	s = *d__;
	for (j = 2; j <= mcap; ++j) {
	    u = (double) j;
	    t = pi[j - 1] * ((u - 1. - *d__) / u);
	    s += t;
	    pi[j] = t;
	}
	s = 1. - s;
	r__ = 0.;
	u = (double) mcap;
	t = u * pi[mcap];

	for (k = mcap1; k <= (dimsfd_2.n); ++k) {
	    km = k - mcap;
	    z__ = 0.;
	    for (j = 1; j <= mcap; ++j) {
		z__ += pi[j] * x[k - j];
	    }
	    if (r__ == 0.) {
		amk[k] = z__;
		ak[k] = s;
	    } else {
		d__1 = u / (double) k;
		v = t * (1. - pow(d__1, *d__)) / *d__;
		amk[k] = z__ + v * r__ / ((double) km - 1.);
		ak[k] = s - v;
	    }
	    r__ += x[km];
	}
    }

/*  form muhat - see formula on W523. */

    r__ = 0.;
    s = 0.;
    for (k = 1; k <= (dimsfd_2.n); ++k) {
	t = ak[k];
	u = (x[k] - amk[k]) * t;
	v = t * t;
	if (k <= mcap) {
	    z__ = vk[k];
	    u /= z__;
	    v /= z__;
	}
	r__ += u;
	s += v;
    }
    filtfd_1.hatmu = r__ / s;

/*  form filtered version */

    s = 0.;
    for (k = 1; k <= mcap; ++k)
	s += log(vk[k]);

    *slogvk = s;
    s = 0.;
    for (k = 1; k <= (dimsfd_2.n); ++k) {
	t = x[k] - amk[k] - filtfd_1.hatmu * ak[k];
	if (k <= mcap) {
	    t /= sqrt(vk[k]);
	}
	s += t;
	y[k] = t;
    }
    if (dimsfd_2.npq == 0) {
	return 0;
    }
    t = (double) dimsfd_2.n;
    u = z__ / t;
    for (k = 1; k <= dimsfd_2.n; ++k)
	y[k] -= u;

    return 0;
} /* fdfilt_ */


/****************************************************************************
 ****************************************************************************
 Subroutine */
int ajqp_(double *qp, double *a, double *ajac,
	  int *lajac, int *iflag, double *y)
{
    /* System generated locals */
    int ajac_dim1, ajac_offset;

    /* Local variables */
    static int i__, k, l;
    static double s, t;
    static int km;

/*     double precision qp(npq), a(nm), ajac(nm,npq), y(n)
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 --------------------------------------------------------------------------
     Parameter adjustments */
    --qp;
    --a;
    ajac_dim1 = *lajac;
    ajac_offset = 1 + ajac_dim1;
    ajac -= ajac_offset;
    --y;

    /* Function Body */
    if (*iflag == 2) {
	goto L200;
    }
    if (*iflag != 1) {
	return 0;
    }

/*  objective calculation */

    for (k = dimsfd_2.maxpq1; k <= (dimsfd_2.n); ++k) {
	km = k - dimsfd_2.maxpq;
	t = 0.;
	if (dimsfd_2.np != 0) {
	    for (l = 1; l <= (dimsfd_2.np); ++l) {
		t -= qp[dimsfd_2.nq + l] * y[k - l];
	    }
	}
	s = 0.;
	if (dimsfd_2.nq != 0) {
	    for (l = 1; l <= (dimsfd_2.nq); ++l) {
		if (km <= l) {
		    goto L101;
		}
		s += qp[l] * a[km - l];
	    }
	}
L101:
	s = y[k] + (t + s);
	if (fabs(s) <= mauxfd_1.bignum) {
	    a[km] = s;
	} else {
	    a[km] = sign(s) * mauxfd_1.bignum;
	}
    }
    ++cntrfd_1.nfun;
    return 0;
L200:

/*  jacobian calculation */

    for (i__ = 1; i__ <= (dimsfd_2.npq); ++i__) {
	for (k = dimsfd_2.maxpq1; k <= (dimsfd_2.n); ++k) {
	    km = k - dimsfd_2.maxpq;
	    t = 0.;
	    if (dimsfd_2.nq != 0) {
		for (l = 1; l <= (dimsfd_2.nq); ++l) {
		    if (km <= l) {
			goto L201;
		    }
		    t += qp[l] * ajac[km - l + i__ * ajac_dim1];
		}
	    }
L201:
	    if (i__ <= dimsfd_2.nq) {
		if (km > i__) {
		    s = a[km - i__] + t;
		} else {
		    s = t;
		}
	    } else {
		s = -y[k - (i__ - dimsfd_2.nq)] + t;
	    }
	    if (fabs(s) <= mauxfd_1.bignum) {
		ajac[km + i__ * ajac_dim1] = s;
	    } else {
		ajac[km + i__ * ajac_dim1] = sign(s) * mauxfd_1.bignum;
	    }
	}
    }
    ++cntrfd_1.ngrd;
    return 0;
} /* ajqp_ */

/****************************************************************************
 ****************************************************************************/
/* Subroutine */
int ajp_(double *p, double *a, double *ajac,
	 int *lajac, int *iflag, double *y)
    /*  p(np), a(nm), ajac(nm,npq), y(n) */
{
/* copyright 1991 Department of Statistics, University of Washington
   written by Chris Fraley
   -------------------------------------------------------------------------
   */

    /* System generated locals */
    int ajac_dim1, ajac_offset;

    /* Local variables */
    static int i__, k, l;
    static double t;

    /* Parameter adjustments */
    --p;
    --a;
    ajac_dim1 = *lajac;
    ajac_offset = 1 + ajac_dim1;
    ajac -= ajac_offset;
    --y;

    /* Function Body */
    if (*iflag == 2) {
	goto L200;
    }
    if (*iflag != 1) {
	return 0;
    }
    if (dimsfd_2.np == 0) {
	return 0;
    }

/*  objective calculation */

    for (k = dimsfd_2.np + 1; k <= (dimsfd_2.n); ++k) {
	t = 0.;
	for (l = 1; l <= (dimsfd_2.np); ++l) {
	    t -= p[l] * y[k - l];
	}
	a[k - dimsfd_2.np] = y[k] + t;
    }
    return 0;
L200:

/*  jacobian calculation */

    for (i__ = 1; i__ <= (dimsfd_2.np); ++i__) {
	for (k = dimsfd_2.np + 1; k <= (dimsfd_2.n); ++k) {
	    ajac[k - dimsfd_2.np + i__ * ajac_dim1] = -y[k - i__];
	}
    }
    return 0;
} /* ajp_
 ****************************************************************************
 ****************************************************************************/

/* Subroutine */ int
ajq_(double *qp, double *a, double *ajac,
     int *lajac, int *iflag, double *y)
     /*     double precision qp(npq), a(nm), ajac(nm,npq), y(n) */
{
/* copyright 1991 Department of Statistics, University of Washington
   written by Chris Fraley
   -------------------------------------------------------------------
   */


    /* System generated locals */
    int ajac_dim1, ajac_offset;

    /* Local variables */
    static int i__, k, l;
    static double s, t;
    static int km;

    /* Parameter adjustments */
    --qp;
    --a;
    ajac_dim1 = *lajac;
    ajac_offset = 1 + ajac_dim1;
    ajac -= ajac_offset;
    --y;

    /* Function Body */
    if (*iflag == 2) {
	goto L200;
    }
    if (*iflag != 1) {
	return 0;
    }
    if (dimsfd_2.nq == 0) {
	return 0;
    }

/*  objective calculation */

    for (k = dimsfd_2.maxpq1; k <= (dimsfd_2.n); ++k) {
	km = k - dimsfd_2.maxpq;
	t = 0.;
	if (dimsfd_2.np != 0) {
	    for (l = 1; l <= (dimsfd_2.np); ++l) {
		t -= qp[dimsfd_2.nq + l] * y[k - l];
	    }
	}
	s = 0.;
	if (dimsfd_2.nq != 0) {
	    for (l = 1; l <= (dimsfd_2.nq); ++l) {
		if (km <= l) {
		    goto L101;
		}
		s += qp[l] * a[km - l];
	    }
	}
L101:
	a[km] = y[k] + (t + s);
    }
    ++cntrfd_1.nfun;
    return 0;
L200:

/*  jacobian calculation */

    for (i__ = 1; i__ <= (dimsfd_2.npq); ++i__) {
	for (k = dimsfd_2.maxpq1; k <= (dimsfd_2.n); ++k) {
	    km = k - dimsfd_2.maxpq;
	    t = 0.;
	    if (dimsfd_2.nq != 0) {
		for (l = 1; l <= (dimsfd_2.nq); ++l) {
		    if (km <= l) {
			goto L201;
		    }
		    t += qp[l] * ajac[km - l + i__ * ajac_dim1];
		}
	    }
L201:
	    if (i__ <= dimsfd_2.nq) {
		if (km > i__) {
		    ajac[km + i__ * ajac_dim1] = a[km - i__] + t;
		} else {
		    ajac[km + i__ * ajac_dim1] = t;
		}
	    } else {
		ajac[km + i__ * ajac_dim1] = -y[k - (i__ - dimsfd_2.nq)] + t;
	    }
	}
    }
    ++cntrfd_1.ngrd;
    return 0;
} /* ajq_ */

