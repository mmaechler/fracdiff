/* fdcore.f -- translated by f2c (version 20031025).
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
    doublereal epsp25, epspt3, epspt5, epsp75, bignum;
} mauxfd_;

#define mauxfd_1 mauxfd_

union {
    struct {
	integer nn, mm, np, nq, npq, npq1, maxpq, maxpq1, minpq, nm;
    } _1;
    struct {
	integer n, m, np, nq, npq, npq1, maxpq, maxpq1, minpq, nm;
    } _2;
} dimsfd_;

#define dimsfd_1 (dimsfd_._1)
#define dimsfd_2 (dimsfd_._2)

struct {
    integer maxopt, maxfun, nopt, nfun, ngrd, ifun, igrd, info;
} cntrfd_;

#define cntrfd_1 cntrfd_

union {
    struct {
	doublereal told, tolf, tolx, tolg, anorm, deltax, gnorm;
    } _1;
    struct {
	doublereal dtol, ftol, xtol, gtol, anorm, deltax, gnorm;
    } _2;
} tolsfd_;

#define tolsfd_1 (tolsfd_._1)
#define tolsfd_2 (tolsfd_._2)

struct {
    integer ly, lamk, lak, lvk, lphi, lpi;
} wfilfd_;

#define wfilfd_1 wfilfd_

struct {
    integer lqp, la, lajac, ipvt, ldiag, lqtf, lwa1, lwa2, lwa3, lwa4;
} woptfd_;

#define woptfd_1 woptfd_

struct {
    integer ilimit, jlimit;
} limsfd_;

#define limsfd_1 limsfd_

struct {
    integer igamma, jgamma;
} gammfd_;

#define gammfd_1 gammfd_

struct {
    integer iminpk, jminpk;
} mnpkfd_;

#define mnpkfd_1 mnpkfd_

struct {
    integer ksvd, kcov, kcor;
} hessfd_;

#define hessfd_1 hessfd_

union {
    struct {
	doublereal hatmu, wnv, cllf;
    } _1;
    struct {
	doublereal hatmu, wnv, hood;
    } _2;
} filtfd_;

#define filtfd_1 (filtfd_._1)
#define filtfd_2 (filtfd_._2)

/* Table of constant values */

static real c_b2 = -99.f;
static integer c__1 = 1;
static integer c__0 = 0;
static doublereal c_b21 = 1.;

/* ****************************************************************************** */
/* ****************************************************************************** */
/* Subroutine */ int fracdf_(doublereal *x, integer *n, integer *m, integer *
	nar, integer *nma, doublereal *dtol, doublereal *drange, doublereal *
	hood, doublereal *d__, doublereal *ar, doublereal *ma, doublereal *w, 
	integer *lenw, integer *inform__, doublereal *flmin, doublereal *
	flmax, doublereal *epmin, doublereal *epmax)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    extern doublereal dopt_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern /* Subroutine */ int fdcom_(integer *, integer *, integer *, 
	    integer *, real *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal delta;
    static integer lfree;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer lwfree, lenthw;

/*     real               x(n) */
/*     double precision   ar(*), ma(*), drange(2) */
/*     double precision   w(*) */
/* ------------------------------------------------------------------------------ */

/*   Input : */

/*  x       double   time series for the ARIMA model */
/*  n       integer  length of the time series */
/*  M       integer  number of terms in the likelihood approximation */
/*                   suggested value 100 (see Haslett and Raftery 1989) */
/*  nar     integer  number of autoregressive parameters */
/*  nma     integer  number of moving average parameters */
/*  dtol    double   desired length of final interval of uncertainty for d */
/*                   suggested value : 4th root of machine precision */
/*                   if dtol < 0 it is automatically set to this value */
/*                   dtol will be altered if necessary by the program */
/*  drange  double   array of length 2 giving minimum and maximum values f */
/*                   for the fractional differencing parameter */
/*  d       double   initial guess for optimal fractional differencing parameter */
/*  w       double   work array */
/*  lenw    integer  length of double precision workspace w, must be at least */
/* 		max( p+q+2*(n+M), 3*n+(n+6.5)*(p+q) +1, (3+2*(p+q+1))*(p+q+1)+1) */
/*   MM:		max( p+q+2*(n+M), 3*n+(n+6.5)*(p+q) +1,     31 * 12) */
/* 	 is what the code below rather checks */

/*  Output : */

/*  dtol    double   value of dtol ultimately used by the algorithm */
/*  d       double   final value optimal fractional differencing parameter */
/*  hood    double   logarithm of the maximum likelihood */
/*  ar      double   optimal autoregressive parameters */
/*  ma      double   optimal moving average parameters */

/* ------------------------------------------------------------------------------ */
/*  copyright 1991 Department of Statistics, University of Washington */
/*  written by Chris Fraley */
/* ----------------------------------------------------------------------------- */
    /* Parameter adjustments */
    --x;
    --ar;
    --ma;
    --drange;
    --w;

    /* Function Body */
    if (*m <= 0) {
	*m = 100;
    }
/*     MM: Using 'fdcom' instead of 'code copy' -- FIXME: use #include in C */
/*     initialize several of the above common blocks: */
    fdcom_(n, m, nar, nma, &c_b2, flmin, flmax, epmin, epmax);
    lfree = woptfd_1.lwa4 + *n - dimsfd_1.minpq;
/* 	= 1+ ipvt + 5.5*npq + n - minpq */
/* 	= 2+ 6.5*npq + 3*n - 2*minpq + (n-maxpq)*npq */
/* and               lvk+M = 1 + npq + 2(n + M) */
/* Computing MAX */
    i__1 = wfilfd_1.lvk + *m, i__1 = max(i__1,lfree);
    lwfree = max(i__1,372);
/*                                 ^^^^^^^ MM: where is this needed? */
    if (lwfree > *lenw + 1) {
	limsfd_1.ilimit = lwfree - *lenw;
/*       write( 6, *) 'insufficient storage : ', */
/*    *               'increase length of w by at least', ILIMIT */
	*inform__ = 1;
/* 	return the *desired* workspace storage: */
	*lenw = lwfree;
	return 0;
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
	tolsfd_1.told = max(*dtol,mauxfd_1.epspt5);
/* Computing MAX */
	d__1 = *dtol / 10.;
	tolsfd_1.tolf = max(d__1,mauxfd_1.epsp75);
    }
    tolsfd_1.tolg = tolsfd_1.tolf;
    tolsfd_1.tolx = tolsfd_1.told;
    *dtol = tolsfd_1.told;
/*     if (npq .ne. 0) call dcopy( npq, zero, 0, w(lqp), 1) */
    if (dimsfd_1.npq != 0) {
	dcopy_(&dimsfd_1.np, &ar[1], &c__1, &w[woptfd_1.lqp + dimsfd_1.nq], &
		c__1);
	dcopy_(&dimsfd_1.nq, &ma[1], &c__1, &w[woptfd_1.lqp], &c__1);
    }
    cntrfd_1.nopt = 0;
    cntrfd_1.nfun = 0;
    cntrfd_1.ngrd = 0;
/* 	  ==== */
    *d__ = dopt_(&x[1], d__, &drange[1], hood, &delta, &w[1]);
/* 	  ==== */
    if (cntrfd_1.nopt >= cntrfd_1.maxopt) {
	limsfd_1.jlimit = 1;
    }
/*       write( 6, *) */
/*       write( 6, *) 'WARNING : optimization limit reached' */
/*     end if */
    if (gammfd_1.igamma != 0 || mnpkfd_1.iminpk != 0) {
	*d__ = machfd_1.fltmax;
	*hood = machfd_1.fltmax;
	dcopy_(&dimsfd_1.np, &machfd_1.fltmax, &c__0, &ar[1], &c__1);
	dcopy_(&dimsfd_1.nq, &machfd_1.fltmax, &c__0, &ma[1], &c__1);
	if (gammfd_1.igamma != 0) {
	    *inform__ = 2;
	}
	if (mnpkfd_1.iminpk != 0) {
	    *inform__ = 3;
	}
	return 0;
    }
    dcopy_(&dimsfd_1.np, &w[woptfd_1.lqp + dimsfd_1.nq], &c__1, &ar[1], &c__1)
	    ;
    dcopy_(&dimsfd_1.nq, &w[woptfd_1.lqp], &c__1, &ma[1], &c__1);
    if (gammfd_1.jgamma != 0) {
	*inform__ = 4;
    }
    if (mnpkfd_1.jminpk != 0) {
	*inform__ = 5;
    }
    if (limsfd_1.jlimit != 0) {
	*inform__ = 6;
    }
    return 0;
/* 900  format( 4h itr, 14h     d          ,   14h    est mean  , */
/*     *                16h     white noise,  17h     log likelihd, */
/*     *                 4h  nf, 3h ng) */
} /* fracdf_ */

/*     fracdf() {main} */
/* ****************************************************************************** */
/* ****************************************************************************** */

/* optimization with respect to d based on Brent's fmin algorithm */

doublereal dopt_(doublereal *x, doublereal *dinit, doublereal *drange, 
	doublereal *hood, doublereal *delta, doublereal *w)
{
    /* Initialized data */

    static doublereal cc = .38196601125011;

    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal d__, aa, bb, fa, dd, fb, ee, hh, fu, fv, fw, fx, rr, ss,
	     tt, uu, vv, ww, xx, eps, tol, tol1, tol2, tol3;
    extern doublereal pqopt_(doublereal *, doublereal *, doublereal *);

/*     real              x(n) */
/*  copyright 1991 Department of Statistics, University of Washington */
/*  written by Chris Fraley */
/* ------------------------------------------------------------------------------ */

/*  cc is the squared inverse of the golden ratio (see data statement) */

/*     cc = half*(three-sqrt(5.0d0)) */
    /* Parameter adjustments */
    --w;
    --drange;
    --x;

    /* Function Body */

/*  eps is approximately the square root of the relative machine */
/*  precision. */

    eps = machfd_1.epsmax;
    tol1 = eps + 1.;
    eps = sqrt(eps);
/* -Wall: */
    ret_val = -1.;
    dd = 0.;

    aa = drange[1];
    bb = drange[2];
    if (*dinit > aa + tolsfd_2.dtol && *dinit < bb - tolsfd_2.dtol) {
	vv = *dinit;
    } else {
	vv = aa + cc * (bb - aa);
    }
    ww = vv;
    xx = vv;
    uu = xx;
    ee = 0.;
    cntrfd_1.nopt = 1;
    fx = pqopt_(&x[1], &xx, &w[1]);
/*             ===== */
    fv = fx;
    fw = fx;
    tol = max(tolsfd_2.dtol,0.);
    tol3 = tol / 3.;

/*  main loop starts here */

L10:
    if (gammfd_1.igamma != 0 || mnpkfd_1.iminpk != 0) {
	d__ = uu;
	*hood = machfd_1.fltmax;
	return ret_val;
    }
    hh = (aa + bb) * .5;
    tol1 = eps * (abs(xx) + 1.) + tol3;
    tol2 = tol1 * 2.;

/*  check stopping criterion */

    *delta = (d__1 = xx - hh, abs(d__1)) + (bb - aa) * .5;
/*     if (abs(xx-hh) .le. (tol2-half*(bb-aa))) goto 100 */
    if (*delta <= tol2) {
	goto L100;
    }
    if (cntrfd_1.nopt >= cntrfd_1.maxopt) {
	goto L100;
    }
/*     if (delpq .le. EPSMAX*(one+pqnorm)) goto 100 */
    rr = 0.;
    ss = 0.;
    tt = 0.;
    if (abs(ee) > tol1) {

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
    if (abs(tt) >= (d__1 = ss * .5 * rr, abs(d__1)) || tt <= ss * (aa - xx) ||
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

    if (abs(dd) >= tol1) {
	uu = xx + dd;
    } else {
	if (dd <= 0.) {
	    uu = xx - tol1;
	} else {
	    uu = xx + tol1;
	}
    }
    ++cntrfd_1.nopt;
    fu = pqopt_(&x[1], &uu, &w[1]);

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
/* 900  format( i4, 2(1pe14.6), 1pe16.7, 1pe17.8, 1x, 2(i3)) */
/* 901  format( i4, 3(1pe10.2), 1pe11.2, 2(i3), 3(1pe8.1), i2) */
} /* dopt_ */

/*     dopt() */
/* ************************************************************************** */
/* ************************************************************************** */
doublereal pqopt_(doublereal *x, doublereal *d__, doublereal *w)
{
    /* Initialized data */

    static integer modelm = 1;
    static doublereal factlm = 100.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static doublereal t, u, bic;
    extern /* Subroutine */ int ajp_(), ajq_();
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int ajqp_();
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), lmder1_(U_fp, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), fdfilt_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal slogvk;

/*     real              x(n) */
/*     work array exactly as in main  fracdf() : */
/*     These are passed to LMDER1() to be optimized: */
    /* Parameter adjustments */
    --w;
    --x;

    /* Function Body */
/*     copyright 1991 Department of Statistics, University of Washington */
/*     written by Chris Fraley */
/* ---------------------------------------------------------------------------- */
    fdfilt_(&x[1], d__, &w[(0 + (0 + (wfilfd_1.ly << 3))) / 8], &slogvk, &w[(
	    0 + (0 + (wfilfd_1.lamk << 3))) / 8], &w[(0 + (0 + (wfilfd_1.lak 
	    << 3))) / 8], &w[(0 + (0 + (wfilfd_1.lvk << 3))) / 8], &w[(0 + (0 
	    + (wfilfd_1.lphi << 3))) / 8], &w[(0 + (0 + (wfilfd_1.lpi << 3))) 
	    / 8]);
    if (gammfd_1.igamma != 0) {
	ret_val = machfd_1.fltmax;
	filtfd_2.wnv = machfd_1.fltmax;
	filtfd_2.hood = -machfd_1.fltmax;
	return ret_val;
    }
    t = (doublereal) dimsfd_2.n;
    if (dimsfd_2.npq == 0) {
/* 	trivial case  p = q = 0 : */
	filtfd_2.wnv = ddot_(&dimsfd_2.n, &w[wfilfd_1.ly], &c__1, &w[
		wfilfd_1.ly], &c__1) / t;
	cntrfd_1.ifun = 0;
	cntrfd_1.igrd = 0;
	cntrfd_1.info = -1;
    } else {

/*     optimize as an unconstrained optimization problem */

	if (modelm == 2) {
	    dcopy_(&dimsfd_2.npq, &c_b21, &c__0, &w[woptfd_1.ldiag], &c__1);
	}
	if (cntrfd_1.nopt < 0) {
	    if (dimsfd_2.np != 0) {
		i__1 = dimsfd_2.n - dimsfd_2.np;
		i__2 = dimsfd_2.n - dimsfd_2.np;
		lmder1_((U_fp)ajp_, &i__1, &dimsfd_2.np, &w[woptfd_1.lqp + 
			dimsfd_2.nq], &w[woptfd_1.la], &w[woptfd_1.lajac], &
			i__2, &tolsfd_2.ftol, &tolsfd_2.xtol, &tolsfd_2.gtol, 
			&cntrfd_1.maxfun, &w[woptfd_1.ldiag], &modelm, &
			factlm, &cntrfd_1.info, &cntrfd_1.ifun, &
			cntrfd_1.igrd, &w[woptfd_1.ipvt], &w[woptfd_1.lqtf], &
			w[woptfd_1.lwa1], &w[woptfd_1.lwa2], &w[woptfd_1.lwa3]
			, &w[woptfd_1.lwa4], &w[wfilfd_1.ly]);
	    }
	    if (dimsfd_2.nq != 0) {
		i__1 = dimsfd_2.n - dimsfd_2.nq;
		i__2 = dimsfd_2.n - dimsfd_2.nq;
		lmder1_((U_fp)ajq_, &i__1, &dimsfd_2.nq, &w[woptfd_1.lqp], &w[
			woptfd_1.la], &w[woptfd_1.lajac], &i__2, &
			tolsfd_2.ftol, &tolsfd_2.xtol, &tolsfd_2.gtol, &
			cntrfd_1.maxfun, &w[woptfd_1.ldiag], &modelm, &factlm,
			 &cntrfd_1.info, &cntrfd_1.ifun, &cntrfd_1.igrd, &w[
			woptfd_1.ipvt], &w[woptfd_1.lqtf], &w[woptfd_1.lwa1], 
			&w[woptfd_1.lwa2], &w[woptfd_1.lwa3], &w[
			woptfd_1.lwa4], &w[wfilfd_1.ly]);
	    }
	}
	lmder1_((U_fp)ajqp_, &dimsfd_2.nm, &dimsfd_2.npq, &w[woptfd_1.lqp], &
		w[woptfd_1.la], &w[woptfd_1.lajac], &dimsfd_2.nm, &
		tolsfd_2.ftol, &tolsfd_2.xtol, &tolsfd_2.gtol, &
		cntrfd_1.maxfun, &w[woptfd_1.ldiag], &modelm, &factlm, &
		cntrfd_1.info, &cntrfd_1.ifun, &cntrfd_1.igrd, &w[
		woptfd_1.ipvt], &w[woptfd_1.lqtf], &w[woptfd_1.lwa1], &w[
		woptfd_1.lwa2], &w[woptfd_1.lwa3], &w[woptfd_1.lwa4], &w[
		wfilfd_1.ly]);
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
/*     call daxpy( npq, (-one), w(lpq), 1, w(lqp), 1 */
/*     delpq  = sqrt(ddot( npq, w(lqp), 1, w(lqp), 1)) */
/*     pqnorm = sqrt(ddot( npq, w(lpq), 1, w(lpq), 1)) */
	filtfd_2.wnv = tolsfd_2.anorm * tolsfd_2.anorm / (doublereal) (
		dimsfd_2.nm - 1);
    }
    u = t * (log(filtfd_2.wnv) + 2.8378) + slogvk;
    ret_val = u / 2.;
    bic = u + (doublereal) (dimsfd_2.np + dimsfd_2.nq + 1) * log(t);
    filtfd_2.hood = -ret_val;
    return ret_val;
} /* pqopt_ */

/*     pqopt() */
/* ************************************************************************** */
/* ************************************************************************** */
/* Subroutine */ int fdfilt_(doublereal *x, doublereal *d__, doublereal *y, 
	doublereal *slogvk, doublereal *amk, doublereal *ak, doublereal *vk, 
	doublereal *phi, doublereal *pi)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal), sqrt(
	    doublereal);

    /* Local variables */
    static integer j, k;
    static doublereal r__, s, t, u, v, z__, g0;
    static integer km, mcap, mcap1;
    extern doublereal dgamr_(doublereal *), dgamma_(doublereal *);

/* called as	 fdfilt( x, d, w(ly), slogvk, */
/* 				     w(lamk), w(lak), w(lvk), w(lphi), w(lpi)) */
/*     real              x(n) */
/*     double precision  y(n), amk(n), ak(n) */
/*     double precision  vk(M), phi(M), pi(M) */
/* ************************************************************************** */
/* input  : */
/*          x       real    original time series */
/*          d       double  estimated value of d */
/* output : */
/*          y       double  filtered series */
/*          slogvk  double  the sum of the logarithms of the vk */
/* notes  : */
/*          y can use the same storage as either ak or amk */
/*          phi and pi can use the same storage */
/*          can be arranged so that phi, pi and vk share the same storage */

/* MM:  Which filtering exactly ???? */
/* --   --> look at ./fdsim.f  which is similar (but simpler) */

/* ************************************************************************** */
/* copyright 1991 Department of Statistics, University of Washington */
/* written by Chris Fraley */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --pi;
    --phi;
    --vk;
    --ak;
    --amk;
    --y;
    --x;

    /* Function Body */
    mcap = min(dimsfd_2.m,dimsfd_2.n);
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

    i__1 = mcap;
    for (k = 3; k <= i__1; ++k) {
	km = k - 1;
	t = (doublereal) km;
	u = t - *d__;

/*  calculate phi() and vk() using the recursion formula on W498 */

	i__2 = km - 1;
	for (j = 1; j <= i__2; ++j) {
	    s = t - (doublereal) j;
	    phi[j] *= t * (s - *d__) / (u * s);
	}
	v = *d__ / u;
	phi[km] = v;
	vk[k] = vk[km] * (1. - v * v);

/*  form amk(k) and ak(k) */

	u = 0.;
	v = 1.;
	i__2 = km;
	for (j = 1; j <= i__2; ++j) {
	    t = phi[j];
	    u += t * x[k - j];
	    v -= t;
	}
	amk[k] = u;
	ak[k] = v;
    }
    if (dimsfd_2.m < dimsfd_2.n) {
/* 	   i.e. mcap = min(M,n) != n */

/*     k = mcap+1, n */

/*     calculate pi(j), j = 1,mcap */

	pi[1] = *d__;
	s = *d__;
	i__1 = mcap;
	for (j = 2; j <= i__1; ++j) {
	    u = (doublereal) j;
	    t = pi[j - 1] * ((u - 1. - *d__) / u);
	    s += t;
	    pi[j] = t;
	}
	s = 1. - s;
	r__ = 0.;
	u = (doublereal) mcap;
	t = u * pi[mcap];

	i__1 = dimsfd_2.n;
	for (k = mcap1; k <= i__1; ++k) {
	    km = k - mcap;
	    z__ = 0.;
	    i__2 = mcap;
	    for (j = 1; j <= i__2; ++j) {
		z__ += pi[j] * x[k - j];
	    }
	    if (r__ == 0.) {
		amk[k] = z__;
		ak[k] = s;
	    } else {
		d__1 = u / (doublereal) k;
		v = t * (1. - pow_dd(&d__1, d__)) / *d__;
		amk[k] = z__ + v * r__ / ((doublereal) km - 1.);
		ak[k] = s - v;
	    }
	    r__ += x[km];
	}
    }

/*  form muhat - see formula on W523. */

    r__ = 0.;
    s = 0.;
    i__1 = dimsfd_2.n;
    for (k = 1; k <= i__1; ++k) {
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
    i__1 = mcap;
    for (k = 1; k <= i__1; ++k) {
	s += log(vk[k]);
    }
    *slogvk = s;
    s = 0.;
    i__1 = dimsfd_2.n;
    for (k = 1; k <= i__1; ++k) {
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
    t = (doublereal) dimsfd_2.n;
    u = z__ / t;
    i__1 = dimsfd_2.n;
    for (k = 1; k <= i__1; ++k) {
	y[k] -= u;
    }
    return 0;
} /* fdfilt_ */

/*     fdfilt() */
/* **************************************************************************** */
/* **************************************************************************** */
/* Subroutine */ int ajqp_(doublereal *qp, doublereal *a, doublereal *ajac, 
	integer *lajac, integer *iflag, doublereal *y)
{
    /* System generated locals */
    integer ajac_dim1, ajac_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, k, l;
    static doublereal s, t;
    static integer km;

/*     double precision qp(npq), a(nm), ajac(nm,npq), y(n) */
/* copyright 1991 Department of Statistics, University of Washington */
/* written by Chris Fraley */
/* -------------------------------------------------------------------------- */
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

/*  objective calculation */

    i__1 = dimsfd_2.n;
    for (k = dimsfd_2.maxpq1; k <= i__1; ++k) {
	km = k - dimsfd_2.maxpq;
	t = 0.;
	if (dimsfd_2.np != 0) {
	    i__2 = dimsfd_2.np;
	    for (l = 1; l <= i__2; ++l) {
		t -= qp[dimsfd_2.nq + l] * y[k - l];
	    }
	}
	s = 0.;
	if (dimsfd_2.nq != 0) {
	    i__2 = dimsfd_2.nq;
	    for (l = 1; l <= i__2; ++l) {
		if (km <= l) {
		    goto L101;
		}
		s += qp[l] * a[km - l];
	    }
	}
L101:
	s = y[k] + (t + s);
	if (abs(s) <= mauxfd_1.bignum) {
	    a[km] = s;
	} else {
	    a[km] = d_sign(&c_b21, &s) * mauxfd_1.bignum;
	}
    }
    ++cntrfd_1.nfun;
    return 0;
L200:

/*  jacobian calculation */

    i__1 = dimsfd_2.npq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = dimsfd_2.n;
	for (k = dimsfd_2.maxpq1; k <= i__2; ++k) {
	    km = k - dimsfd_2.maxpq;
	    t = 0.;
	    if (dimsfd_2.nq != 0) {
		i__3 = dimsfd_2.nq;
		for (l = 1; l <= i__3; ++l) {
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
	    if (abs(s) <= mauxfd_1.bignum) {
		ajac[km + i__ * ajac_dim1] = s;
	    } else {
		ajac[km + i__ * ajac_dim1] = d_sign(&c_b21, &s) * 
			mauxfd_1.bignum;
	    }
	}
    }
    ++cntrfd_1.ngrd;
    return 0;
} /* ajqp_ */

/* **************************************************************************** */
/* **************************************************************************** */
/* Subroutine */ int ajp_(doublereal *p, doublereal *a, doublereal *ajac, 
	integer *lajac, integer *iflag, doublereal *y)
{
    /* System generated locals */
    integer ajac_dim1, ajac_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, l;
    static doublereal t;

/*     double precision p(np), a(nm), ajac(nm,npq), y(n) */
/* copyright 1991 Department of Statistics, University of Washington */
/* written by Chris Fraley */
/* -------------------------------------------------------------------------- */
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

    i__1 = dimsfd_2.n;
    for (k = dimsfd_2.np + 1; k <= i__1; ++k) {
	t = 0.;
	i__2 = dimsfd_2.np;
	for (l = 1; l <= i__2; ++l) {
	    t -= p[l] * y[k - l];
	}
/* L101: */
	a[k - dimsfd_2.np] = y[k] + t;
    }
    return 0;
L200:

/*  jacobian calculation */

    i__1 = dimsfd_2.np;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = dimsfd_2.n;
	for (k = dimsfd_2.np + 1; k <= i__2; ++k) {
	    ajac[k - dimsfd_2.np + i__ * ajac_dim1] = -y[k - i__];
	}
    }
    return 0;
} /* ajp_ */

/* **************************************************************************** */
/* **************************************************************************** */
/* Subroutine */ int ajq_(doublereal *qp, doublereal *a, doublereal *ajac, 
	integer *lajac, integer *iflag, doublereal *y)
{
    /* System generated locals */
    integer ajac_dim1, ajac_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, k, l;
    static doublereal s, t;
    static integer km;

/*     double precision qp(npq), a(nm), ajac(nm,npq), y(n) */
/* copyright 1991 Department of Statistics, University of Washington */
/* written by Chris Fraley */
/* -------------------------------------------------------------------------- */
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

    i__1 = dimsfd_2.n;
    for (k = dimsfd_2.maxpq1; k <= i__1; ++k) {
	km = k - dimsfd_2.maxpq;
	t = 0.;
	if (dimsfd_2.np != 0) {
	    i__2 = dimsfd_2.np;
	    for (l = 1; l <= i__2; ++l) {
		t -= qp[dimsfd_2.nq + l] * y[k - l];
	    }
	}
	s = 0.;
	if (dimsfd_2.nq != 0) {
	    i__2 = dimsfd_2.nq;
	    for (l = 1; l <= i__2; ++l) {
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

    i__1 = dimsfd_2.npq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = dimsfd_2.n;
	for (k = dimsfd_2.maxpq1; k <= i__2; ++k) {
	    km = k - dimsfd_2.maxpq;
	    t = 0.;
	    if (dimsfd_2.nq != 0) {
		i__3 = dimsfd_2.nq;
		for (l = 1; l <= i__3; ++l) {
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

