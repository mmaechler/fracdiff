/*-*- mode: C; kept-old-versions: 12;  kept-new-versions: 20; -*-
 *
 * fdcore.f -- translated by f2c (version 20031025).
 * and produced by  f2c-clean,v 1.10 2002/03/28 16:37:27 maechler
 *
 * and manually pretty edited by Martin Maechler, 2004-09-18
 */

#include <Rmath.h>
// for warning(), and "monitoring" output
#include <R.h>


/* dcopy() and ddot() only:*/
#include <R_ext/BLAS.h>

#include "fracdiff.h"

extern double dgamr_(double *);
extern double dgamma_(double *);

static
double dopt(double *x, double dinit, double *drange, int verbose,
	    double *hood, double *delta, double *w, int *iw, double *min_fnorm);

static
double pqopt(double *x, double d__, double *w, int *iw, double *min_fnorm);


/* These + ajqp_(..)  are passed to LMDER1() to be optimized: */
typedef void /* Unknown procedure type */ (*U_fp)();

static void
ajp_(double *p, double *a, double *ajac, int lajac, int op_code, double *y);

static void
ajq_(double *qp, double *a, double *ajac, int lajac, int op_code, double *y);

/* Common Block Declarations */

/* 1 - local ones  --- MM: maybe get rid of (some of) them : */
static
struct { int maxopt, maxfun, nopt, nfun, ngrd, ifun, igrd, info; } OP;

static struct { double d, f, x, g; } TOL;

static struct { int iminpk, jminpk; } MinPck;

static struct { int ilimit, jlimit; } limsfd_;

/* 2 - global ones ---
 * all defined here :*/
#define FD_EXTERNAL

#include "mach_comm.h"
#include "maux_comm.h"

#include "gamm_comm.h"

#include "hess_comm.h"


/* Table of constant values (used as pointers) */

static double c_m99 = -99.;
static int ic__1 = 1;
static int ic__0 = 0;
static double c__1 = 1.;

/*****************************************************************************
 ******************************************************************************/
void fracdf(double *x, int *n, int *m, int *nar, int *nma,
	    double *dtol, double *drange, double *hood_etc,
	    double *d__, double *ar, double *ma,
	    double *w, int *lenw, int *iw,
	    int *inform, // <- also use as input for verbose
	    double *flmin, double *flmax, double *epmin, double *epmax)
{
/* ----------------------------------------------------------------------------
   Input :

  x(n)    double   time series for the ARIMA model
  n       int  length of the time series
  m       int  number of terms in the likelihood approximation
                   suggested value 100 (see Haslett and Raftery 1989)
  nar     int  number of autoregressive parameters
  nma     int  number of moving average parameters
  dtol    double   desired length of final interval of uncertainty for d
                   suggested value : 4th root of machine precision
                   if dtol < 0 it is automatically set to this value
                   dtol will be altered if necessary by the program
  drange(2) double array of length 2 giving minimum and maximum values f
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
  hood_etc double[3] [1]: logarithm of the maximum likelihood
                     [2]: minimal objective value
		     [3]: estimated noise variance
  ar      double   optimal autoregressive parameters
  ma      double   optimal moving average parameters

  ----------------------------------------------------------------------------
  copyright 1991 Department of Statistics, University of Washington
  written by Chris Fraley
  ----------------------------------------------------------------------------*/

#define MIN_FNORM &hood_etc[1]
#define NS_VAR     hood_etc[2]

    /* Local variables */
    double delta;
    int lfree, lwfree,
	verbose = inform[0],
	w_lqp = w_opt.lqp - 1;// '-1' : so we do *not* need 'w--'

    if (*m <= 0) /* default: */
	*m = 100;

/* MM: Using 'fdcom' instead of 'code copy' -- FIXME: use #include in C
 *     initialize several of the above common blocks: */
    fdcom(n, m, nar, nma, &c_m99, flmin, flmax, epmin, epmax);

    lfree = w_opt.lwa4 + *n - Dims.minpq;
/*	= 1+ ipvt + 5.5*npq + n - minpq
	= 2+ 6.5*npq + 3*n - 2*minpq + (n-maxpq)*npq
 and               lvk+M = 1 + npq + 2(n + M)
 */

    lwfree = imax2((12*31), imax2(w_fil.lvk + *m, lfree));
/*                 ^^^^^^^ MM: where is this needed? */
    if (lwfree > *lenw + 1) {
	limsfd_.ilimit = lwfree - *lenw;
	REprintf("** Insufficient storage : Increase length of w by at least %d\n",
		 limsfd_.ilimit);
	*inform = 1;
	/* return the *desired* workspace storage: */
	*lenw = lwfree;
	return;
    }
    OP.maxopt = 100;
    OP.maxfun = 100;
/* set error and warning flags */
    *inform = 0;
    gammfd_.igamma = 0;
    MinPck.iminpk = 0;
    limsfd_.ilimit = 0;
    gammfd_.jgamma = 0;
    MinPck.jminpk = 0;
    limsfd_.jlimit = 0;

    if (*dtol > .1)
	*dtol = .1;
    if (*dtol <= 0.) {
	TOL.d = mauxfd_.epsp25;
	TOL.f = mauxfd_.epspt3;
    } else {
	TOL.d = fmax2(*dtol, mauxfd_.epspt5);
	TOL.f = fmax2(*dtol / 10., mauxfd_.epsp75);
    }
    TOL.g = TOL.f;
    TOL.x = TOL.d;
    *dtol = TOL.d;
/*     if (npq != 0) call dcopy( npq, zero, 0, w(lqp), 1) */
    if (Dims.pq != 0) {
	F77_CALL(dcopy)(&Dims.p, ar, &ic__1, &w[w_lqp + Dims.q], &ic__1);
	F77_CALL(dcopy)(&Dims.q, ma, &ic__1, &w[w_lqp], &ic__1);
    }
    OP.nopt = 0;
    OP.nfun = 0;
    OP.ngrd = 0;
/* 	   ==== */
    *d__ = dopt(x, *d__, drange, verbose,
/* 	   ===*/
		&hood_etc[0], &delta, w, iw, /* min_fnorm = */&hood_etc[1]);

    hood_etc[2] = filtfd_.wnv;
    if (OP.nopt >= OP.maxopt) {
	limsfd_.jlimit = 1;
	warning("fracdf(): optimization iteration limit %d reached", OP.maxopt);
    }

    if (gammfd_.igamma != 0 || MinPck.iminpk != 0) {
	*d__ = machfd_.fltmax;
	hood_etc[0] = machfd_.fltmax;
	F77_CALL(dcopy)(&Dims.p, &machfd_.fltmax, &ic__0, ar, &ic__1);
	F77_CALL(dcopy)(&Dims.q, &machfd_.fltmax, &ic__0, ma, &ic__1);

	if (gammfd_.igamma != 0) { *inform = 2; return; }
	if (MinPck.iminpk != 0) { *inform = 3; return; }
    }
    F77_CALL(dcopy)(&Dims.p, &w[w_lqp + Dims.q], &ic__1, ar, &ic__1);
    F77_CALL(dcopy)(&Dims.q, &w[w_lqp],          &ic__1, ma, &ic__1);

    if (gammfd_.jgamma != 0) { *inform = 4; return; }
    if (MinPck.jminpk != 0) { *inform = 5; return; }
    if (limsfd_.jlimit != 0) { *inform = 6; }
    return;
/* 900  format( 4h itr, 14h     d          ,   14h    est mean  ,
     *                16h     white noise,  17h     log likelihd,
     *                 4h  nf, 3h ng) */

} /* fracdf() {main} */


/******************************************************************************
 *****************************************************************************

 optimization with respect to d based on Brent's fmin algorithm */

static
double dopt(double *x, double dinit, double *drange, int verbose,
	    double *hood, double *delta, double *w, int *iw, double *min_fnorm)

{
/*     float              x(n) */

    /* cc is the squared inverse of the golden ratio, cc := (3-sqrt(5.))/2 : */
    static double cc = .38196601125011;

    static double aa, bb, dd, ee, hh, fu, fv, fw, fx, rr, ss,
	     tt, uu, vv, ww, xx, eps, tol, tol1, tol2, tol3;

/*  copyright 1991 Department of Statistics, University of Washington
  written by Chris Fraley
 ------------------------------------------------------------------------------
*/

    /* eps is approximately the square root of the relative machine precision. */
    eps = machfd_.epsmax;
    tol1 = eps + 1.;
    eps = sqrt(eps);

    aa = drange[0];
    bb = drange[1];
    if (dinit > aa + TOL.d &&
	dinit < bb - TOL.d) {
	vv = dinit;
    } else {
	vv = aa + cc * (bb - aa);
    }
    ww = vv;
    xx = vv;
    uu = xx;
    dd = 0.;
    ee = 0.;
    OP.nopt = 1;
    fx = pqopt(x, xx, w, iw, min_fnorm);
/*       ===== */
    if(verbose) {
	REprintf("dopt() debugging: dinit = %g ==> xx = %g, fx = pqopt(x[], xx) = %g; min_fnorm = %g\n",
		 dinit, xx, fx, *min_fnorm);
	REprintf("                 it. |           uu |    pqopt(uu) |      delta |\n");
	//                             | 123456789012 | 123456789012 | 1234567890 |\n",
//      REprintf(" .. DBG dopt() [%2d]:| %12g | %12g | %10.6e |\n",
    }
    fv = fx;
    fw = fx;
    tol = fmax2(TOL.d, 0.);
    tol3 = tol / 3.;

/* main loop starts here ======================================================*/
L10:
    if (gammfd_.igamma != 0 || MinPck.iminpk != 0) {
	*hood = machfd_.fltmax;
	warning("** dopt() ERROR:  invalid gamma (%d) or Minpack (%d) codes",
		gammfd_.igamma, MinPck.iminpk);
	return -1.;
    }
    hh = (aa + bb) * .5;
    tol1 = eps * (fabs(xx) + 1.) + tol3;
    tol2 = tol1 * 2.;

/*  check stopping criterion */

    *delta = fabs(xx - hh) + (bb - aa) * .5;
    /*     if (abs(xx-hh) .le. (tol2-half*(bb-aa))) goto 100 */
    if(verbose && OP.nopt > 1)
	REprintf(" .. DBG dopt() [%2d]:| %12g | %12g | %10.6e |\n",
		 OP.nopt, uu, fu, *delta);
    if (*delta <= tol2) {
	goto L_end;
    }
    if (OP.nopt >= OP.maxopt) {
	goto L_end;
    }
    /* Maybe another check :
     *     if (delpq <= EPSMAX*(one+pqnorm)) goto 100 */

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
    if (fabs(tt) >= fabs(ss * .5 * rr) || tt <= ss * (aa - xx) ||
	     tt >= ss * (bb - xx))
    { /*---  a golden-section step ---*/

	if (xx >= hh) {
	    ee = aa - xx;
	} else {
	    ee = bb - xx;
	}
	dd = cc * ee;

    }
    else { /*--- a parabolic-interpolation step ---*/

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
    ++OP.nopt;
    fu = pqopt(x, uu, w, iw, min_fnorm);

/*  update  aa, bb, vv, ww, and xx */

    if (fx >= fu) {
	if (uu >= xx) {
	    aa = xx;
	} else {
	    bb = xx;
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
	} else {
	    aa = uu;
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

L_end:
    *hood = -fx;
    filtfd_.cllf = *hood;
    return xx;
/* 900  format( i4, 2(1pe14.6), 1pe16.7, 1pe17.8, 1x, 2(i3))
 901  format( i4, 3(1pe10.2), 1pe11.2, 2(i3), 3(1pe8.1), i2) */
} /* dopt */

/* ****************************************************************************
******************************************************************************/

void fdcom(int *n, int *m, int *nar, int *nma,
	   double *hood, double *flmin, double *flmax,
	   double *epmin, double *epmax)
    /* is also called from R --> need all pointers */
{
/* Fill "parameter"s into global variables (Common blocks) needed later:
 *
 * copyright 1991 Department of Statistics, University of Washington
   written by Chris Fraley
 -----------------------------------------------------------------------------*/

    filtfd_.cllf = *hood;

/* machine constants */
    machfd_.fltmin = *flmin;
    machfd_.fltmax = *flmax;
    machfd_.epsmin = *epmin;
    machfd_.epsmax = *epmax;
    mauxfd_.epspt5 = sqrt(machfd_.epsmin);
    mauxfd_.epsp25 = sqrt(mauxfd_.epspt5);
    mauxfd_.epspt3 = pow(machfd_.epsmin, 0.3);
    mauxfd_.epsp75 = pow(machfd_.epsmin, 0.75);
    mauxfd_.bignum = 1. / machfd_.epsmin;

/* useful quantities -- integer "dimensions" : */
    Dims.n = *n;
    Dims.m = *m;
    Dims.p = *nar;
    Dims.q = *nma;
    Dims.pq = Dims.p + Dims.q;
    Dims.pq1 = Dims.pq + 1;
    if(Dims.p >= Dims.q) {
	Dims.maxpq = Dims.p;
	Dims.minpq = Dims.q;
    } else {
	Dims.maxpq = Dims.q;
	Dims.minpq = Dims.p;
    }
    Dims.maxpq1 = Dims.maxpq + 1;
    Dims.nm = *n - Dims.maxpq;

/* workspace allocation */
    w_opt.lqp = 1;
    w_fil.ly = w_opt.lqp + Dims.pq;
    w_fil.lamk = w_fil.ly;
    w_fil.lak = w_fil.lamk + *n;
    w_fil.lphi= w_fil.lak  + *n;
    w_fil.lvk = w_fil.lphi + *m; /* = lamk + 2*n + M = 1 + npq + 2n + M */
    w_fil.lpi = w_fil.lphi;
    w_opt.la  = w_fil.ly   + *n;
    w_opt.lajac = w_opt.la + *n - Dims.minpq;
    /* old   ipvt = lajac  +  max( (n-np)*np, (n-nq)*nq, (n-maxpq)*npq) */
    w_opt.ipvt = w_opt.lajac + (*n - Dims.maxpq) * Dims.pq;
    w_opt.ldiag= w_opt.ipvt + Dims.pq / 2 + 1;
    w_opt.lqtf = w_opt.ldiag + Dims.pq;
    w_opt.lwa1 = w_opt.lqtf + Dims.pq;
    w_opt.lwa2 = w_opt.lwa1 + Dims.pq;
    w_opt.lwa3 = w_opt.lwa2 + Dims.pq;
    w_opt.lwa4 = w_opt.lwa3 + Dims.pq;
/*      lfree  = lwa4   +  n - minpq */
    return;
} /* fdcom */



/**************************************************************************
 ************************************************************************** */
static
double pqopt(double *x, double d__, double *w, int *iw, double *min_fnorm)
{
    /* x: double x(n) */
    /* w: work array exactly as in main  fracdf() */

    /* 'const' (but need to pass pointers of these): */
    static int modelm = 1;
    static double factlm = 100.;

    /* Local variables */
    double t, u, slogvk;

    /* Parameter adjustments */
    --w;

    /* copyright 1991 Department of Statistics, University of Washington
     * written by Chris Fraley
 ---------------------------------------------------------------------------- */
    fdfilt(x, d__,
	   &w[(0 + (0 + (w_fil.ly   << 3))) / 8], &slogvk,
	   &w[(0 + (0 + (w_fil.lamk << 3))) / 8],
	   &w[(0 + (0 + (w_fil.lak  << 3))) / 8],
	   &w[(0 + (0 + (w_fil.lvk  << 3))) / 8],
	   &w[(0 + (0 + (w_fil.lphi << 3))) / 8],
	   &w[(0 + (0 + (w_fil.lpi  << 3))) / 8]);
    if (gammfd_.igamma != 0) {
	filtfd_.wnv  =  machfd_.fltmax;
	filtfd_.cllf = -machfd_.fltmax;
	warning("** pqopt() gamma error (%d)", gammfd_.igamma);
	return machfd_.fltmax;
    }
    t = (double) Dims.n;

    if (Dims.pq == 0) { /* trivial case ---  p = q = 0 : */

	filtfd_.wnv = F77_CALL(ddot)(&Dims.n,
				     &w[w_fil.ly], &ic__1,
				     &w[w_fil.ly], &ic__1) / t;
	OP.ifun = 0;
	OP.igrd = 0;
	OP.info = -1;
    }
    else {

/*     optimize as an unconstrained optimization problem */

	if (modelm == 2) {
	    F77_CALL(dcopy)(&Dims.pq, &c__1, &ic__0,
			    &w[w_opt.ldiag], &ic__1);
	}
	if (OP.nopt < 0) { // (never used ??)
	    REprintf("pqopt() -- nopt < 0 case --- should never happen.  Please report!");
	    if (Dims.p != 0) {
		int n_p = Dims.n - Dims.p;
		lmder1((U_fp)ajp_, n_p, Dims.p,
		       &w[w_opt.lqp + Dims.q], &w[w_opt.la], &w[w_opt.lajac], n_p,
		       TOL.f, TOL.x, TOL.g, OP.maxfun, &w[w_opt.ldiag],
		       modelm, factlm, &OP.info, &OP.ifun, &OP.igrd,
		       iw /* was &w[w_opt.ipvt] */, &w[w_opt.lqtf],
		       &w[w_opt.lwa1], &w[w_opt.lwa2], &w[w_opt.lwa3], &w[w_opt.lwa4],
		       &w[w_fil.ly]);
	    }
	    if (Dims.q != 0) {
		int n_q = Dims.n - Dims.q;
		lmder1((U_fp)ajq_, n_q, Dims.q,
		       &w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac], n_q,
		       TOL.f, TOL.x, TOL.g, OP.maxfun, &w[w_opt.ldiag],
		       modelm, factlm, &OP.info, &OP.ifun, &OP.igrd,
		       iw /* was &w[w_opt.ipvt] */, &w[w_opt.lqtf],
		       &w[w_opt.lwa1], &w[w_opt.lwa2], &w[w_opt.lwa3], &w[w_opt.lwa4],
		       &w[w_fil.ly]);
	    }
	}
	*min_fnorm =
	    lmder1((U_fp)ajqp_, Dims.nm, Dims.pq,
		   &w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac], Dims.nm,
		   TOL.f, TOL.x, TOL.g, OP.maxfun, &w[w_opt.ldiag],
		   modelm, factlm, &OP.info, &OP.ifun, &OP.igrd,
		   iw /* was &w[w_opt.ipvt] */, &w[w_opt.lqtf],
		   &w[w_opt.lwa1], &w[w_opt.lwa2], &w[w_opt.lwa3], &w[w_opt.lwa4],
		   &w[w_fil.ly]);

	if (OP.info == 0) { /* 'MINPACK : improper input parameters */
	    MinPck.iminpk = 10;
	    filtfd_.wnv = machfd_.fltmax;
	    filtfd_.cllf = -machfd_.fltmax;
	    return machfd_.fltmax;
	}

	if(OP.info== 5) MinPck.jminpk = 5; /* MINPACK : function evaluation limit reached */
	if(OP.info== 6) MinPck.jminpk = 6; /* MINPACK : ftol is too small */
	if(OP.info== 7) MinPck.jminpk = 7; /* MINPACK : xtol is too small */
	if(OP.info== 8) MinPck.jminpk = 8; /* MINPACK : gtol is too small */


/*     call daxpy( npq, (-one), w(lpq), 1, w(lqp), 1
     delpq  = sqrt(ddot( npq, w(lqp), 1, w(lqp), 1))
     pqnorm = sqrt(ddot( npq, w(lpq), 1, w(lpq), 1)) */

	filtfd_.wnv = *min_fnorm * *min_fnorm / (double) (Dims.nm - 1);
    }
    u = t * (log(filtfd_.wnv) + 2.8378) + slogvk;
    /* unused: BIC = u + (double) (Dims.p + Dims.q + 1) * log(t); */
    filtfd_.cllf = -u / 2.;
    return u / 2;
} /* End pqopt() */

/*************************************************************************** */

void
fdfilt(double *x, double d__,
       /* -> output */
       double *y, double *slogvk,
       /* using */
       double *amk, double *ak, double *vk,
       double *phi, double *pi)
{
/* called as	 fdfilt( x, d, w(ly), slogvk,
  		        w(lamk), w(lak), w(lvk), w(lphi), w(lpi))
     float             x(n)
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
 --   --> look at ./fdsim.c  which is similar (but simpler)
          and ../filters.R

**************************************************************************
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 -----------------------------------------------------------------------*/

    /* System generated locals */
    double d__1;

    /* Local variables */
    int j, k, km, mcap;
    double r__, s, t, u, v, z__, g0;

    /* Parameter adjustments */
    --pi;
    --phi;
    --vk;
    --ak;
    --amk;
    --y;
    --x;

    /* Function Body */
    mcap = imin2(Dims.m, Dims.n);

/* calculate amk(k), vk(k), and ak(k) for k=1,n (see W522-4 for notation). */


/*  k = 1 */

    amk[1] = 0.;
    ak[1] = 1.;

/*  k = 2 ;  initialize phi(1) */

    z__ = d__ / (1. - d__);
    amk[2] = z__ * x[1];
    ak[2] = 1. - z__;
    phi[1] = z__;
    d__1 = 1. - d__;
    t = dgamr_(&d__1);
    if (gammfd_.igamma != 0) {
	return;
    }
    d__1 = 1. - d__ * 2.;
    g0 = dgamma_(&d__1) * (t * t);
    if (gammfd_.igamma != 0) {
	return;
    }
    vk[1] = g0;
    vk[2] = g0 * (1. - z__ * z__);

/*  k = 3, mcap */

    for (k = 3; k <= mcap; ++k) {
	km = k - 1;
	t = (double) km;
	u = t - d__;

	/*  calculate phi() and vk() using the recursion formula on W498 */

	for (j = 1; j <= (km - 1); ++j) {
	    s = t - (double) j;
	    phi[j] *= t * (s - d__) / (u * s);
	}
	v = d__ / u;
	phi[km] = v;
	vk[k] = vk[km] * (1. - v * v);

	/*  form amk(k) and ak(k) */

	u = 0.;
	v = 1.;
	for (j = 1; j <= km; ++j) {
	    t = phi[j];
	    u += t * x[k - j];
	    v -= t;
	}
	amk[k] = u;
	ak[k] = v;
    }

/*     k = mcap+1, n */

    if (Dims.m < Dims.n) { /* i.e. mcap = min(M,n) != n */

      /* calculate pi(j), j = 1,mcap */

	pi[1] = d__;
	s = d__;
	for (j = 2; j <= mcap; ++j) {
	    u = (double) j;
	    t = pi[j - 1] * ((u - 1. - d__) / u);
	    s += t;
	    pi[j] = t;
	}
	s = 1. - s;
	r__ = 0.;
	u = (double) mcap;
	t = u * pi[mcap];

	for (k = mcap+1; k <= Dims.n; ++k) {
	    km = k - mcap;
	    z__ = 0.;
	    for (j = 1; j <= mcap; ++j) {
		z__ += pi[j] * x[k - j];
	    }
	    if (r__ == 0.) {
		amk[k] = z__;
		ak[k] = s;
	    } else {
		v = t * (1. - pow(u / k, d__)) / d__;
		amk[k] = z__ + v * r__ / ((double) km - 1.);
		ak[k] = s - v;
	    }
	    r__ += x[km];
	}
    }

/*  form muhat - see formula on W523. */

    r__ = 0.;
    s = 0.;
    for (k = 1; k <= (Dims.n); ++k) {
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
    filtfd_.hatmu = r__ / s;

/*  form filtered version */

    s = 0.;
    for (k = 1; k <= mcap; ++k)
	s += log(vk[k]);

    *slogvk = s;
    s = 0.;
    for (k = 1; k <= (Dims.n); ++k) {
	t = x[k] - amk[k] - filtfd_.hatmu * ak[k];
	if (k <= mcap)
	    t /= sqrt(vk[k]);

	s += t;
	y[k] = t;
    }
    if (Dims.pq == 0) {
	return;
    }
    t = (double) Dims.n;
    u = z__ / t;
    for (k = 1; k <= Dims.n; ++k)
	y[k] -= u;

    return;
} /* fdfilt */


/****************************************************************************
*****************************************************************************/
/* Passed to lmder1() minimizer, but also called from
 * hesdpq()  in ./fdhess.c : */
void ajqp_(double *qp, double *a, double *ajac, int lajac, int op_code, double *y)
{
    /* System generated locals */
    int ajac_dim1, ajac_offset;

    /* Local variables */
    static int i, k, l;
    static double s, t;
    static int km;

/*     double precision qp(npq), a(nm), ajac(nm,npq), y(n)
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 --------------------------------------------------------------------------
     Parameter adjustments */
    --qp;
    --a;
    ajac_dim1 = lajac;
    ajac_offset = 1 + ajac_dim1;
    ajac -= ajac_offset;
    --y;

    if (op_code == 1) { /*  objective calculation */

	for (k = Dims.maxpq1; k <= (Dims.n); ++k) {
	    km = k - Dims.maxpq;
	    t = 0.;
	    if (Dims.p != 0) {
		for (l = 1; l <= (Dims.p); ++l) {
		    t -= qp[Dims.q + l] * y[k - l];
		}
	    }
	    s = 0.;
	    if (Dims.q != 0) {
		for (l = 1; l <= (Dims.q); ++l) {
		    if (km <= l)
			break;
		    s += qp[l] * a[km - l];
		}
	    }

	    s = y[k] + (t + s);
	    if (fabs(s) <= mauxfd_.bignum) {
		a[km] = s;
	    } else {
		a[km] = sign(s) * mauxfd_.bignum;
	    }
	}
	++OP.nfun;

    } else if (op_code == 2) { /*  jacobian calculation */

	for (i = 1; i <= (Dims.pq); ++i) {
	    for (k = Dims.maxpq1; k <= (Dims.n); ++k) {
		km = k - Dims.maxpq;
		t = 0.;
		if (Dims.q != 0) {
		    for (l = 1; l <= (Dims.q); ++l) {
			if (km <= l)
			    break;
			t += qp[l] * ajac[km - l + i * ajac_dim1];
		    }
		}

		if (i <= Dims.q) {
		    if (km > i) {
			s = a[km - i] + t;
		    } else {
			s = t;
		    }
		} else {
		    s = -y[k - (i - Dims.q)] + t;
		}
		if (fabs(s) <= mauxfd_.bignum) {
		    ajac[km + i * ajac_dim1] = s;
		} else {
		    ajac[km + i * ajac_dim1] = sign(s) * mauxfd_.bignum;
		}
	    }
	}
	++OP.ngrd;
    }
    else { // invalid op_code
	warning("ajqp_(): invalid op_code = %d", op_code);
    }
    return;
} /* ajqp_ */

/****************************************************************************
 ****************************************************************************/

static void
ajp_(double *p, double *a, double *ajac, int lajac, int op_code, double *y)
    /*  p(np), a(nm), ajac(nm,npq), y(n) */
{
/* copyright 1991 Department of Statistics, University of Washington
   written by Chris Fraley
   -------------------------------------------------------------------------
   */

    /* Local variables */
    int i, k;

    /* Parameter adjustments */
    --p;
    --a;
    --y;

    /* Function Body */

    if (op_code == 1) { /*  objective calculation */

	if (Dims.p == 0) {

	}
	for (k = Dims.p + 1; k <= (Dims.n); ++k) {
	    double t = 0;
	    for (i = 1; i <= (Dims.p); ++i)
		t -= p[i] * y[k - i];

	    a[k - Dims.p] = y[k] + t;
	}
    }
    else if (op_code == 2) { /*  jacobian calculation */
	/* L200: */

	/* Matrix 1-indexing adjustments (System generated): */
	int ajac_dim1 =  lajac;
	ajac -= (1 + ajac_dim1);

	for (i = 1; i <= Dims.p; ++i)
	    for (k = Dims.p + 1; k <= (Dims.n); ++k)
		ajac[k - Dims.p + i * ajac_dim1] = - y[k - i];
    }
    return;
} /* ajp_
****************************************************************************
****************************************************************************/

static void
ajq_(double *qp, double *a, double *ajac, int lajac, int op_code, double *y)
     /*     double precision qp(npq), a(nm), ajac(nm,npq), y(n) */
{
/* copyright 1991 Department of Statistics, University of Washington
   written by Chris Fraley
   -------------------------------------------------------------------
*/

    /* Local variables */
    int i, k, l, km;
    double s, t;

    /* Parameter adjustments */
    --qp;
    --a;
    --y;

    if (op_code == 1) { /*---  objective calculation ---*/

	if (Dims.q == 0)
	    return;

	for (k = Dims.maxpq1; k <= (Dims.n); ++k) {
	    km = k - Dims.maxpq;
	    t = 0.;
	    if (Dims.p != 0) {
		for (l = 1; l <= (Dims.p); ++l) {
		    t -= qp[Dims.q + l] * y[k - l];
		}
	    }
	    s = 0.;
	    if (Dims.q != 0) {
		for (l = 1; l <= (Dims.q); ++l) {
		    if (km <= l)
			break;
		    s += qp[l] * a[km - l];
		}
	    }
	    a[km] = y[k] + (t + s);
	}
	++OP.nfun;
    }
    else if (op_code == 2) { /*---  jacobian calculation  ---*/
	/* L200: */

	/* Matrix 1-indexing adjustments (System generated): */
	int ajac_dim1 =  lajac;
	ajac -= (1 + ajac_dim1);

	for (i = 1; i <= (Dims.pq); ++i) {
	    for (k = Dims.maxpq1; k <= (Dims.n); ++k) {
		km = k - Dims.maxpq;
		t = 0.;
		if (Dims.q != 0) {
		    for (l = 1; l <= (Dims.q); ++l) {
			if (km <= l)
			    break;
			t += qp[l] * ajac[km - l + i * ajac_dim1];
		    }
		}
		if (i <= Dims.q) {
		    if (km > i) {
			ajac[km + i * ajac_dim1] = a[km - i] + t;
		    } else {
			ajac[km + i * ajac_dim1] = t;
		    }
		} else {
		    ajac[km + i * ajac_dim1] = -y[k - (i - Dims.q)] + t;
		}
	    }
	}
	++OP.ngrd;
    }
    return;
} /* ajq_ */

