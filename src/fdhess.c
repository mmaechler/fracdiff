/* fdhess.f -- translated by f2c (version 20031025).
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

/* ddot(), daxpy(), dcopy(), dscal() : */
#include <R_ext/BLAS.h>

/* dsvdc: */
#include <R_ext/Linpack.h>

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
    struct {
	int nn, mm, np, nq, npq, npq1, maxpq, maxpq1, minpq, nm;
    } _1;
    struct {
	int n, m, np, nq, npq, npq1, maxpq, maxpq1, minpq, nm;
    } _2;
} dimsfd_;

#define dimsfd_1 (dimsfd_._1)
#define dimsfd_2 (dimsfd_._2)

struct {
    double hatmu, wnv, cllf;
} filtfd_;

#define filtfd_1 filtfd_

struct {
    int ly, lamk, lak, lvk, lphi, lpi;
} wfilfd_;

#define wfilfd_1 wfilfd_

struct {
    int lqp, la, lajac, ipvt, ldiag, lqtf, lwa1, lwa2, lwa3, lwa4;
} woptfd_;

#define woptfd_1 woptfd_

struct {
    int igamma, jgamma;
} gammfd_;

#define gammfd_1 gammfd_

struct {
    int ksvd, kcov, kcor;
} hessfd_;

#define hessfd_1 hessfd_

/* Table of constant values */

static double c_b2 = .3;
static double c_b3 = .75;
static int c__1 = 1;
static int c__11 = 11;
static double c_b8 = 0.;
static int c__0 = 0;
static int c__2 = 2;
static double c_b78 = -1.;

/* ******************************************************************************
 ******************************************************************************
 Fill "parameter"s into global variables (Common blocks) called later:
 Subroutine */ int fdcom_(int *n, int *m, int *nar, int *
	nma, double *hood, double *flmin, double *flmax,
	double *epmin, double *epmax)
{
    /* Builtin functions */
    double sqrt(double), pow_dd(double *, double *);

/*  copyright 1991 Department of Statistics, University of Washington
  written by Chris Fraley
 ----------------------------------------------------------------------------- */
    filtfd_1.cllf = *hood;
/* machine constants */
    machfd_1.fltmin = *flmin;
    machfd_1.fltmax = *flmax;
    machfd_1.epsmin = *epmin;
    machfd_1.epsmax = *epmax;
    mauxfd_1.epspt5 = sqrt(machfd_1.epsmin);
    mauxfd_1.epsp25 = sqrt(mauxfd_1.epspt5);
    mauxfd_1.epspt3 = pow_dd(&machfd_1.epsmin, &c_b2);
    mauxfd_1.epsp75 = pow_dd(&machfd_1.epsmin, &c_b3);
    mauxfd_1.bignum = 1. / machfd_1.epsmin;
/* useful quantities */
    dimsfd_1.nn = *n;
    dimsfd_1.mm = *m;
    dimsfd_1.np = *nar;
    dimsfd_1.nq = *nma;
    dimsfd_1.npq = dimsfd_1.np + dimsfd_1.nq;
    dimsfd_1.npq1 = dimsfd_1.npq + 1;
    dimsfd_1.maxpq = max(dimsfd_1.np,dimsfd_1.nq);
    dimsfd_1.minpq = min(dimsfd_1.np,dimsfd_1.nq);
    dimsfd_1.maxpq1 = dimsfd_1.maxpq + 1;
    dimsfd_1.nm = *n - dimsfd_1.maxpq;
/* workspace allocation */
    woptfd_1.lqp = 1;
    wfilfd_1.ly = woptfd_1.lqp + dimsfd_1.npq;
    wfilfd_1.lamk = wfilfd_1.ly;
    wfilfd_1.lak = wfilfd_1.lamk + *n;
    wfilfd_1.lphi = wfilfd_1.lak + *n;
    wfilfd_1.lvk = wfilfd_1.lphi + *m;
/* 	= lamk  + 2*n +  M = 1 + npq + 2n + M */
    wfilfd_1.lpi = wfilfd_1.lphi;
    woptfd_1.la = wfilfd_1.ly + *n;
    woptfd_1.lajac = woptfd_1.la + *n - dimsfd_1.minpq;
/* old ipvt   = lajac  +  max( (n-np)*np, (n-nq)*nq, (n-maxpq)*npq) */
    woptfd_1.ipvt = woptfd_1.lajac + (*n - dimsfd_1.maxpq) * dimsfd_1.npq;
    woptfd_1.ldiag = woptfd_1.ipvt + dimsfd_1.npq / 2 + 1;
    woptfd_1.lqtf = woptfd_1.ldiag + dimsfd_1.npq;
    woptfd_1.lwa1 = woptfd_1.lqtf + dimsfd_1.npq;
    woptfd_1.lwa2 = woptfd_1.lwa1 + dimsfd_1.npq;
    woptfd_1.lwa3 = woptfd_1.lwa2 + dimsfd_1.npq;
    woptfd_1.lwa4 = woptfd_1.lwa3 + dimsfd_1.npq;
/*      lfree  = lwa4   +  n - minpq */
    return 0;
} /* fdcom_

 ******************************************************************************
 ******************************************************************************
 Subroutine */ int fdhpq_(double *x, double *h__, int *lh,

	double *w)
{
    /* System generated locals */
    int h_dim1, h_offset;

    /* Local variables */
    extern /* Subroutine */ int hesspq_(double *, double *,
	    double *, int *, double *, int *, double *,
	    double *);

/*     float		x(n)
     double precision	H(lH, npq1)
  copyright 1991 Department of Statistics, University of Washington
  written by Chris Fraley
 -----------------------------------------------------------------------------
     Parameter adjustments */
    --x;
    h_dim1 = *lh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --w;

    /* Function Body */
    hesspq_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
	    dimsfd_2.nm, &h__[h_offset], lh, &w[woptfd_1.lwa4], &w[
	    woptfd_1.lwa1]);
/*     call dcopy( npq1, zero, 0, H(1,1), lH) */
/*     call dcopy( npq , zero, 0, H(2,1), 1) */
    return 0;
} /* fdhpq_

 ******************************************************************************
 ******************************************************************************
 Subroutine */ int fdcov_(double *x, double *d__, double *hh,

	double *hd, double *cov, int *lcov, double *cor,
	int *lcor, double *se, double *w, int *info)
{
    /* System generated locals */
    int cov_dim1, cov_offset, cor_dim1, cor_offset, i__1, i__2;
    double d__1, d__2;

    /* Builtin functions */
    double sqrt(double);

    /* Local variables */
    static int i__, j, k, le, ls, lu, lv;
    static int lwork;
    static double temp;
    extern /* Subroutine */ int hesdpq_(double *, double *,
	    double *, double *, double *), invsvd_(double *,
	    double *, int *, double *, int *, double *,
	    int *);

/*     float               x(n)
     double precision   d, hh, hd(npq1), cov(lcov,npq1),
    *                   cor(lcor,npq1), se(npq1)
  copyright 1991 Department of Statistics, University of Washington
  written by Chris Fraley
 -----------------------------------------------------------------------------
     Parameter adjustments */
    --x;
    --hd;
    cov_dim1 = *lcov;
    cov_offset = 1 + cov_dim1;
    cov -= cov_offset;
    cor_dim1 = *lcor;
    cor_offset = 1 + cor_dim1;
    cor -= cor_offset;
    --se;
    --w;

    /* Function Body */
    hesdpq_(&x[1], d__, hh, &hd[1], &w[1]);
    F77_CALL(dcopy)(&dimsfd_2.npq1, &hd[1], &c__1, &cov[cov_offset], lcov);
    gammfd_1.igamma = 0;
    gammfd_1.jgamma = 0;
    hessfd_1.ksvd = 0;
    hessfd_1.kcov = 0;
    hessfd_1.kcor = 0;
    *info = 0;
    temp = 1.;
    i__1 = dimsfd_2.npq1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = dimsfd_2.npq1;
	for (j = i__ + 1; j <= i__2; ++j) {
	    cov[j + i__ * cov_dim1] = cov[i__ + j * cov_dim1];
	}
    }
    ls = wfilfd_1.ly;
    lu = ls + dimsfd_2.npq1 + 1;
    lv = lu + dimsfd_2.npq1 * dimsfd_2.npq1;
    le = lv + dimsfd_2.npq1 * dimsfd_2.npq1;
    lwork = le + dimsfd_2.npq1;
/*      lfree = lwork + npq1 */
    F77_CALL(dsvdc)(&cov[cov_offset], lcov, &dimsfd_2.npq1, &dimsfd_2.npq1, &w[ls], &w[
	    le], &w[lu], &dimsfd_2.npq1, &w[lv], &dimsfd_2.npq1, &w[lwork], &
	    c__11, info);
    if (*info != 0) {
	F77_CALL(dcopy)(&dimsfd_2.npq1, &c_b8, &c__0, &se[1], &c__1);
	i__1 = dimsfd_2.npq1;
	for (j = 1; j <= i__1; ++j) {
	    F77_CALL(dcopy)(&dimsfd_2.npq1, &c_b8, &c__0, &cov[j * cov_dim1 + 1], &
		    c__1);
	}
	hessfd_1.ksvd = 1;
	*info = 3;
	return 0;
    }
    invsvd_(&w[ls], &w[lu], &dimsfd_2.npq1, &w[lv], &dimsfd_2.npq1, &cov[
	    cov_offset], lcov);
    i__1 = dimsfd_2.npq1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = dimsfd_2.npq1;
	for (j = i__ + 1; j <= i__2; ++j) {
	    cov[j + i__ * cov_dim1] = cov[i__ + j * cov_dim1];
	}
    }
    temp = 1.;
    i__1 = dimsfd_2.npq1;
    for (j = 1; j <= i__1; ++j) {
	if (cov[j + j * cov_dim1] > 0.) {
	    se[j] = sqrt(cov[j + j * cov_dim1]);
	} else {
/* Computing MIN */
	    d__1 = temp, d__2 = cov[j + j * cov_dim1];
	    temp = min(d__1,d__2);
	    se[j] = 0.;
	}
    }
    if (temp == 1.) {
	i__1 = dimsfd_2.npq1;
	for (k = 1; k <= i__1; ++k) {
	    F77_CALL(dcopy)(&k, &cov[k * cov_dim1 + 1], &c__1, &cor[k * cor_dim1 + 1],
		    &c__1);
	}
	i__1 = dimsfd_2.npq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    i__2 = dimsfd_2.npq1 - i__ + 1;
	    d__1 = 1. / se[i__];
	    F77_CALL(dscal)(&i__2, &d__1, &cor[i__ + i__ * cor_dim1], lcor);
	}
	i__1 = dimsfd_2.npq1;
	for (j = 1; j <= i__1; ++j) {
	    d__1 = 1. / se[j];
	    F77_CALL(dscal)(&j, &d__1, &cor[j * cor_dim1 + 1], &c__1);
	}
    } else {
	hessfd_1.kcor = 1;
	i__1 = dimsfd_2.npq1;
	for (j = 1; j <= i__1; ++j) {
	    F77_CALL(dcopy)(&dimsfd_2.npq1, &c_b8, &c__0, &cor[j * cor_dim1 + 1], &
		    c__1);
	}
    }
    i__1 = dimsfd_2.npq1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = dimsfd_2.npq1;
	for (j = i__ + 1; j <= i__2; ++j) {
	    cor[j + i__ * cor_dim1] = cor[i__ + j * cor_dim1];
	}
    }
    if (gammfd_1.igamma != 0) {
	*info = 4;
    }
    if (gammfd_1.jgamma != 0) {
	*info = 1;
    }
    if (hessfd_1.ksvd != 0) {
	*info = 3;
    }
    if (hessfd_1.kcov != 0) {
	*info = 2;
    }
    if (hessfd_1.kcor != 0) {
	*info = 3;
    }
    return 0;
} /* fdcov_

 ******************************************************************************
 ******************************************************************************
 Subroutine */ int invsvd_(double *s, double *u, int *lu,

	double *v, int *lv, double *cov, int *lcov)
{
    /* System generated locals */
    int u_dim1, u_offset, v_dim1, v_offset, cov_dim1, cov_offset, i__1,
	    i__2;
    double d__1;

    /* Local variables */
    static int i__, j, k;
    static double ss;
    static int krank;

/*     double precision   s(npq1), u(lu,npq1), v(lv,npq1), cov(lcov,npq1)
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 -----------------------------------------------------------------------------
     Parameter adjustments */
    --s;
    u_dim1 = *lu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *lv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    cov_dim1 = *lcov;
    cov_offset = 1 + cov_dim1;
    cov -= cov_offset;

    /* Function Body */
    krank = dimsfd_2.npq1;
    i__1 = dimsfd_2.npq1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ss = s[i__];
	i__2 = dimsfd_2.npq1;
	for (j = 1; j <= i__2; ++j) {
	    if (ss < 1.) {
		if ((d__1 = u[i__ + j * u_dim1], abs(d__1)) > ss *
			machfd_1.fltmax) {
		    krank = i__ - 1;
		    hessfd_1.kcov = 1;
		    goto L100;
		}
	    }
	}
    }
L100:
    i__1 = dimsfd_2.npq1;
    for (k = 1; k <= i__1; ++k) {
	F77_CALL(dcopy)(&k, &c_b8, &c__0, &cov[k * cov_dim1 + 1], &c__1);
    }
    if (krank == 0) {
	return 0;
    }
/*      do k = 1, npq1 */
/*        do i = 1, npq1 */
/*          do j = i, npq1 */
/*            H(i,j) =  H(i,j) + s(k)*u(i,k)*v(j,k) */
/*          end do */
/*        end do */
/*      end do */
/*      do k = 1, npq1 */
/*        ss = s(k) */
/*        do j = 1, npq1 */
/*          call daxpy( j, ss*v(j,k), u(1,k), 1, H(1,j), 1) */
/*        end do */
/*      end do */
    i__1 = krank;
    for (k = 1; k <= i__1; ++k) {
	ss = -1. / s[k];
	i__2 = dimsfd_2.npq1;
	for (j = 1; j <= i__2; ++j) {
	    d__1 = ss * u[j + k * u_dim1];
	    F77_CALL(daxpy)(&j, &d__1, &v[k * v_dim1 + 1], &c__1, &cov[j * cov_dim1 +
		    1], &c__1);
	}
    }
    return 0;
} /* invsvd_

 ******************************************************************************
 ******************************************************************************
 Subroutine */ int hesspq_(double *qp, double *a, double *ajac,

	int *lajac, double *h__, int *lh, double *aij,
	double *g)
{
    /* System generated locals */
    int ajac_dim1, ajac_offset, h_dim1, h_offset, i__1, i__2, i__3, i__4;

    /* Local variables */
    static int i__, j, k, l;
    static double s, t, u;
    static int km;
    static double fac;


/*     double precision	qp(npq), a(nm), ajac(nm,npq)
     double precision	H(lH,npq1), aij(nm), g(npq)
 analytic Hessian with respect to p and q variables
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 -----------------------------------------------------------------------------
     Parameter adjustments */
    --qp;
    --a;
    ajac_dim1 = *lajac;
    ajac_offset = 1 + ajac_dim1;
    ajac -= ajac_offset;
    h_dim1 = *lh;
    h_offset = 1 + h_dim1;
    h__ -= h_offset;
    --aij;
    --g;

    /* Function Body */
    fac = 1. / (filtfd_1.wnv * (double) (dimsfd_2.nm - 1));
    if (dimsfd_2.nq != 0 && dimsfd_2.np != 0) {
	i__1 = dimsfd_2.npq;
	for (k = 1; k <= i__1; ++k) {
	    g[k] = F77_CALL(ddot)(&dimsfd_2.nm, &a[1], &c__1, &ajac[k * ajac_dim1 + 1],
		     &c__1);
	}
	i__1 = dimsfd_2.np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    u = g[dimsfd_2.nq + i__];
	    i__2 = dimsfd_2.nq;
	    for (j = 1; j <= i__2; ++j) {
		u = g[j] * u;
		i__3 = dimsfd_2.n;
		for (k = dimsfd_2.maxpq1; k <= i__3; ++k) {
		    km = k - dimsfd_2.maxpq;
		    t = 0.;
		    i__4 = dimsfd_2.nq;
		    for (l = 1; l <= i__4; ++l) {
			if (km <= l) {
			    goto L301;
			}
			t += qp[l] * aij[km - l];
		    }
L301:
		    if (km > j) {
			aij[km] = ajac[km - j + (dimsfd_2.nq + i__) *
				ajac_dim1] + t;
		    } else {
			aij[km] = t;
		    }
		}
		s = F77_CALL(ddot)(&dimsfd_2.nm, &ajac[(dimsfd_2.nq + i__) * ajac_dim1
			+ 1], &c__1, &ajac[j * ajac_dim1 + 1], &c__1);
		t = F77_CALL(ddot)(&dimsfd_2.nm, &a[1], &c__1, &aij[1], &c__1);
		h__[i__ + 1 + (dimsfd_2.np + j + 1) * h_dim1] = -((double)
			 dimsfd_2.n) * (s + t - fac * 2. * u) * fac;
	    }
	}
    }
    if (dimsfd_2.nq != 0) {
	i__1 = dimsfd_2.nq;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    u = g[i__];
	    i__2 = dimsfd_2.nq;
	    for (j = i__; j <= i__2; ++j) {
		u = g[j] * u;
		i__3 = dimsfd_2.n;
		for (k = dimsfd_2.maxpq1; k <= i__3; ++k) {
		    km = k - dimsfd_2.maxpq;
		    t = 0.;
		    i__4 = dimsfd_2.nq;
		    for (l = 1; l <= i__4; ++l) {
			if (km <= l) {
			    goto L302;
			}
			t += qp[l] * aij[km - l];
		    }
L302:
		    s = 0.;
		    if (km > i__) {
			s += ajac[km - i__ + j * ajac_dim1];
		    }
		    if (km > j) {
			s += ajac[km - j + i__ * ajac_dim1];
		    }
		    aij[km] = s + t;
		}
		s = F77_CALL(ddot)(&dimsfd_2.nm, &ajac[i__ * ajac_dim1 + 1], &c__1, &
			ajac[j * ajac_dim1 + 1], &c__1);
		t = F77_CALL(ddot)(&dimsfd_2.nm, &a[1], &c__1, &aij[1], &c__1);
		h__[dimsfd_2.np + i__ + 1 + (dimsfd_2.np + j + 1) * h_dim1] =
			-((double) dimsfd_2.n) * (s + t - fac * 2. * u) *
			fac;
	    }
	}
    }
    if (dimsfd_2.np != 0) {
	i__1 = dimsfd_2.np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    u = g[dimsfd_2.nq + i__];
	    i__2 = dimsfd_2.np;
	    for (j = i__; j <= i__2; ++j) {
		u = g[dimsfd_2.nq + j] * u;
/*            do k = maxpq1, n */
/*              km  =  k - maxpq */
/*              t  = zero */
/*              if (nq .ne. 0) then */
/*               do l = 1, nq */
/*                  if (km .le. l) goto 303 */
/*                  t  = t + qp(l)*aij(km-l) */
/*               end do */
/*              end if */
/* 303          continue */
/*              aij(km) = t */
/*            end do */
		s = F77_CALL(ddot)(&dimsfd_2.nm, &ajac[(dimsfd_2.nq + i__) * ajac_dim1
			+ 1], &c__1, &ajac[(dimsfd_2.nq + j) * ajac_dim1 + 1],
			 &c__1);
/*            t = ddot( nm, a             , 1, aij           , 1) */
/*            H(i+1,j+1) = -dble(n)*((s + t) - two*fac*u)*fac */
		h__[i__ + 1 + (j + 1) * h_dim1] = -((double) dimsfd_2.n) *
			 (s - fac * 2. * u) * fac;
	    }
	}
    }
    return 0;
} /* hesspq_

     hesspq
 ******************************************************************************
 ******************************************************************************
 Subroutine */ int hesdpq_(double *x, double *d__, double *hh,

	double *hd, double *w)
{
    /* System generated locals */
    double d__1;

    /* Builtin functions */
    double log(double);

    /* Local variables */
    static double fa, fb;

    extern /* Subroutine */ int ajqp_(double *, double *, double *
	    , int *, int *, double *),
    fdfilt_(double *, double *, double *, double *,
	    double *, double *, double *, double *,
	    double *),
    gradpq_(double *, double *, double *, int *);
    static double slogvk;

/*     float		 x(n)
     double precision	 d, hh, hd(npq1), w(*)
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 -----------------------------------------------------------------------------
     Parameter adjustments */
    --w;
    --hd;
    --x;

    /* Function Body */
    if (*hh <= 0.) {
	*hh = (abs(filtfd_1.cllf) + 1.) * mauxfd_1.epspt5;
    }
    *hh = min(*hh,.1);
    if (*d__ - *hh > 0.) {
	d__1 = *d__ - *hh;
	fdfilt_(&x[1], &d__1, &w[wfilfd_1.ly], &slogvk, &w[wfilfd_1.lamk], &w[
		wfilfd_1.lak], &w[wfilfd_1.lvk], &w[wfilfd_1.lphi], &w[
		wfilfd_1.lpi]);
	if (dimsfd_2.npq != 0) {
	    ajqp_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
		    dimsfd_2.nm, &c__1, &w[wfilfd_1.ly]);
	    ajqp_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
		    dimsfd_2.nm, &c__2, &w[wfilfd_1.ly]);
	    gradpq_(&w[woptfd_1.lwa1], &w[woptfd_1.la], &w[woptfd_1.lajac], &
		    dimsfd_2.nm);
	    filtfd_1.wnv = F77_CALL(ddot)(&dimsfd_2.nm, &w[woptfd_1.la], &c__1, &w[
		    woptfd_1.la], &c__1);
	    d__1 = 1. / filtfd_1.wnv;
	    F77_CALL(dscal)(&dimsfd_2.npq, &d__1, &w[woptfd_1.lwa1], &c__1);
	    filtfd_1.wnv /= (double) (dimsfd_2.nm - 1);
	} else {
	    filtfd_1.wnv = F77_CALL(ddot)(&dimsfd_2.nm, &w[wfilfd_1.ly], &c__1, &w[
		    wfilfd_1.ly], &c__1) / (double) (dimsfd_2.nm - 1);
	}
	fa = -((double) dimsfd_2.n * (log(filtfd_1.wnv) + 2.8378) +
		slogvk) / 2.;
	if (*d__ + *hh < .5) {
	    d__1 = *d__ + *hh;
	    fdfilt_(&x[1], &d__1, &w[wfilfd_1.ly], &slogvk, &w[wfilfd_1.lamk],
		     &w[wfilfd_1.lak], &w[wfilfd_1.lvk], &w[wfilfd_1.lphi], &
		    w[wfilfd_1.lpi]);
	    if (dimsfd_2.npq != 0) {
		ajqp_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
			dimsfd_2.nm, &c__1, &w[wfilfd_1.ly]);
		ajqp_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
			dimsfd_2.nm, &c__2, &w[wfilfd_1.ly]);
		gradpq_(&w[woptfd_1.lwa2], &w[woptfd_1.la], &w[woptfd_1.lajac]
			, &dimsfd_2.nm);
		filtfd_1.wnv = F77_CALL(ddot)(&dimsfd_2.nm, &w[woptfd_1.la], &c__1, &w[
			woptfd_1.la], &c__1);
		d__1 = 1. / filtfd_1.wnv;
		F77_CALL(dscal)(&dimsfd_2.npq, &d__1, &w[woptfd_1.lwa2], &c__1);
		filtfd_1.wnv /= (double) (dimsfd_2.nm - 1);
	    } else {
		filtfd_1.wnv = F77_CALL(ddot)(&dimsfd_2.nm, &w[wfilfd_1.ly], &c__1, &w[
			wfilfd_1.ly], &c__1) / (double) (dimsfd_2.nm - 1);
	    }
	    fb = -((double) dimsfd_2.n * (log(filtfd_1.wnv) + 2.8378) +
		    slogvk) / 2.;
	    hd[1] = (fa + fb - filtfd_1.cllf * 2.) / (*hh * *hh);
	} else {
	    d__1 = *d__ - *hh * 2.;
	    fdfilt_(&x[1], &d__1, &w[wfilfd_1.ly], &slogvk, &w[wfilfd_1.lamk],
		     &w[wfilfd_1.lak], &w[wfilfd_1.lvk], &w[wfilfd_1.lphi], &
		    w[wfilfd_1.lpi]);
	    if (dimsfd_2.npq != 0) {
		ajqp_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
			dimsfd_2.nm, &c__1, &w[wfilfd_1.ly]);
		ajqp_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
			dimsfd_2.nm, &c__2, &w[wfilfd_1.ly]);
		gradpq_(&w[woptfd_1.lwa2], &w[woptfd_1.la], &w[woptfd_1.lajac]
			, &dimsfd_2.nm);
		filtfd_1.wnv = F77_CALL(ddot)(&dimsfd_2.nm, &w[woptfd_1.la], &c__1, &w[
			woptfd_1.la], &c__1);
		d__1 = 1. / filtfd_1.wnv;
		F77_CALL(dscal)(&dimsfd_2.npq, &d__1, &w[woptfd_1.lwa2], &c__1);
		filtfd_1.wnv /= (double) (dimsfd_2.nm - 1);
	    } else {
		filtfd_1.wnv = F77_CALL(ddot)(&dimsfd_2.nm, &w[wfilfd_1.ly], &c__1, &w[
			wfilfd_1.ly], &c__1) / (double) (dimsfd_2.nm - 1);
	    }
	    fb = -((double) dimsfd_2.n * (log(filtfd_1.wnv) + 2.8378) +
		    slogvk) / 2.;
	    hd[1] = (filtfd_1.cllf + fb - fa * 2.) / (*hh * 2. * *hh);
	}
    } else {
	d__1 = *d__ + *hh;
	fdfilt_(&x[1], &d__1, &w[wfilfd_1.ly], &slogvk, &w[wfilfd_1.lamk], &w[
		wfilfd_1.lak], &w[wfilfd_1.lvk], &w[wfilfd_1.lphi], &w[
		wfilfd_1.lpi]);
	if (dimsfd_2.npq != 0) {
	    ajqp_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
		    dimsfd_2.nm, &c__1, &w[wfilfd_1.ly]);
	    ajqp_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
		    dimsfd_2.nm, &c__2, &w[wfilfd_1.ly]);
	    gradpq_(&w[woptfd_1.lwa1], &w[woptfd_1.la], &w[woptfd_1.lajac], &
		    dimsfd_2.nm);
	    filtfd_1.wnv = F77_CALL(ddot)(&dimsfd_2.nm, &w[woptfd_1.la], &c__1, &w[
		    woptfd_1.la], &c__1);
	    d__1 = 1. / filtfd_1.wnv;
	    F77_CALL(dscal)(&dimsfd_2.npq, &d__1, &w[woptfd_1.lwa1], &c__1);
	    filtfd_1.wnv /= (double) (dimsfd_2.nm - 1);
	} else {
	    filtfd_1.wnv = F77_CALL(ddot)(&dimsfd_2.nm, &w[wfilfd_1.ly], &c__1, &w[
		    wfilfd_1.ly], &c__1) / (double) (dimsfd_2.nm - 1);
	}
	fa = -((double) dimsfd_2.n * (log(filtfd_1.wnv) + 2.8378) +
		slogvk) / 2.;
	d__1 = *d__ + *hh * 2.;
	fdfilt_(&x[1], &d__1, &w[wfilfd_1.ly], &slogvk, &w[wfilfd_1.lamk], &w[
		wfilfd_1.lak], &w[wfilfd_1.lvk], &w[wfilfd_1.lphi], &w[
		wfilfd_1.lpi]);
	if (dimsfd_2.npq != 0) {
	    ajqp_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
		    dimsfd_2.nm, &c__1, &w[wfilfd_1.ly]);
	    ajqp_(&w[woptfd_1.lqp], &w[woptfd_1.la], &w[woptfd_1.lajac], &
		    dimsfd_2.nm, &c__2, &w[wfilfd_1.ly]);
	    gradpq_(&w[woptfd_1.lwa1], &w[woptfd_1.la], &w[woptfd_1.lajac], &
		    dimsfd_2.nm);
	    filtfd_1.wnv = F77_CALL(ddot)(&dimsfd_2.nm, &w[woptfd_1.la], &c__1, &w[
		    woptfd_1.la], &c__1);
	    d__1 = 1. / filtfd_1.wnv;
	    F77_CALL(dscal)(&dimsfd_2.npq, &d__1, &w[woptfd_1.lwa1], &c__1);
	    filtfd_1.wnv /= (double) (dimsfd_2.nm - 1);
	} else {
	    filtfd_1.wnv = F77_CALL(ddot)(&dimsfd_2.nm, &w[wfilfd_1.ly], &c__1, &w[
		    wfilfd_1.ly], &c__1) / (double) (dimsfd_2.nm - 1);
	}
	fb = -((double) dimsfd_2.n * (log(filtfd_1.wnv) + 2.8378) +
		slogvk) / 2.;
	hd[1] = (filtfd_1.cllf + fb - fa * 2.) / (*hh * 2. * *hh);
    }
    if (dimsfd_2.npq == 0) {
	return 0;
    }
    F77_CALL(daxpy)(&dimsfd_2.npq, &c_b78, &w[woptfd_1.lwa2], &c__1, &w[woptfd_1.lwa1],
	     &c__1);
    d__1 = (double) dimsfd_2.n / (*hh * 2.);
    F77_CALL(dscal)(&dimsfd_2.npq, &d__1, &w[woptfd_1.lwa1], &c__1);
    F77_CALL(dcopy)(&dimsfd_2.npq, &w[woptfd_1.lwa1], &c__1, &hd[2], &c__1);
    return 0;
} /* hesdpq_

     hesdpq
 ******************************************************************************
 ******************************************************************************
 Subroutine */ int gradpq_(double *g, double *a, double *ajac,

	int *ljac)
{
    /* System generated locals */
    int ajac_dim1, ajac_offset, i__1;

    /* Local variables */
    static int i__, j;

/*     double precision	g(npq), a(nm), ajac(nm,npq)
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 ------------------------------------------------------------------------------
     Parameter adjustments */
    --g;
    --a;
    ajac_dim1 = *ljac;
    ajac_offset = 1 + ajac_dim1;
    ajac -= ajac_offset;

    /* Function Body */
    if (dimsfd_2.np != 0) {
	i__1 = dimsfd_2.np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    g[i__] = F77_CALL(ddot)(&dimsfd_2.nm, &a[1], &c__1, &ajac[(dimsfd_2.nq +
		    i__) * ajac_dim1 + 1], &c__1);
	}
    }
    if (dimsfd_2.nq != 0) {
	i__1 = dimsfd_2.nq;
	for (j = 1; j <= i__1; ++j) {
	    g[dimsfd_2.np + j] = F77_CALL(ddot)(&dimsfd_2.nm, &a[1], &c__1, &ajac[j *
		    ajac_dim1 + 1], &c__1);
	}
    }
    return 0;
} /* gradpq_ */

