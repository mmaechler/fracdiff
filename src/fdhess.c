/*-*- mode: C; kept-old-versions: 12;  kept-new-versions: 20; -*-
 *
 * fdhess.f -- translated by f2c (version 20031025).
 * and produced by  f2c-clean,v 1.10 2002/03/28 16:37:27 maechler
 *
 * and manually pretty edited by Martin Maechler, 2004-10-01
*/

#include <Rmath.h>

#include "fracdiff.h"

/* ddot(), daxpy(), dcopy(), dscal() : */
#include <R_ext/BLAS.h>

/* dsvdc: */
#include <R_ext/Linpack.h>

#ifndef min
# define	min(a, b)		((a) > (b) ? (b) : (a))
#endif

/* called from R : */

void fdcov(double *x, double *d__, double *hh,
	   double *hd, double *cov, int *lcov, double *cor,
	   int *lcor, double *se, double *w, int *info);

void fdhpq(/*double *x, */double *h__, int *lh, double *w);

/* local to this file: */
static
int hesdpq_(double *, double *,
	    double *, double *, double *);
static
void hesspq_(double *qp, double *a, double *ajac,
	     int *lajac, double *h__, int *lh, double *aij, double *g);

static
int invsvd_(double *, double *, int *,
	    double *, int *, double *, int *);

static
void gradpq(double *g, double a[], double ajac[], int l_ajac);


/* Common Block Declarations --- included as "extern" */
#define FD_EXTERNAL extern

#include "mach_comm.h"
#include "maux_comm.h"
#include "Dims_comm.h"

#include "filt_comm.h"
#include "wfil_comm.h"
#include "wopt_comm.h"

#include "gamm_comm.h"
#include "hess_comm.h"


/* Table of constant values */
static int c__0 = 0;
static int c__1 = 1;
static int c__2 = 2;
static int c__11 = 11;

static double c_0d = 0.;
static double c_m1 = -1.;

/*******************************************************************************
 *******************************************************************************/

void fdhpq(/*double *x, */double *h__, int *lh, double *w)
{
/* double precision	H(lH, npq1)
*/

/*  copyright 1991 Department of Statistics, University of Washington
  written by Chris Fraley
 -----------------------------------------------------------------------------
     Parameter adjustments */
    --w;

    /* Function Body */
    hesspq_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac], &Dims.nm,
	    h__, lh, &w[woptfd_.lwa4], &w[woptfd_.lwa1]);
/*     call dcopy( npq1, zero, 0, H(1,1), lH) */
/*     call dcopy( npq , zero, 0, H(2,1), 1) */
    return;
} /* fdhpq */

/*******************************************************************************
 ****************************************************************************** */

void fdcov(double *x, double *d__, double *hh,
	   double *hd, double *cov, int *lcov, double *cor,
	   int *lcor, double *se, double *w, int *info)
{
    /* System generated locals */
    int cov_dim1, cov_offset, cor_dim1, cor_offset, npq1, i__2;
    double d__1;

    /* Local variables */
    int i__, j, k, le, ls, lu, lv, lwork;
    double temp;

/*     float		   x(n)
     double precision	d, hh, hd(npq1), cov(lcov,npq1),
    *			cor(lcor,npq1), se(npq1)
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
    npq1 = Dims.npq1;

    hesdpq_(&x[1], d__, hh, &hd[1], &w[1]);
    F77_CALL(dcopy)(&npq1, &hd[1], &c__1, &cov[cov_offset], lcov);

    gammfd_.igamma = 0;
    gammfd_.jgamma = 0;
    hessfd_.ksvd = 0;
    hessfd_.kcov = 0;
    hessfd_.kcor = 0;
    *info = 0;
    temp = 1.;

    for (i__ = 1; i__ <= npq1; ++i__) {
	for (j = i__ + 1; j <= npq1; ++j) {
	    cov[j + i__ * cov_dim1] = cov[i__ + j * cov_dim1];
	}
    }
    ls = wfilfd_.ly;
    lu = ls + npq1 + 1;
    lv = lu + npq1 * npq1;
    le = lv + npq1 * npq1;
    lwork = le + npq1;
/*	lfree = lwork + npq1 */
    F77_CALL(dsvdc)(&cov[cov_offset], lcov, &npq1, &npq1, &w[ls],
		    &w[le], &w[lu], &npq1, &w[lv], &npq1, &w[lwork],
		    &c__11, info);
    if (*info != 0) {
	F77_CALL(dcopy)(&npq1, &c_0d, &c__0, &se[1], &c__1);
	for (j = 1; j <= npq1; ++j) {
	    F77_CALL(dcopy)(&npq1, &c_0d, &c__0,
			    &cov[j * cov_dim1 + 1], & c__1);
	}
	hessfd_.ksvd = 1;
	*info = 3;
	return;
    }
    invsvd_(&w[ls], &w[lu], &npq1, &w[lv], &npq1, &cov[
	    cov_offset], lcov);
    for (i__ = 1; i__ <= npq1; ++i__) {
	for (j = i__ + 1; j <= npq1; ++j) {
	    cov[j + i__ * cov_dim1] = cov[i__ + j * cov_dim1];
	}
    }
    temp = 1.;
    for (j = 1; j <= npq1; ++j) {
	if (cov[j + j * cov_dim1] > 0.) {
	    se[j] = sqrt(cov[j + j * cov_dim1]);
	} else {
	    temp = fmin2(temp, cov[j + j * cov_dim1]);
	    se[j] = 0.;
	}
    }
    if (temp == 1.) {
	for (k = 1; k <= npq1; ++k) {
	    F77_CALL(dcopy)(&k, &cov[k * cov_dim1 + 1], &c__1,
				&cor[k * cor_dim1 + 1], &c__1);
	}
	for (i__ = 1; i__ <= npq1; ++i__) {
	    i__2 = npq1 - i__ + 1;
	    d__1 = 1. / se[i__];
	    F77_CALL(dscal)(&i__2, &d__1, &cor[i__ + i__ * cor_dim1], lcor);
	}
	for (j = 1; j <= npq1; ++j) {
	    d__1 = 1. / se[j];
	    F77_CALL(dscal)(&j, &d__1, &cor[j * cor_dim1 + 1], &c__1);
	}
    } else {
	hessfd_.kcor = 1;
	for (j = 1; j <= npq1; ++j) {
	    F77_CALL(dcopy)(&npq1, &c_0d, &c__0,
			    &cor[j * cor_dim1 + 1], &c__1);
	}
    }
    for (i__ = 1; i__ <= npq1; ++i__)
	for (j = i__ + 1; j <= npq1; ++j)
	    cor[j + i__ * cor_dim1] = cor[i__ + j * cor_dim1];

    if (gammfd_.igamma != 0) *info = 4;
    if (gammfd_.jgamma != 0) *info = 1;
    if (hessfd_.ksvd != 0)   *info = 3;
    if (hessfd_.kcov != 0)   *info = 2;
    if (hessfd_.kcor != 0)   *info = 3;
    return;
} /* fdcov */

/******************************************************************************
 ******************************************************************************
 Subroutine */ int
invsvd_(double *s, double *u, int *lu,
	double *v, int *lv, double *cov, int *lcov)
{
/*     double precision   s(npq1), u(lu,npq1), v(lv,npq1), cov(lcov,npq1)
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 ---------------------------------------------------------------------------*/

    /* System generated locals */
    int u_dim1, u_offset, v_dim1, v_offset, cov_dim1, cov_offset;
    double d__1;

    /* Local variables */
    int i__, j, k, krank, npq1 = Dims.npq1;
    double ss;

    /* Parameter adjustments */
    --s;
    u_dim1 = *lu;	u_offset = 1 + u_dim1;		  u -= u_offset;
    v_dim1 = *lv;	v_offset = 1 + v_dim1;		  v -= v_offset;
    cov_dim1 = *lcov; cov_offset = 1 + cov_dim1;	cov -= cov_offset;

    /* Function Body */
    krank = npq1;
    for (i__ = 1; i__ <= npq1; ++i__) {
	ss = s[i__];
	for (j = 1; j <= npq1; ++j) {
	    if (ss < 1.) {
		if (fabs(u[i__ + j * u_dim1]) > ss * machfd_.fltmax) {
		    krank = i__ - 1;
		    hessfd_.kcov = 1;
		    goto L100;
		}
	    }
	}
    }
L100:
    for (k = 1; k <= npq1; ++k) {
	F77_CALL(dcopy)(&k, &c_0d, &c__0, &cov[k * cov_dim1 + 1], &c__1);
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
    for (k = 1; k <= krank; ++k) {
	ss = -1. / s[k];
	for (j = 1; j <= npq1; ++j) {
	    d__1 = ss * u[j + k * u_dim1];
	    F77_CALL(daxpy)(&j, &d__1, &v[k * v_dim1 + 1], &c__1,
			    &cov[j * cov_dim1 + 1], &c__1);
	}
    }
    return 0;
} /* invsvd_

 ******************************************************************************
 *****************************************************************************/

void hesspq_(double *qp, double *a, double *ajac,
	     int *lajac, double *h__, int *lh, double *aij, double *g)
{
/*   double precision	qp(npq), a(nm), ajac(nm,npq)
     double precision	H(lH,npq1), aij(nm), g(npq)

 analytic Hessian with respect to p and q variables
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 ----------------------------------------------------------------------------*/

    /* System generated locals */
    int ajac_dim1, ajac_offset, h_dim1, h_offset, i__1, i__2, i__3;

    /* Local variables */
    int i__, j, k, l, km;
    double s, t, u, fac;

    /* Parameter adjustments */
    --qp;
    --a;
    ajac_dim1 = *lajac; ajac_offset = 1 + ajac_dim1; ajac -= ajac_offset;
    h_dim1 = *lh;	h_offset    = 1 + h_dim1;     h__ -= h_offset;
    --aij;
    --g;

    fac = 1. / (filtfd_.wnv * (double) (Dims.nm - 1));
    if (Dims.nq != 0 && Dims.np != 0) {
	for (k = 1; k <= Dims.npq; ++k) {
	    g[k] = F77_CALL(ddot)(&Dims.nm, &a[1], &c__1,
				  &ajac[k * ajac_dim1 + 1], &c__1);
	}
	i__1 = Dims.np;
	i__2 = Dims.nq;
	i__3 = Dims.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    u = g[Dims.nq + i__];
	    int i__a = (Dims.nq + i__)* ajac_dim1;
	    for (j = 1; j <= i__2; ++j) {
		u *= g[j];
		for (k = Dims.maxpq1; k <= i__3; ++k) {
		    km = k - Dims.maxpq;
		    t = 0.;
		    for (l = 1; l < km && l <= i__2; ++l)
			t += qp[l] * aij[km - l];

		    if (km > j)
			aij[km] = ajac[km - j + i__a] + t;
		    else
			aij[km] = t;
		}
		s = F77_CALL(ddot)(&Dims.nm, &ajac[i__a + 1], &c__1,
				   &ajac[j * ajac_dim1 + 1], &c__1);
		t = F77_CALL(ddot)(&Dims.nm, &a[1], &c__1, &aij[1], &c__1);
		h__[i__ + 1 + (Dims.np + j + 1) * h_dim1] =
		    -((double) Dims.n) * (s + t - fac * 2. * u) * fac;
	    }
	}
    }
    if (Dims.nq != 0) {
	i__1 = Dims.nq;
	i__3 = Dims.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    u = g[i__];
	    int i__a = i__ * ajac_dim1;
	    for (j = i__; j <= i__1; ++j) {
		u *= g[j];
		int j__a = j * ajac_dim1;
		for (k = Dims.maxpq1; k <= i__3; ++k) {
		    km = k - Dims.maxpq;
		    t = 0.;
		    for (l = 1; l < km && l <= i__1; ++l)
			t += qp[l] * aij[km - l];

		    s = 0.;
		    if (km > i__) s += ajac[km - i__ + j__a];
		    if (km > j)	  s += ajac[km - j + i__a];

		    aij[km] = s + t;
		}
		s = F77_CALL(ddot)(&Dims.nm, &ajac[i__a + 1], &c__1,
				   &ajac[j__a + 1], &c__1);
		t = F77_CALL(ddot)(&Dims.nm, &a[1], &c__1, &aij[1], &c__1);
		h__[Dims.np + i__ + 1 + (Dims.np + j + 1) * h_dim1] =
			-((double) Dims.n) * (s + t - fac * 2. * u) * fac;
	    }
	}
    }
    if (Dims.np != 0) {
	i__1 = Dims.np;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    u = g[Dims.nq + i__];
	    for (j = i__; j <= i__1; ++j) {
		u = g[Dims.nq + j] * u;
/*	      do k = maxpq1, n */
/*		km  =  k - maxpq */
/*		t  = zero */
/*		if (nq .ne. 0) then */
/*		 do l = 1, nq */
/*		    if (km .le. l) goto 303 */
/*		    t  = t + qp(l)*aij(km-l) */
/*		 end do */
/*		end if */
/* 303		continue */
/*		aij(km) = t */
/*	      end do */
		s = F77_CALL(ddot)(&Dims.nm, &ajac[(Dims.nq+i__)*ajac_dim1 + 1],
				   &c__1, &ajac[(Dims.nq + j) * ajac_dim1 + 1],
				   &c__1);
/*	      t = ddot( nm, a		  , 1, aij	     , 1) */
/*	      H(i+1,j+1) = -dble(n)*((s + t) - two*fac*u)*fac */
		h__[i__ + 1 + (j + 1) * h_dim1] =
		    -((double) Dims.n) * (s - fac * 2. * u) * fac;
	    }
	}
    }
    return;
} /* hesspq_ */


/******************************************************************************
 ******************************************************************************
 Subroutine */ int
hesdpq_(double *x, double *d__, double *hh, double *hd, double *w)
{
/*     float		 x(n)
     double precision	 d, hh, hd(npq1), w(*)
*/
    /* System generated locals */
    double d__1;

    /* Local variables */
    double fa, fb, slogvk;

/* copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 -----------------------------------------------------------------------------
     Parameter adjustments */
    --w;
    --hd;

    /* Function Body */
    if (*hh <= 0.) {
	*hh = (fabs(filtfd_.cllf) + 1.) * mauxfd_.epspt5;
    }
    if(*hh > 0.1) *hh = 0.1;
    if (*d__ - *hh > 0.) {
	fdfilt(x, *d__ - *hh, &w[wfilfd_.ly], &slogvk,
	       &w[wfilfd_.lamk], &w[wfilfd_.lak], &w[wfilfd_.lvk],
	       &w[wfilfd_.lphi], &w[wfilfd_.lpi]);
	if (Dims.npq != 0) {
	    ajqp_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac],
		  &Dims.nm, &c__1, &w[wfilfd_.ly]);
	    ajqp_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac],
		  &Dims.nm, &c__2, &w[wfilfd_.ly]);
	    gradpq(&w[woptfd_.lwa1], &w[woptfd_.la], &w[woptfd_.lajac],Dims.nm);
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[woptfd_.la], &c__1,
					  &w[woptfd_.la], &c__1);
	    d__1 = 1. / filtfd_.wnv;
	    F77_CALL(dscal)(&Dims.npq, &d__1, &w[woptfd_.lwa1], &c__1);
	    filtfd_.wnv /= (double) (Dims.nm - 1);
	} else {
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[wfilfd_.ly], &c__1,
					  &w[wfilfd_.ly], &c__1) /
		(double) (Dims.nm - 1);
	}
	fa = -((double) Dims.n * (log(filtfd_.wnv) + 2.8378) + slogvk) / 2.;
	if (*d__ + *hh < .5) {
	    fdfilt(x, *d__ + *hh, &w[wfilfd_.ly], &slogvk,
		   &w[wfilfd_.lamk], &w[wfilfd_.lak], &w[wfilfd_.lvk],
		   &w[wfilfd_.lphi], &w[wfilfd_.lpi]);
	    if (Dims.npq != 0) {
		ajqp_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac],
		      &Dims.nm, &c__1, &w[wfilfd_.ly]);
		ajqp_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac],
		      &Dims.nm, &c__2, &w[wfilfd_.ly]);
		gradpq(&w[woptfd_.lwa2], &w[woptfd_.la], &w[woptfd_.lajac],
		       Dims.nm);
		filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[woptfd_.la], &c__1,
					      &w[woptfd_.la], &c__1);
		d__1 = 1. / filtfd_.wnv;
		F77_CALL(dscal)(&Dims.npq, &d__1, &w[woptfd_.lwa2], &c__1);
		filtfd_.wnv /= (double) (Dims.nm - 1);
	    } else {
		filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[wfilfd_.ly], &c__1,
					      &w[wfilfd_.ly], &c__1) /
		    (double) (Dims.nm - 1);
	    }
	    fb = -((double) Dims.n * (log(filtfd_.wnv) + 2.8378) + slogvk)/ 2.;
	    hd[1] = (fa + fb - filtfd_.cllf * 2.) / (*hh * *hh);
	}
	else {
	    fdfilt(x, *d__ - *hh * 2., &w[wfilfd_.ly], &slogvk,
		   &w[wfilfd_.lamk], &w[wfilfd_.lak], &w[wfilfd_.lvk],
		   &w[wfilfd_.lphi], &w[wfilfd_.lpi]);
	    if (Dims.npq != 0) {
		ajqp_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac],
		      &Dims.nm, &c__1, &w[wfilfd_.ly]);
		ajqp_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac],
		      &Dims.nm, &c__2, &w[wfilfd_.ly]);
		gradpq(&w[woptfd_.lwa2], &w[woptfd_.la], &w[woptfd_.lajac],
		       Dims.nm);
		filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[woptfd_.la], &c__1,
					      &w[woptfd_.la], &c__1);
		d__1 = 1. / filtfd_.wnv;
		F77_CALL(dscal)(&Dims.npq, &d__1, &w[woptfd_.lwa2], &c__1);
		filtfd_.wnv /= (double) (Dims.nm - 1);
	    } else {
		filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[wfilfd_.ly], &c__1,
					      &w[wfilfd_.ly], &c__1) /
		    (double) (Dims.nm - 1);
	    }
	    fb = -((double)Dims.n * (log(filtfd_.wnv) + 2.8378) + slogvk) / 2.;
	    hd[1] = (filtfd_.cllf + fb - fa * 2.) / (*hh * 2. * *hh);
	}
    } else {
	fdfilt(x, *d__ + *hh, &w[wfilfd_.ly], &slogvk,
	       &w[wfilfd_.lamk], &w[wfilfd_.lak], &w[wfilfd_.lvk],
	       &w[wfilfd_.lphi], &w[wfilfd_.lpi]);
	if (Dims.npq != 0) {
	    ajqp_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac],
		  &Dims.nm, &c__1, &w[wfilfd_.ly]);
	    ajqp_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac],
		  &Dims.nm, &c__2, &w[wfilfd_.ly]);
	    gradpq(&w[woptfd_.lwa1], &w[woptfd_.la], &w[woptfd_.lajac],Dims.nm);
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[woptfd_.la], &c__1,
					  &w[woptfd_.la], &c__1);
	    d__1 = 1. / filtfd_.wnv;
	    F77_CALL(dscal)(&Dims.npq, &d__1, &w[woptfd_.lwa1], &c__1);
	    filtfd_.wnv /= (double) (Dims.nm - 1);
	} else {
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[wfilfd_.ly], &c__1,
					  &w[wfilfd_.ly], &c__1) /
		(double) (Dims.nm - 1);
	}
	fa = -((double) Dims.n * (log(filtfd_.wnv) + 2.8378) + slogvk) / 2.;
	fdfilt(x, *d__ + *hh * 2., &w[wfilfd_.ly], &slogvk,
	       &w[wfilfd_.lamk], &w[wfilfd_.lak], &w[wfilfd_.lvk],
	       &w[wfilfd_.lphi], &w[wfilfd_.lpi]);
	if (Dims.npq != 0) {
	    ajqp_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac],
		  &Dims.nm, &c__1, &w[wfilfd_.ly]);
	    ajqp_(&w[woptfd_.lqp], &w[woptfd_.la], &w[woptfd_.lajac],
		  &Dims.nm, &c__2, &w[wfilfd_.ly]);
	    gradpq(&w[woptfd_.lwa1], &w[woptfd_.la], &w[woptfd_.lajac],Dims.nm);
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[woptfd_.la], &c__1,
					  &w[woptfd_.la], &c__1);
	    d__1 = 1. / filtfd_.wnv;
	    F77_CALL(dscal)(&Dims.npq, &d__1, &w[woptfd_.lwa1], &c__1);
	    filtfd_.wnv /= (double) (Dims.nm - 1);
	} else {
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[wfilfd_.ly], &c__1,
					  &w[wfilfd_.ly], &c__1) /
		(double) (Dims.nm - 1);
	}
	fb = -((double) Dims.n * (log(filtfd_.wnv) + 2.8378) + slogvk) / 2.;
	hd[1] = (filtfd_.cllf + fb - fa * 2.) / (*hh * 2. * *hh);
    }
    if (Dims.npq == 0) {
	return 0;
    }
    F77_CALL(daxpy)(&Dims.npq, &c_m1, &w[woptfd_.lwa2], &c__1,
		    &w[woptfd_.lwa1], &c__1);
    d__1 = (double) Dims.n / (*hh * 2.);
    F77_CALL(dscal)(&Dims.npq, &d__1, &w[woptfd_.lwa1], &c__1);
    F77_CALL(dcopy)(&Dims.npq, &w[woptfd_.lwa1], &c__1, &hd[2], &c__1);
    return 0;
} /* hesdpq_ */

/******************************************************************************
 *****************************************************************************/

void gradpq(double *g, double a[], double ajac[], int l_ajac)
{
/*     double precision	g(npq), a(nm), ajac(nm,npq)
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 -----------------------------------------------------------------------------*/

    int i, j;

    for (i = 0; i < Dims.np; ++i)
	g[i] = F77_CALL(ddot)(&Dims.nm, a, &c__1,
			      &ajac[(Dims.nq + i) * l_ajac], &c__1);

    for (j = 0; j < Dims.nq; ++j)
	g[Dims.np + j] = F77_CALL(ddot)(&Dims.nm, a, &c__1,
					&ajac[j * l_ajac], &c__1);
    return;
} /* gradpq */

