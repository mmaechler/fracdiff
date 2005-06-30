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
void hesdpq(double *, double, double *, double *, double *);
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
/* double precision	H(lH, pq1)
*/

/*  copyright 1991 Department of Statistics, University of Washington
  written by Chris Fraley
 -----------------------------------------------------------------------------
     Parameter adjustments */
    --w;

    /* Function Body */
    hesspq_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac], &Dims.nm,
	    h__, lh, &w[w_opt.lwa4], &w[w_opt.lwa1]);
/*     call dcopy( pq1, zero, 0, H(1,1), lH) */
/*     call dcopy( pq , zero, 0, H(2,1), 1) */
    return;
} /* fdhpq */

/*******************************************************************************
 ****************************************************************************** */

void fdcov(double *x, double *d__, double *hh,
	   double *hd, double *cov, int *lcov, double *cor,
	   int *lcor, double *se, double *w, int *info)
{
    /* System generated locals */
    int cov_dim1, cov_offset, cor_dim1, cor_offset, pq1, i2;
    double d__1;

    /* Local variables */
    int i, j, k, le, ls, lu, lv, lwork;
    double temp;

/*     float		   x(n)
     double precision	d, hh, hd(pq1), cov(lcov,pq1),
    *			cor(lcor,pq1), se(pq1)
  copyright 1991 Department of Statistics, University of Washington
  written by Chris Fraley
  ----------------------------------------------------------------------------*/

    /* Parameter adjustments */
    cov_dim1 = *lcov;    cov_offset = 1 + cov_dim1;    cov -= cov_offset;
    cor_dim1 = *lcor;    cor_offset = 1 + cor_dim1;    cor -= cor_offset;
    --se;
    --w;

    /* Function Body */
    pq1 = Dims.pq1;

    hesdpq(x, *d__, hh, hd, &w[1]);
    F77_CALL(dcopy)(&pq1, hd, &c__1, &cov[cov_offset], lcov);

    gammfd_.igamma = 0;
    gammfd_.jgamma = 0;
    hessfd_.ksvd = 0;
    hessfd_.kcov = 0;
    hessfd_.kcor = 0;
    *info = 0;
    temp = 1.;

    for (i = 1; i <= pq1; ++i) {
	for (j = i + 1; j <= pq1; ++j) {
	    cov[j + i * cov_dim1] = cov[i + j * cov_dim1];
	}
    }
    ls = w_fil.ly;
    lu = ls + pq1 + 1;
    lv = lu + pq1 * pq1;
    le = lv + pq1 * pq1;
    lwork = le + pq1;
/*	lfree = lwork + pq1 */
    F77_CALL(dsvdc)(&cov[cov_offset], lcov, &pq1, &pq1, &w[ls],
		    &w[le], &w[lu], &pq1, &w[lv], &pq1, &w[lwork],
		    &c__11, info);
    if (*info != 0) {
	F77_CALL(dcopy)(&pq1, &c_0d, &c__0, &se[1], &c__1);
	for (j = 1; j <= pq1; ++j) {
	    F77_CALL(dcopy)(&pq1, &c_0d, &c__0,
			    &cov[j * cov_dim1 + 1], & c__1);
	}
	hessfd_.ksvd = 1;
	*info = 3;
	return;
    }
    invsvd_(&w[ls], &w[lu], &pq1, &w[lv], &pq1, &cov[cov_offset], lcov);
    for (i = 1; i <= pq1; ++i) {
	for (j = i + 1; j <= pq1; ++j) {
	    cov[j + i * cov_dim1] = cov[i + j * cov_dim1];
	}
    }
    temp = 1.;
    for (j = 1; j <= pq1; ++j) {
	if (cov[j + j * cov_dim1] > 0.) {
	    se[j] = sqrt(cov[j + j * cov_dim1]);
	} else {
	    temp = fmin2(temp, cov[j + j * cov_dim1]);
	    se[j] = 0.;
	}
    }
    if (temp == 1.) {
	for (k = 1; k <= pq1; ++k) {
	    F77_CALL(dcopy)(&k, &cov[k * cov_dim1 + 1], &c__1,
				&cor[k * cor_dim1 + 1], &c__1);
	}
	for (i = 1; i <= pq1; ++i) {
	    i2 = pq1 - i + 1;
	    d__1 = 1. / se[i];
	    F77_CALL(dscal)(&i2, &d__1, &cor[i + i * cor_dim1], lcor);
	}
	for (j = 1; j <= pq1; ++j) {
	    d__1 = 1. / se[j];
	    F77_CALL(dscal)(&j, &d__1, &cor[j * cor_dim1 + 1], &c__1);
	}
    } else {
	hessfd_.kcor = 1;
	for (j = 1; j <= pq1; ++j) {
	    F77_CALL(dcopy)(&pq1, &c_0d, &c__0,
			    &cor[j * cor_dim1 + 1], &c__1);
	}
    }
    for (i = 1; i <= pq1; ++i)
	for (j = i + 1; j <= pq1; ++j)
	    cor[j + i * cor_dim1] = cor[i + j * cor_dim1];

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
/*     double precision   s(pq1), u(lu,pq1), v(lv,pq1), cov(lcov,pq1)
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 ---------------------------------------------------------------------------*/

    /* System generated locals */
    int u_dim1, u_offset, v_dim1, v_offset, cov_dim1, cov_offset;
    double d__1;

    /* Local variables */
    int i__, j, k, krank, pq1 = Dims.pq1;
    double ss;

    /* Parameter adjustments */
    --s;
    u_dim1 = *lu;	u_offset = 1 + u_dim1;		  u -= u_offset;
    v_dim1 = *lv;	v_offset = 1 + v_dim1;		  v -= v_offset;
    cov_dim1 = *lcov; cov_offset = 1 + cov_dim1;	cov -= cov_offset;

    /* Function Body */
    krank = pq1;
    for (i__ = 1; i__ <= pq1; ++i__) {
	ss = s[i__];
	for (j = 1; j <= pq1; ++j) {
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
    for (k = 1; k <= pq1; ++k) {
	F77_CALL(dcopy)(&k, &c_0d, &c__0, &cov[k * cov_dim1 + 1], &c__1);
    }
    if (krank == 0) {
	return 0;
    }
/*      do k = 1, pq1 */
/*        do i = 1, pq1 */
/*          do j = i, pq1 */
/*            H(i,j) =  H(i,j) + s(k)*u(i,k)*v(j,k) */
/*          end do */
/*        end do */
/*      end do */
/*      do k = 1, pq1 */
/*        ss = s(k) */
/*        do j = 1, pq1 */
/*          call daxpy( j, ss*v(j,k), u(1,k), 1, H(1,j), 1) */
/*        end do */
/*      end do */
    for (k = 1; k <= krank; ++k) {
	ss = -1. / s[k];
	for (j = 1; j <= pq1; ++j) {
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
/*   double precision	qp(pq), a(nm), ajac(nm,pq)
     double precision	H(lH,pq1), aij(nm), g(pq)

 analytic Hessian with respect to p and q variables
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 ----------------------------------------------------------------------------*/

    /* System generated locals */
    int ajac_dim1 = *lajac, ajac_offset;
    int h_dim1 = *lh;
    int i1, i2, i3;

    /* Local variables */
    int i, j, k, l, km;
    double s, t, u, fac;

    /* Parameter adjustments */
    --qp;
    ajac_offset = 1 + ajac_dim1; ajac -= ajac_offset;
    --aij;
    --g;

    fac = 1. / (filtfd_.wnv * (double) (Dims.nm - 1));
    if (Dims.q != 0 && Dims.p != 0) {
	for (k = 1; k <= Dims.pq; ++k) {
	    g[k] = F77_CALL(ddot)(&Dims.nm, a, &c__1,
				  &ajac[k * ajac_dim1 + 1], &c__1);
	}
	i1 = Dims.p;
	i2 = Dims.q;
	i3 = Dims.n;
	for (i = 1; i <= i1; ++i) {
	    u = g[Dims.q + i];
	    int i_aj = (Dims.q + i)* ajac_dim1;
	    for (j = 1; j <= i2; ++j) {
		u *= g[j];
		for (k = Dims.maxpq1; k <= i3; ++k) {
		    km = k - Dims.maxpq;
		    t = 0.;
		    for (l = 1; l < km && l <= i2; ++l)
			t += qp[l] * aij[km - l];

		    if (km > j)
			aij[km] = ajac[km - j + i_aj] + t;
		    else
			aij[km] = t;
		}
		s = F77_CALL(ddot)(&Dims.nm, &ajac[i_aj + 1], &c__1,
				   &ajac[j * ajac_dim1 + 1], &c__1);
		t = F77_CALL(ddot)(&Dims.nm, a, &c__1, &aij[1], &c__1);
		h__[i + (Dims.p + j) * h_dim1] =
		    - Dims.n * (s + t - 2 * fac * u) * fac;
	    }
	}
    }
    if (Dims.q != 0) {
	i1 = Dims.q;
	i3 = Dims.n;
	for (i = 1; i <= i1; ++i) {
	    u = g[i];
	    int i_aj = i * ajac_dim1;
	    for (j = i; j <= i1; ++j) {
		u *= g[j];
		int j_aj = j * ajac_dim1;
		for (k = Dims.maxpq1; k <= i3; ++k) {
		    km = k - Dims.maxpq;
		    t = 0.;
		    for (l = 1; l < km && l <= i1; ++l)
			t += qp[l] * aij[km - l];

		    s = 0.;
		    if (km > i) s += ajac[km - i + j_aj];
		    if (km > j) s += ajac[km - j + i_aj];

		    aij[km] = s + t;
		}
		s = F77_CALL(ddot)(&Dims.nm, &ajac[i_aj + 1], &c__1,
				   &ajac[j_aj + 1], &c__1);
		t = F77_CALL(ddot)(&Dims.nm, a, &c__1, &aij[1], &c__1);
		h__[Dims.p + i + (Dims.p + j) * h_dim1] =
			-Dims.n * (s + t - 2 * fac * u) * fac;
	    }
	}
    }
    if (Dims.p != 0) {
	i1 = Dims.p;
	for (i = 1; i <= i1; ++i) {
	    u = g[Dims.q + i];
	    for (j = i; j <= i1; ++j) {
		u = g[Dims.q + j] * u;
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

		/* t = ddot( nm, a , 1, aij , 1) */
		s = F77_CALL(ddot)(&Dims.nm, &ajac[(Dims.q+i)*ajac_dim1 + 1],
				   &c__1, &ajac[(Dims.q + j) * ajac_dim1 + 1],
				   &c__1);

		/* H(i+1,j+1) = -dble(n)*((s + t) - two*fac*u)*fac */
		h__[i + (j) * h_dim1] = - Dims.n * (s - 2 * fac * u) * fac;
	    }
	}
    }
    return;
} /* hesspq_ */


/******************************************************************************
 *****************************************************************************/
void
hesdpq(double *x, double d_, double *hh, double *hd, double *w)
{
/*     float		 x(n)
       double precision	 d, hh, hd(pq1), w(*)

 * copyright 1991 Department of Statistics, University of Washington
   written by Chris Fraley
   ---------------------------------------------------------------------------*/

    /* System generated locals */
    double d__1;
    /* Local variables */
    double fa, fb, slogvk;

    /* Parameter adjustments */
    --w;

    /* Function Body */
    if (*hh <= 0.) {
	*hh = (fabs(filtfd_.cllf) + 1.) * mauxfd_.epspt5;
    }
    if(*hh > 0.1) *hh = 0.1;
    if (d_ - *hh > 0.) {
	fdfilt(x, d_ - *hh, &w[w_fil.ly], &slogvk,
	       &w[w_fil.lamk], &w[w_fil.lak], &w[w_fil.lvk],
	       &w[w_fil.lphi], &w[w_fil.lpi]);
	if (Dims.pq != 0) {
	    ajqp_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac],
		  &Dims.nm, &c__1, &w[w_fil.ly]);
	    ajqp_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac],
		  &Dims.nm, &c__2, &w[w_fil.ly]);
	    gradpq(&w[w_opt.lwa1], &w[w_opt.la], &w[w_opt.lajac],Dims.nm);
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[w_opt.la], &c__1,
					 &w[w_opt.la], &c__1);
	    d__1 = 1. / filtfd_.wnv;
	    F77_CALL(dscal)(&Dims.pq, &d__1, &w[w_opt.lwa1], &c__1);
	    filtfd_.wnv /= (Dims.nm - 1);
	} else {
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[w_fil.ly], &c__1,
					 &w[w_fil.ly], &c__1) / (Dims.nm - 1);
	}
	fa = -(Dims.n * (log(filtfd_.wnv) + 2.8378) + slogvk) / 2.;
	if (d_ + *hh < .5) {
	    fdfilt(x, d_ + *hh, &w[w_fil.ly], &slogvk,
		   &w[w_fil.lamk], &w[w_fil.lak], &w[w_fil.lvk],
		   &w[w_fil.lphi], &w[w_fil.lpi]);
	    if (Dims.pq != 0) {
		ajqp_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac],
		      &Dims.nm, &c__1, &w[w_fil.ly]);
		ajqp_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac],
		      &Dims.nm, &c__2, &w[w_fil.ly]);
		gradpq(&w[w_opt.lwa2], &w[w_opt.la], &w[w_opt.lajac],
		       Dims.nm);
		filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[w_opt.la], &c__1,
					     &w[w_opt.la], &c__1);
		d__1 = 1. / filtfd_.wnv;
		F77_CALL(dscal)(&Dims.pq, &d__1, &w[w_opt.lwa2], &c__1);
		filtfd_.wnv /= (Dims.nm - 1);
	    } else {
		filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[w_fil.ly], &c__1,
					     &w[w_fil.ly], &c__1) / (Dims.nm - 1);
	    }
	    fb = -(Dims.n * (log(filtfd_.wnv) + 2.8378) + slogvk)/ 2.;
	    hd[0] = (fa + fb - filtfd_.cllf * 2.) / (*hh * *hh);
	}
	else {
	    fdfilt(x, d_ - *hh * 2., &w[w_fil.ly], &slogvk,
		   &w[w_fil.lamk], &w[w_fil.lak], &w[w_fil.lvk],
		   &w[w_fil.lphi], &w[w_fil.lpi]);
	    if (Dims.pq != 0) {
		ajqp_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac],
		      &Dims.nm, &c__1, &w[w_fil.ly]);
		ajqp_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac],
		      &Dims.nm, &c__2, &w[w_fil.ly]);
		gradpq(&w[w_opt.lwa2], &w[w_opt.la], &w[w_opt.lajac],
		       Dims.nm);
		filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[w_opt.la], &c__1,
					     &w[w_opt.la], &c__1);
		d__1 = 1. / filtfd_.wnv;
		F77_CALL(dscal)(&Dims.pq, &d__1, &w[w_opt.lwa2], &c__1);
		filtfd_.wnv /= (Dims.nm - 1);
	    } else {
		filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[w_fil.ly], &c__1,
					     &w[w_fil.ly], &c__1) / (Dims.nm - 1);
	    }
	    fb = -(Dims.n * (log(filtfd_.wnv) + 2.8378) + slogvk) / 2.;
	    hd[0] = (filtfd_.cllf + fb - fa * 2.) / (*hh * 2. * *hh);
	}
    }
    else { /* (d_ <= *hh ) : */

	fdfilt(x, d_ + *hh, &w[w_fil.ly], &slogvk,
	       &w[w_fil.lamk], &w[w_fil.lak], &w[w_fil.lvk],
	       &w[w_fil.lphi], &w[w_fil.lpi]);
	if (Dims.pq != 0) {
	    ajqp_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac],
		  &Dims.nm, &c__1, &w[w_fil.ly]);
	    ajqp_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac],
		  &Dims.nm, &c__2, &w[w_fil.ly]);
	    gradpq(&w[w_opt.lwa1], &w[w_opt.la], &w[w_opt.lajac],Dims.nm);
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[w_opt.la], &c__1,
					 &w[w_opt.la], &c__1);
	    d__1 = 1. / filtfd_.wnv;
	    F77_CALL(dscal)(&Dims.pq, &d__1, &w[w_opt.lwa1], &c__1);
	    filtfd_.wnv /= (Dims.nm - 1);
	} else {
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[w_fil.ly], &c__1,
					 &w[w_fil.ly], &c__1) / (Dims.nm - 1);
	}
	fa = -(Dims.n * (log(filtfd_.wnv) + 2.8378) + slogvk) / 2.;
	fdfilt(x, d_ + *hh * 2., &w[w_fil.ly], &slogvk,
	       &w[w_fil.lamk], &w[w_fil.lak], &w[w_fil.lvk],
	       &w[w_fil.lphi], &w[w_fil.lpi]);
	if (Dims.pq != 0) {
	    ajqp_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac],
		  &Dims.nm, &c__1, &w[w_fil.ly]);
	    ajqp_(&w[w_opt.lqp], &w[w_opt.la], &w[w_opt.lajac],
		  &Dims.nm, &c__2, &w[w_fil.ly]);
	    gradpq(&w[w_opt.lwa1], &w[w_opt.la], &w[w_opt.lajac],Dims.nm);
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[w_opt.la], &c__1,
					 &w[w_opt.la], &c__1);
	    d__1 = 1. / filtfd_.wnv;
	    F77_CALL(dscal)(&Dims.pq, &d__1, &w[w_opt.lwa1], &c__1);
	    filtfd_.wnv /= (Dims.nm - 1);
	} else {
	    filtfd_.wnv = F77_CALL(ddot)(&Dims.nm, &w[w_fil.ly], &c__1,
					 &w[w_fil.ly], &c__1) / (Dims.nm - 1);
	}
	fb = -(Dims.n * (log(filtfd_.wnv) + 2.8378) + slogvk) / 2.;
	hd[0] = (filtfd_.cllf + fb - fa * 2.) / (*hh * 2. * *hh);
    }
    if (Dims.pq == 0) {
	return;
    }
    F77_CALL(daxpy)(&Dims.pq, &c_m1, &w[w_opt.lwa2], &c__1,
		    &w[w_opt.lwa1], &c__1);
    d__1 = Dims.n / (*hh * 2.);
    F77_CALL(dscal)(&Dims.pq, &d__1, &w[w_opt.lwa1], &c__1);
    F77_CALL(dcopy)(&Dims.pq, &w[w_opt.lwa1], &c__1, &hd[+1], &c__1);
    return;
} /* hesdpq */

/******************************************************************************
 *****************************************************************************/

void gradpq(double *g, double a[], double ajac[], int l_ajac)
{
/*     double precision	g(pq), a(nm), ajac(nm,pq)
 copyright 1991 Department of Statistics, University of Washington
 written by Chris Fraley
 -----------------------------------------------------------------------------*/

    int i, j;

    for (i = 0; i < Dims.p; ++i)
	g[i] = F77_CALL(ddot)(&Dims.nm, a, &c__1,
			      &ajac[(Dims.q + i) * l_ajac], &c__1);

    for (j = 0; j < Dims.q; ++j)
	g[Dims.p + j] = F77_CALL(ddot)(&Dims.nm, a, &c__1,
				       &ajac[j * l_ajac], &c__1);
    return;
} /* gradpq */

