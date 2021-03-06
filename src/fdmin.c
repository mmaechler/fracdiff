/*-*- mode: C; kept-old-versions: 12;  kept-new-versions: 20; -*-
 *
 * fdmin.f -- translated by f2c (version 20031025).
 * and produced by f2c-clean,v 1.10 2002/03/28 16:37:27 maechler
 *
 * and manually pretty edited by Martin Maechler, 2004-10-01, ff.
 */

#include <Rmath.h>
// for warning():
#include <R.h>

#include "fracdiff.h"


/* Common Block Declarations --- included as "extern" */
#define FD_EXTERNAL extern

#include "mach_comm.h"
#include "maux_comm.h"
// #include "tols_comm.h"


/* Constant (used to pass pointer) */
static /*logical*/int c_true = (1);

static void qrfac(int *, int *, double *, int *,
		  /*logical*/int *, int *, int *,
		  double *, double *, double *);
static void qrsolv(int, double *, int *, int *,
		   double *, double *, double *, double *, double *);

/* --------- EXPORTS ------------------- */

static double enorm(int, double *);

static double
lmpar(int, double *, int *, int *, double *, double *, double *,
      double,
      double *, double *, double *, double *);
/* and this : */

double lmder1(S_fp fcn, int m, int n,
	      double *x, double *fvec, double *fjac, int ldfjac,
	      double ftol, double xtol, double gtol, int maxfev, double *diag,
	      int mode, double factor,
	      int *info, int *nfev, int *njev,
	      int *ipvt, double *qtf,
	      double *wa1, double *wa2, double *wa3, double *wa4, double *y)
{
    // THE return value (since 2011-08-08):
    double fd_min_fnorm = -99.; // Wall
/*
       subroutine lmder

       the purpose of lmder is to minimize the sum of the squares of
       m nonlinear functions in n variables by a modification of
       the levenberg-marquardt algorithm. the user must provide a
       subroutine which calculates the functions and the jacobian.

       the subroutine statement is
       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
       maxfev,diag,mode,factor,nprint,info,nfev,
       njev,ipvt,qtf,wa1,wa2,wa3,wa4)

       where

       fcn is the name of the user-supplied subroutine which
       calculates the functions and the jacobian. fcn must
       be declared in an external statement in the user
       calling program, and should be written as follows.

       subroutine fcn(m, n, x,fvec,fjac, ldfjac,iflag)
       int m,n, ldfjac,iflag
       double precision x(n), fvec(m), fjac(ldfjac,n)
       ----------
       if iflag = 1 calculate the functions at x and
       return this vector in fvec. do not alter fjac.
       if iflag = 2 calculate the jacobian at x and
       return this matrix in fjac. do not alter fvec.
       ----------
       return
       end

       the value of iflag should not be changed by fcn unless
       the user wants to terminate execution of lmder.
       in this case set iflag to a negative int.

       m is a positive int input variable set to the number
       of functions.

       n is a positive int input variable set to the number
       of variables. n must not exceed m.

       x is an array of length n. on input x must contain
       an initial estimate of the solution vector. on output x
       contains the final estimate of the solution vector.

       fvec is an output array of length m which contains
       the functions evaluated at the output x.

       fjac is an output m by n array. the upper n by n submatrix
       of fjac contains an upper triangular matrix r with
       diagonal elements of nonincreasing magnitude such that

       t     t		  t
       p *(jac *jac)*p = r *r,

       where p is a permutation matrix and jac is the final
       calculated jacobian. column j of p is column ipvt(j)
       (see below) of the identity matrix. the lower trapezoidal
       part of fjac contains information generated during
       the computation of r.

       ldfjac is a positive int input variable not less than m
       which specifies the leading dimension of the array fjac.

       ftol is a nonnegative input variable. termination
       occurs when both the actual and predicted relative
       reductions in the sum of squares are at most ftol.
       therefore, ftol measures the relative error desired
       in the sum of squares.

       xtol is a nonnegative input variable. termination
       occurs when the relative error between two consecutive
       iterates is at most xtol. therefore, xtol measures the
       relative error desired in the approximate solution.

       gtol is a nonnegative input variable. termination
       occurs when the cosine of the angle between fvec and
       any column of the jacobian is at most gtol in absolute
       value. therefore, gtol measures the orthogonality
       desired between the function vector and the columns
       of the jacobian.

       maxfev is a positive int input variable. termination
       occurs when the number of calls to fcn with iflag = 1
       has reached maxfev.

       diag is an array of length n. if mode = 1 (see
       below), diag is internally set. if mode = 2, diag
       must contain positive entries that serve as
       multiplicative scale factors for the variables.

       mode is an int input variable. if mode = 1, the
       variables will be scaled internally. if mode = 2,
       the scaling is specified by the input diag. other
       values of mode are equivalent to mode = 1.

       factor is a positive input variable used in determining the
       initial step bound. this bound is set to the product of
       factor and the euclidean norm of diag*x if nonzero, or else
       to factor itself. in most cases factor should lie in the
       interval (.1,100.).100. is a generally recommended value.

       nprint is an int input variable that enables controlled
       printing of iterates if it is positive. in this case,
       fcn is called with iflag = 0 at the beginning of the first
       iteration and every nprint iterations thereafter and
       immediately prior to return, with x, fvec, and fjac
       available for printing. fvec and fjac should not be
       altered. if nprint is not positive, no special calls
       of fcn with iflag = 0 are made.

       info is an int output variable. if the user has
       terminated execution, info is set to the (negative)
       value of iflag. see description of fcn. otherwise,
       info is set as follows.

       info = 0  improper input parameters.

       info = 1  both actual and predicted relative reductions
       in the sum of squares are at most ftol.

       info = 2  relative error between two consecutive iterates
       is at most xtol.

       info = 3  conditions for info = 1 and info = 2 both hold.

       info = 4  the cosine of the angle between fvec and any
       column of the jacobian is at most gtol in
       absolute value.

       info = 5  number of calls to fcn with iflag = 1 has
       reached maxfev.

       info = 6  ftol is too small. no further reduction in
       the sum of squares is possible.

       info = 7  xtol is too small. no further improvement in
       the approximate solution x is possible.

       info = 8  gtol is too small. fvec is orthogonal to the
       columns of the jacobian to machine precision.

       nfev is an int output variable set to the number of
       calls to fcn with iflag = 1.

       njev is an int output variable set to the number of
       calls to fcn with iflag = 2.

       ipvt is an int output array of length n. ipvt
       defines a permutation matrix p such that jac*p = q*r,
       where jac is the final calculated jacobian, q is
       orthogonal (not stored), and r is upper triangular
       with diagonal elements of nonincreasing magnitude.
       column j of p is column ipvt(j) of the identity matrix.

       qtf is an output array of length n which contains
       the first n elements of the vector (q transpose)*fvec.

       wa1, wa2, and wa3 are work arrays of length n.

       wa4 is a work array of length m.

       subprograms called

       user-supplied ...... fcn

       minpack-supplied ... dpmpar,enorm,lmpar,qrfac

       fortran-supplied ... fabs,dmax1,dmin1,dsqrt,mod

       argonne national laboratory. minpack project. march 1980.
       burton s. garbow, kenneth e. hillstrom, jorge j. more

       epsmch is the machine precision.
       epsmch = dpmpar(1)
*/

    /* Initialized data */
    static double p1 = .1;
    static double p5 = .5;
    static double p25 = .25;
    static double p75 = .75;
    static double p0001 = 1e-4;

    /* System generated locals */
    int fjac_offset;
    double d__1;

    /* Local variables */
    int i__, j, l, iter, iflag, nprint;
    double par, sum, temp, temp1, temp2,
	ratio, enorm_n, xnorm, fnorm1, actred, dirder, prered,
	T_gnorm, delta = 0.; /* Wall*/

    /* Parameter adjustments */
    --wa4;
    --fvec;
    --wa3;
    --wa2;
    --wa1;
    --qtf;
    --ipvt;
    --diag;
    --x;
    --y;
    fjac_offset = 1 + ldfjac;	fjac -= fjac_offset;

    /* Function Body */
    temp = 0.;
    nprint = 0;
    *info = 0;
    iflag = 0;
    *nfev = 0;
    *njev = 0;
/*     check the input parameters for errors. */

    if (n <= 0 || m < n || ldfjac < m || ftol < 0. || xtol < 0. ||
	gtol < 0. || maxfev <= 0 || factor <= 0.) {
	warning("lmder1(): invalid (scalar) input");
	goto L_end;
    }
    if (mode == 2) { /* check diag[] */
	for (j = 1; j <= n; ++j)
	    if (diag[j] <= 0.) goto L_end;
    }

/* evaluate the function at the starting point and calculate its norm. */

    iflag = 1;
    (*fcn)(&x[1], &fvec[1], &fjac[fjac_offset], ldfjac, iflag, &y[1]);
    *nfev = 1;
    if (iflag < 0) {
	warning("lmder1(): problem in function evaluation at starting point");
	goto L_end;
    }

    fd_min_fnorm = fmin2(enorm(m, &fvec[1]), mauxfd_1.bignum);

/*     initialize levenberg-marquardt parameter and iteration counter. */

    par = 0.;
    iter = 1;

/* ==== beginning of the outer loop. ==========================================*/
L30:

/*	  calculate the jacobian matrix. */

    iflag = 2;
    (*fcn)(&x[1], &fvec[1], &fjac[fjac_offset], ldfjac, iflag, &y[1]);
    ++(*njev);
    if (iflag < 0)
	goto L_end;

/*	  if requested, call fcn to enable printing of iterates. */
    if (nprint > 0) {
	iflag = 0;
	if ((iter - 1) % nprint == 0)
	    (*fcn)(&x[1], &fvec[1], &fjac[fjac_offset], ldfjac, iflag, &y[1]);

	if (iflag < 0)
	    goto L_end;
    }
    /* L40: */

/*	  compute the qr factorization of the jacobian. */

    qrfac(&m, &n, &fjac[fjac_offset], &ldfjac, &c_true, &ipvt[1], &n,
	  &wa1[1], &wa2[1], &wa3[1]);

    /* on the first iteration -- do a some initializations : */
    if (iter == 1) {
	/* if mode is 1, scale according
	   to the norms of the columns of the initial jacobian. */
	if (mode == 1) {
	    for (j = 1; j <= n; ++j)
		diag[j] = ((wa2[j] != 0.)? wa2[j] : 1.);
	}

	/* calculate the norm of the scaled x and
	   initialize the step bound delta. */

	for (j = 1; j <= n; ++j)
	    wa3[j] = diag[j] * x[j];

	xnorm = enorm(n, &wa3[1]);
	delta = factor * xnorm;
	if (delta == 0.) {
	    delta = factor;
	}
    }
    /* L80: */

/*	  form (q transpose)*fvec and store the first n components in qtf. */

    for (i__ = 1; i__ <= m; ++i__) {
	wa4[i__] = fvec[i__];
    }
    for (j = 1; j <= n; ++j) {
	if (fjac[j + j * ldfjac] != 0.) {
	    sum = 0.;
	    for (i__ = j; i__ <= m; ++i__)
		sum += fjac[i__ + j * ldfjac] * wa4[i__];

	    temp = -sum / fjac[j + j * ldfjac];
	    for (i__ = j; i__ <= m; ++i__)
		wa4[i__] += fjac[i__ + j * ldfjac] * temp;
	}
	/* L120: */
	fjac[j + j * ldfjac] = wa1[j];
	qtf[j] = wa4[j];
    }

/*	  compute the norm of the scaled gradient. */

    T_gnorm = 0.;
    if (fd_min_fnorm != 0.) {
	for (j = 1; j <= n; ++j) {
	    l = ipvt[j];
	    if (wa2[l] != 0.) {
		sum = 0.;
		for (i__ = 1; i__ <= j; ++i__)
		    sum += fjac[i__ + j * ldfjac] * (qtf[i__] / fd_min_fnorm);
		T_gnorm = fmax2(T_gnorm, fabs(sum / wa2[l]));
	    }
	}
    }
    /* L170: */

    /* test for convergence of the gradient norm. */

    if (T_gnorm <= gtol)	*info = 4;

    if (*info != 0)
	goto L_end;

    /* rescale if necessary. */

    if (mode == 1) {
	for (j = 1; j <= n; ++j)
	    diag[j] = fmax2(diag[j], wa2[j]);

    }
    /* L190: */

    do { // ------------- the inner loop. ------------------------------------

/*	     determine the levenberg-marquardt parameter. */

	par = lmpar(n, &fjac[fjac_offset], &ldfjac, &ipvt[1],
		    &diag[1], &qtf[1], &delta,
		    par,
		    &wa1[1], &wa2[1], &wa3[1], &wa4[1]);

/*	     store the direction p and x + p. calculate the norm of p. */

	for (j = 1; j <= n; ++j) {
	    wa1[j] = -wa1[j];
	    wa2[j] = x[j] + wa1[j];
	    wa3[j] = diag[j] * wa1[j];
	}
	enorm_n = enorm(n, &wa3[1]);

/*	     on the first iteration, adjust the initial step bound. */

	if (iter == 1) {
	    delta = fmin2(delta,enorm_n);
	}

/*	     evaluate the function at x + p and calculate its norm. */

	iflag = 1;
	(*fcn)(&wa2[1], &wa4[1], &fjac[fjac_offset], ldfjac, iflag, &y[1]);
	++(*nfev);
	if (iflag < 0)
	    goto L_end;

	fnorm1 = fmin2(enorm(m, &wa4[1]), mauxfd_1.bignum);

/*	     compute the scaled actual reduction. */

	actred = -1.;
	if (p1 * fnorm1 < fd_min_fnorm) {
	    d__1 = fnorm1 / fd_min_fnorm;
	    actred = 1. - d__1 * d__1;
	}
	/* actred = (fnorm*fnorm - fnorm1*fnorm1) */

	/* compute the scaled predicted reduction and
	   the scaled directional derivative. */

	for (j = 1; j <= n; ++j) {
	    wa3[j] = 0.;
	    l = ipvt[j];
	    temp = wa1[l];
	    for (i__ = 1; i__ <= j; ++i__) {
		wa3[i__] += fjac[i__ + j * ldfjac] * temp;
	    }
	}
	temp1 = enorm(n, &wa3[1]) / fd_min_fnorm;
	temp2 = sqrt(par) * enorm_n / fd_min_fnorm;

	prered = temp1 * temp1 + temp2 * temp2 / p5;
/*	     temp1  = enorm(n,wa3)
	     temp2  = (dsqrt(par)*enorm_n)
	     prered = (temp1**2 + 2.d0*temp2**2)
*/
	dirder = -(temp1 * temp1 + temp2 * temp2);

	/* compute the ratio of the actual to the predicted reduction. */

	if (prered != 0.)
	    ratio = actred / prered;
	else
	    ratio = 0.;

/*	     update the step bound. */

	if (ratio <= p25) {
	    if (actred >= 0.)
		temp = p5;
	    else /* (actred < 0.) */
		temp = p5 * dirder / (dirder + p5 * actred);

	    if (p1 * fnorm1 >= fd_min_fnorm || temp < p1)
		temp = p1;

	    delta = temp * fmin2(delta, enorm_n / p1);
	    par /= temp;
	}
	else { /* ratio > p25 */

	    if (par == 0. || ratio >= p75) {
		delta = enorm_n / p5;
		par = p5 * par;
	    }
	}
	/* L260: */

/*	     test for successful iteration. */

	if (ratio >= p0001) {
/*	     successful iteration. update x, fvec, and their norms. */

	    for (j = 1; j <= n; ++j) {
		x[j] = wa2[j];
		wa2[j] = diag[j] * x[j];
	    }
	    for (i__ = 1; i__ <= m; ++i__)
		fvec[i__] = wa4[i__];

	    xnorm = enorm(n, &wa2[1]);
	    fd_min_fnorm = fnorm1;
	    ++iter;
	}
 /* L290:	tests for convergence. */

	if((fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.) ||
	   (fd_min_fnorm <= ftol))
	    *info = 1;

	if (delta <= xtol) {
	    *info = 2;
	    if (fabs(actred) <= ftol &&
		prered	 <= ftol &&
		p5 * ratio <= 1.)
		*info = 3;
	}

	if (*info != 0)
	    goto L_end;


/* tests for termination and stringent tolerances. */

	if (*nfev >= maxfev)			*info = 5;

	if (fabs(actred) <= machfd_.epsmax &&
	    prered	 <= machfd_.epsmax &&
	    p5 * ratio <= 1.)			*info = 6;

	if (delta <= machfd_.epsmax)		*info = 7;

	if (T_gnorm <= machfd_.epsmax)		*info = 8;

	if (*info != 0)
	    goto L_end;

/*	     end of the inner loop. repeat if iteration unsuccessful. */
    } while (ratio < p0001);

/*	  end of the outer loop. */
    goto L30;

L_end: //  termination, either normal or user imposed.

    if (iflag < 0) {
	*info = iflag;
    }
    iflag = 0;
    if (nprint > 0) {
	(*fcn)(&x[1], &fvec[1], &fjac[fjac_offset], ldfjac, iflag, &y[1]);
    }
    return fd_min_fnorm;
} /* lmder1 */


double enorm(int n, double *x)
{
    /* Initialized data */

    static double rdwarf = 3.834e-20;
    static double rgiant = 1.304e19;

    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    static int i__;
    static double s1, s2, s3, xabs, x1max, x3max, agiant, floatn;

/*     **********

     function enorm

     given an n-vector x, this function calculates the
     euclidean norm of x.

     the euclidean norm is computed by accumulating the sum of
     squares in three different sums. the sums of squares for the
     small and large components are scaled so that no overflows
     occur. non-destructive underflows are permitted. underflows
     and overflows do not occur in the computation of the unscaled
     sum of squares for the intermediate components.
     the definitions of small, intermediate and large components
     depend on two constants, rdwarf and rgiant. the main
     restrictions on these constants are that rdwarf**2 not
     underflow and rgiant**2 not overflow. the constants
     given here are suitable for every known computer.

     the function statement is

       double precision function enorm(n,x)

     where

       n is a positive int input variable.

       x is an input array of length n.

     subprograms called

       fortran-supplied ... fabs,dsqrt

     argonne national laboratory. minpack project. march 1980.
     burton s. garbow, kenneth e. hillstrom, jorge j. more

     **********
     Parameter adjustments */
    --x;

    /* Function Body */
    ret_val = -1.;
    s1 = 0.;
    s2 = 0.;
    s3 = 0.;
    x1max = 0.;
    x3max = 0.;
    floatn = (double) (n);
    agiant = rgiant / floatn;

    for (i__ = 1; i__ <= n; ++i__) {
	xabs = fabs(x[i__]);
	if (xabs > rdwarf && xabs < agiant) {/* sum for intermediate components.*/

	    s2 += xabs * xabs;

	} else if (xabs > rdwarf) { /*		sum for large components. */

	    if (xabs <= x1max) {

		/* Computing 2nd power */
		d__1 = xabs / x1max;
		s1 += d__1 * d__1;

	    } else {
		/* Computing 2nd power */
		d__1 = x1max / xabs;
		s1 = 1. + s1 * (d__1 * d__1);
		x1max = xabs;
	    }

	} else { /*				sum for small components. */

	    if (xabs <= x3max) {

		if (xabs != 0.) {
		    /* Computing 2nd power */
		    d__1 = xabs / x3max;
		    s3 += d__1 * d__1;
		}

	    } else {
		/* Computing 2nd power */
		d__1 = x3max / xabs;
		s3 = 1. + s3 * (d__1 * d__1);
		x3max = xabs;
	    }
	}

    } /* for(i ) */

/*     calculation of norm. */

    if (s1 == 0.) {

	if (s2 == 0.) {
	    ret_val = x3max * sqrt(s3);
	}
	else {
	    if (s2 >= x3max)
		ret_val = sqrt(s2 * (1. + x3max / s2 * (x3max * s3)));
	    else /* (s2 < x3max) */
		ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
	}
    }
    else {
	ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
    }

    return ret_val;
} /* enorm */


static
void qrfac(int *m, int *n, double *a, int *lda,
	   /*logical*/int *pivot, int *ipvt, int *lipvt, double *rdiag,
	   double *acnorm, double *wa)
{
    /* Initialized data */

    static double p05 = .05;

    /* System generated locals */
    int a_dim1, a_offset;
    double d__1;

    /* Local variables */
    int i__, j, k, jp1, minmn;
    double sum, temp, ajnorm;

/*     **********

     subroutine qrfac

     this subroutine uses householder transformations with column
     pivoting (optional) to compute a qr factorization of the
     m by n matrix a. that is, qrfac determines an orthogonal
     matrix q, a permutation matrix p, and an upper trapezoidal
     matrix r with diagonal elements of nonincreasing magnitude,
     such that a*p = q*r. the householder transformation for
     column k, k = 1,2,...,min(m,n), is of the form

                           t
           i - (1/u(k))*u*u

     where u has zeros in the first k-1 positions. the form of
     this transformation and the method of pivoting first
     appeared in the corresponding linpack subroutine.

     the subroutine statement is

       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa)

     where

       m is a positive int input variable set to the number
         of rows of a.

       n is a positive int input variable set to the number
         of columns of a.

       a is an m by n array. on input a contains the matrix for
         which the qr factorization is to be computed. on output
         the strict upper trapezoidal part of a contains the strict
         upper trapezoidal part of r, and the lower trapezoidal
         part of a contains a factored form of q (the non-trivial
         elements of the u vectors described above).

       lda is a positive int input variable not less than m
         which specifies the leading dimension of the array a.

       pivot is a *logical* input variable. if pivot is set true,
         then column pivoting is enforced. if pivot is set false,
         then no column pivoting is done.

       ipvt is an int output array of length lipvt. ipvt
         defines the permutation matrix p such that a*p = q*r.
         column j of p is column ipvt(j) of the identity matrix.
         if pivot is false, ipvt is not referenced.

       lipvt is a positive int input variable. if pivot is false,
         then lipvt may be as small as 1. if pivot is true, then
         lipvt must be at least n.

       rdiag is an output array of length n which contains the
         diagonal elements of r.

       acnorm is an output array of length n which contains the
         norms of the corresponding columns of the input matrix a.
         if this information is not needed, then acnorm can coincide
         with rdiag.

       wa is a work array of length n. if pivot is false, then wa
         can coincide with rdiag.

     subprograms called

       minpack-supplied ... dpmpar,enorm

       fortran-supplied ... dmax1,dsqrt,min0

     argonne national laboratory. minpack project. march 1980.
     burton s. garbow, kenneth e. hillstrom, jorge j. more

     **********
     double precision dpmpar,enorm
     Parameter adjustments */
    --wa;
    --acnorm;
    --rdiag;
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;

    // compute the initial column norms and initialize several arrays.

    for (j = 1; j <= *n; ++j) {
	acnorm[j] = enorm(*m, &a[j * a_dim1 + 1]);
	rdiag[j] = acnorm[j];
	wa[j] = rdiag[j];
	if (*pivot) ipvt[j] = j;
    }

/*     reduce a to r with householder transformations. */

    minmn = imin2(*m,*n);
    for (j = 1; j <= minmn; ++j) {
	if (*pivot) { // bring the column of largest norm into the pivot position.
	    int kmax = j;
	    for (k = j; k <= *n; ++k) {
		if (rdiag[k] > rdiag[kmax])
		    kmax = k;
	    }
	    if (kmax != j) {
		// swap  a[,j]  and  a[,kmax] :
		for (i__ = 1; i__ <= *m; ++i__) {
		    double t = a[i__ + j * a_dim1];
		    a[i__ + j * a_dim1] = a[i__ + kmax * a_dim1];
		    a[i__ + kmax * a_dim1] = t;
		}
		rdiag[kmax] = rdiag[j];
		wa[kmax] = wa[j];
		k = ipvt[j]; ipvt[j] = ipvt[kmax]; ipvt[kmax] = k;
	    }
	}
// L40:

/*        compute the householder transformation to reduce the
        j-th column of a to a multiple of the j-th unit vector. */

	ajnorm = enorm(*m - j + 1, &a[j + j * a_dim1]);

	if (ajnorm == 0.) {
	    goto L100;
	}
	if (a[j + j * a_dim1] < 0.) {
	    ajnorm = -ajnorm;
	}
	for (i__ = j; i__ <= *m; ++i__) {
	    a[i__ + j * a_dim1] /= ajnorm;
	}
	a[j + j * a_dim1] += 1.;

/*        apply the transformation to the remaining columns
        and update the norms. */

	jp1 = j + 1;
	for (k = jp1; k <= *n; ++k) {
	    sum = 0.;
	    for (i__ = j; i__ <= *m; ++i__) {
		sum += a[i__ + j * a_dim1] * a[i__ + k * a_dim1];
	    }
	    temp = sum / a[j + j * a_dim1];
	    for (i__ = j; i__ <= *m; ++i__) {
		a[i__ + k * a_dim1] -= temp * a[i__ + j * a_dim1];
	    }
	    if (*pivot && rdiag[k] != 0.) {
		temp = a[j + k * a_dim1] / rdiag[k];
		rdiag[k] *= sqrt((fmax2(0., 1. - temp * temp)));
		/* Computing 2nd power */
		d__1 = rdiag[k] / wa[k];
		if (p05 * (d__1 * d__1) < machfd_.epsmax) {
		    rdiag[k] = enorm(*m - j, &a[jp1 + k * a_dim1]);
		    wa[k] = rdiag[k];
		}
	    }
	}
L100:
	rdiag[j] = -ajnorm;
    }
    return;
} /* qrfac */

double lmpar(int n, double *r__, int *ldr, int *ipvt,
	     double *diag, double *qtb, double *delta,
	     double par_init,
	     double *x, double *sdiag, double *wa1, double *wa2)
{
    double par = par_init;
    //     ---- the return_value
/*
     subroutine lmpar

     given an m by n matrix a, an n by n nonsingular diagonal
     matrix d, an m-vector b, and a positive number delta,
     the problem is to determine a value for the parameter
     par such that if x solves the system

           a*x = b ,     sqrt(par)*d*x = 0 ,

     in the least squares sense, and dxnorm is the euclidean
     norm of d*x, then either par is 0. and

           (dxnorm-delta) <= 0.1*delta ,

     or par is positive and

           abs(dxnorm-delta) <= 0.1*delta .

     this subroutine completes the solution of the problem
     if it is provided with the necessary information from the
     qr factorization, with column pivoting, of a. that is, if
     a*p = q*r, where p is a permutation matrix, q has orthogonal
     columns, and r is an upper triangular matrix with diagonal
     elements of nonincreasing magnitude, then lmpar expects
     the full upper triangle of r, the permutation matrix p,
     and the first n components of (q transpose)*b. on output
     lmpar also provides an upper triangular matrix s such that

            t   t                   t
           p *(a *a + par*d*d)*p = s *s .

     s is employed within lmpar and may be of separate interest.

     only a few iterations are generally needed for convergence
     of the algorithm. if, however, the limit of 10 iterations
     is reached, then the output par will contain the best
     value obtained so far.

     the subroutine statement is

       subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag,
                        wa1,wa2)

     where

       n is a positive int input variable set to the order of r.

       r is an n by n array. on input the full upper triangle
         must contain the full upper triangle of the matrix r.
         on output the full upper triangle is unaltered, and the
         strict lower triangle contains the strict upper triangle
         (transposed) of the upper triangular matrix s.

       ldr is a positive int input variable not less than n
         which specifies the leading dimension of the array r.

       ipvt is an int input array of length n which defines the
         permutation matrix p such that a*p = q*r. column j of p
         is column ipvt(j) of the identity matrix.

       diag is an input array of length n which must contain the
         diagonal elements of the matrix d.

       qtb is an input array of length n which must contain the first
         n elements of the vector (q transpose)*b.

       delta is a positive input variable which specifies an upper
         bound on the euclidean norm of d*x.

       par is a nonnegative variable. on input par contains an
         initial estimate of the levenberg-marquardt parameter.
         on output par contains the final estimate.

       x is an output array of length n which contains the least
         squares solution of the system a*x = b, sqrt(par)*d*x = 0,
         for the output par.

       sdiag is an output array of length n which contains the
         diagonal elements of the upper triangular matrix s.

       wa1 and wa2 are work arrays of length n.

     subprograms called

       minpack-supplied ... dpmpar,enorm,qrsolv

       fortran-supplied ... fabs,dmax1,dmin1,dsqrt

     argonne national laboratory. minpack project. march 1980.
     burton s. garbow, kenneth e. hillstrom, jorge j. more

***********/

    /* Initialized data */
    static double p1 = .1;
    static double p001 = .001;

    /* System generated locals */
    int r_dim1, r_offset;

    /* Local variables */
    int i__, j, k, l, jp1, iter, nsing;
    double fp, sum, parc, parl, temp, paru, dwarf, gnorm, dxnorm;

    /* Parameter adjustments */
    --wa2;
    --wa1;
    --sdiag;
    --x;
    --qtb;
    --diag;
    --ipvt;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;

    // dwarf is the smallest positive magnitude :
    dwarf = machfd_.fltmin;

/*     compute and store in x the gauss-newton direction. if the
     jacobian is rank-deficient, obtain a least squares solution. */

    nsing = n;
    for (j = 1; j <= n; ++j) {
	wa1[j] = qtb[j];
	if (r__[j + j * r_dim1] == 0. && nsing == n) {
	    nsing = j - 1;
	}
	if (nsing < n) {
	    wa1[j] = 0.;
	}
    }
    for (k = 1; k <= nsing; ++k) {
	j = nsing - k + 1;
	wa1[j] /= r__[j + j * r_dim1];
	temp = wa1[j];
	for (i__ = 1; i__ <= j-1; ++i__) {
	    wa1[i__] -= r__[i__ + j * r_dim1] * temp;
	}
    }
    for (j = 1; j <= n; ++j) {
	l = ipvt[j];
	x[l] = wa1[j];
    }

/*     initialize the iteration counter.
     evaluate the function at the origin, and test
     for acceptance of the gauss-newton direction. */

    iter = 0;
    for (j = 1; j <= n; ++j) {
	wa2[j] = diag[j] * x[j];
    }
    dxnorm = enorm(n, &wa2[1]);
    fp = dxnorm - *delta;
    if (fp <= p1 * *delta) {
	goto L220;
    }

/*     if the jacobian is not rank deficient, the newton
     step provides a lower bound, parl, for the zero of
     the function. Otherwise set this bound to 0. */

    parl = 0.;
    if (nsing >= n) {
	for (j = 1; j <= n; ++j) {
	    l = ipvt[j];
	    wa1[j] = diag[l] * (wa2[l] / dxnorm);
	}
	for (j = 1; j <= n; ++j) {
	    sum = 0.;
	    for (i__ = 1; i__ <= j-1; ++i__) {
		sum += r__[i__ + j * r_dim1] * wa1[i__];
	    }
	    wa1[j] = (wa1[j] - sum) / r__[j + j * r_dim1];
	}
	temp = enorm(n, &wa1[1]);
	parl = fp / *delta / temp / temp;
    }

// L120:
/*     calculate an upper bound, paru, for the 0. of the function. */

    for (j = 1; j <= n; ++j) {
	sum = 0.;
	for (i__ = 1; i__ <= j; ++i__) {
	    sum += r__[i__ + j * r_dim1] * qtb[i__];
	}
	l = ipvt[j];
	wa1[j] = sum / diag[l];
    }
    gnorm = enorm(n, &wa1[1]);
    paru = gnorm / *delta;
    if (paru == 0.) {
	paru = dwarf / fmin2(*delta,p1);
    }

/*     if the input par lies outside of the interval (parl,paru),
     set par to the closer endpoint. */

    par = fmax2(par, parl);
    par = fmin2(par, paru);
    if (par == 0.) {
	par = gnorm / dxnorm;
    }

/*     beginning of an iteration. */

L150:
    ++iter;

/*        evaluate the function at the current value of par. */

    if (par == 0.)
	par = fmax2(dwarf, p001 * paru);

    temp = sqrt(par);
    for (j = 1; j <= n; ++j) {
	wa1[j] = temp * diag[j];
    }
    qrsolv(n, &r__[r_offset], ldr, &ipvt[1], &wa1[1], &qtb[1], &x[1],
	   &sdiag[1], &wa2[1]);
    for (j = 1; j <= n; ++j) {
	wa2[j] = diag[j] * x[j];
    }
    dxnorm = enorm(n, &wa2[1]);
    temp = fp;
    fp = dxnorm - *delta;

/*        if the function is small enough, accept the current value
        of par. also test for the exceptional cases where parl
        is 0. or the number of iterations has reached 10. */

    if (fabs(fp) <= p1 * *delta || (parl == 0. && fp <= temp && temp < 0.) ||
	iter == 10) { // << FIXME: give warning for  iter == 10 !!
	goto L220;
    }

/*        compute the newton correction. */

    for (j = 1; j <= n; ++j) {
	l = ipvt[j];
	wa1[j] = diag[l] * (wa2[l] / dxnorm);
    }
    for (j = 1; j <= n; ++j) {
	wa1[j] /= sdiag[j];
	temp = wa1[j];
	jp1 = j + 1;
	for (i__ = jp1; i__ <= n; ++i__) {
	    wa1[i__] -= r__[i__ + j * r_dim1] * temp;
	}
    }
    temp = enorm(n, &wa1[1]);
    parc = fp / *delta / temp / temp;

/*        depending on the sign of the function, update parl or paru. */

    if (fp > 0.) {
	parl = fmax2(parl,par);
    }
    if (fp < 0.) {
	paru = fmin2(paru,par);
    }

//        compute an improved estimate for par.

    par = fmax2(parl, par + parc);

/*        end of an iteration. */

    goto L150;

L220: //  termination.

    if (iter == 0) {
	par = 0.;
    }
    return par;
} /* lmpar */

/* Subroutine */
static
void qrsolv(int n, double *r__, int *ldr,
	    int *ipvt, double *diag, double *qtb, double *x,
	    double *sdiag, double *wa)
{
    /* Initialized data */

    static double p5 = .5;
    static double p25 = .25;

    /* System generated locals */
    int r_dim1, r_offset;

    /* Local variables */
    int i__, j, k, l, nsing;
    double tan__, cos__, sin__, temp, cotan, qtbpj;

/*     **********

     subroutine qrsolv

     given an m by n matrix a, an n by n diagonal matrix d,
     and an m-vector b, the problem is to determine an x which
     solves the system

           a*x = b ,     d*x = 0 ,

     in the least squares sense.

     this subroutine completes the solution of the problem
     if it is provided with the necessary information from the
     qr factorization, with column pivoting, of a. that is, if
     a*p = q*r, where p is a permutation matrix, q has orthogonal
     columns, and r is an upper triangular matrix with diagonal
     elements of nonincreasing magnitude, then qrsolv expects
     the full upper triangle of r, the permutation matrix p,
     and the first n components of (q transpose)*b. the system
     a*x = b, d*x = 0, is then equivalent to

                  t       t
           r*z = q *b ,  p *d*p*z = 0 ,

     where x = p*z. if this system does not have full rank,
     then a least squares solution is obtained. on output qrsolv
     also provides an upper triangular matrix s such that

            t   t               t
           p *(a *a + d*d)*p = s *s .

     s is computed within qrsolv and may be of separate interest.

     the subroutine statement is

       subroutine qrsolv(n,r,ldr,ipvt,diag,qtb,x,sdiag,wa)

     where

       n is a positive int input variable set to the order of r.

       r is an n by n array. on input the full upper triangle
         must contain the full upper triangle of the matrix r.
         on output the full upper triangle is unaltered, and the
         strict lower triangle contains the strict upper triangle
         (transposed) of the upper triangular matrix s.

       ldr is a positive int input variable not less than n
         which specifies the leading dimension of the array r.

       ipvt is an int input array of length n which defines the
         permutation matrix p such that a*p = q*r. column j of p
         is column ipvt(j) of the identity matrix.

       diag is an input array of length n which must contain the
         diagonal elements of the matrix d.

       qtb is an input array of length n which must contain the first
         n elements of the vector (q transpose)*b.

       x is an output array of length n which contains the least
         squares solution of the system a*x = b, d*x = 0.

       sdiag is an output array of length n which contains the
         diagonal elements of the upper triangular matrix s.

       wa is a work array of length n.

     subprograms called

       fortran-supplied ... fabs,dsqrt

     argonne national laboratory. minpack project. march 1980.
     burton s. garbow, kenneth e. hillstrom, jorge j. more

     **********
     Parameter adjustments */
    --wa;
    --sdiag;
    --x;
    --qtb;
    --diag;
    --ipvt;
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;

    /* Function Body

     copy r and (q transpose)*b to preserve input and initialize s.
     in particular, save the diagonal elements of r in x. */

    for (j = 1; j <= n; ++j) {
	for (i__ = j; i__ <= n; ++i__) {
	    r__[i__ + j * r_dim1] = r__[j + i__ * r_dim1];
	}
	x[j] = r__[j + j * r_dim1];
	wa[j] = qtb[j];
    }

/*     eliminate the diagonal matrix d using a givens rotation. */

    for (j = 1; j <= n; ++j) {

        /* prepare the row of d to be eliminated, locating the
	   diagonal element using p from the qr factorization. */

	l = ipvt[j];
	if (diag[l] == 0.) {
	    goto L90;
	}
	for (k = j; k <= n; ++k) {
	    sdiag[k] = 0.;
	}
	sdiag[j] = diag[l];

/*        the transformations to eliminate the row of d
        modify only a single element of (q transpose)*b
        beyond the first n, which is initially 0.. */

	qtbpj = 0.;
	for (k = j; k <= n; ++k)
	    if(sdiag[k] != 0.) {

		/* determine a givens rotation which eliminates the
		   appropriate element in the current row of d. */

		if (fabs(r__[k + k * r_dim1]) < fabs(sdiag[k])) {
		    cotan = r__[k + k * r_dim1] / sdiag[k];
		    sin__ = p5 / sqrt(p25 + p25 * (cotan * cotan));
		    cos__ = sin__ * cotan;
		} else {
		    tan__ = sdiag[k] / r__[k + k * r_dim1];
		    cos__ = p5 / sqrt(p25 + p25 * (tan__ * tan__));
		    sin__ = cos__ * tan__;
		}

/*           compute the modified diagonal element of r and
	     the modified element of ((q transpose)*b,0). */

		r__[k + k * r_dim1] = cos__ * r__[k + k * r_dim1] + sin__ * sdiag[k];
		temp  =  cos__ * wa[k] + sin__ * qtbpj;
		qtbpj = -sin__ * wa[k] + cos__ * qtbpj;
		wa[k] = temp;

/*           accumulate the tranformation in the row of s. */

		for (i__ = k+1; i__ <= n; ++i__) {
		    double r_n =  cos__ * r__[i__ + k * r_dim1] + sin__ * sdiag[i__];
		    sdiag[i__] = -sin__ * r__[i__ + k * r_dim1] + cos__ * sdiag[i__];
		    r__[i__ + k * r_dim1] = r_n;
		}
	    }

L90:

/*        store the diagonal element of s and restore
        the corresponding diagonal element of r. */

	sdiag[j] = r__[j + j * r_dim1];
	r__[j + j * r_dim1] = x[j];
    } // for( j = 1 .. n )

/*     solve the triangular system for z. if the system is
     singular, then obtain a least squares solution. */

    nsing = n;
    for (j = 1; j <= n; ++j) {
	if (sdiag[j] == 0. && nsing == n) {
	    nsing = j - 1;
	}
	if (nsing < n) {
	    wa[j] = 0.;
	}
    }

    for (k = 1; k <= nsing; ++k) {
	j = nsing - k + 1;
	double sum = 0.;
	for (i__ = j + 1; i__ <= nsing; ++i__)
	    sum += r__[i__ + j * r_dim1] * wa[i__];
	wa[j] = (wa[j] - sum) / sdiag[j];
    }

/*     permute the components of z back to components of x. */

    for (j = 1; j <= n; ++j) {
	l = ipvt[j];
	x[l] = wa[j];
    }
    return;
} /* qrsolv */

