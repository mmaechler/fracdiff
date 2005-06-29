/* fdmin.f -- translated by f2c (version 20031025).
 *
 * and produced by
 * $Id: f2c-clean,v 1.10 2002/03/28 16:37:27 maechler Exp $
 *
 * and manually pretty edited by Martin Maechler, 2004-10-01
 */

#include <Rmath.h>

#include "fracdiff.h"

#ifndef max
# define	max(a, b) 		((a) < (b) ? (b) : (a))
#endif
#ifndef min
# define	min(a, b)		((a) > (b) ? (b) : (a))
#endif
#ifndef abs
# define	abs(x)			((x) >= 0 ? (x) : -(x))
#endif


/* Common Block Declarations --- included as "extern" */
#define FD_EXTERNAL extern

#include "mach_comm.h"
#include "maux_comm.h"

struct {
    double told, tolf, tolx, tolg, fnorm, delta, gnorm;
} tolsfd_;

#define tolsfd_1 tolsfd_

/* Table of constant values */

static /*logical*/int c_true = (1);

/* --------- EXPORTS (need all ??) ------------------- */

double enorm_(int *, double *);
/* Subroutine */
int lmpar_(int *, double *, int *, int *, double *, double *, double *,
	   double *, double *, double *, double *, double *);
static void qrfac_(int *, int *, double *, int *,
		   /*logical*/int *, int *, int *,
		   double *, double *, double *);
static void qrsolv_(int *, double *, int *, int *,
		    double *, double *, double *, double *, double *);

/* ------------------------------- */


/* Subroutine */
int lmder1_(S_fp fcn, int *m, int *n, double *x,
	    double *fvec, double *fjac, int *ldfjac, double *ftol,
	    double *xtol, double *gtol, int *maxfev, double *diag,
	    int *mode, double *factor, int *info,
	    int *nfev, int *njev, int *ipvt, double *qtf,
	    double *wa1, double *wa2, double *wa3, double *wa4, double *y)
{
    /* Initialized data */

    static double one = 1.;
    static double p1 = .1;
    static double p5 = .5;
    static double p25 = .25;
    static double p75 = .75;
    static double p0001 = 1e-4;
    static double zero = 0.;

    /* System generated locals */
    int fjac_dim1, fjac_offset;
    double d__1, d__2;

    /* Local variables */
    int i__, j, l, iter, iflag, nprint;
    double par, sum, temp, temp1, temp2;
    double ratio, enorm_n, xnorm, fnorm1, actred, dirder, prered;

/*     **********

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

         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
         int m,n,ldfjac,iflag
         double precision x(n),fvec(m),fjac(ldfjac,n)
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

                t     t           t
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

     **********
     double precision dpmpar,enorm

     epsmch is the machine precision.

     epsmch = dpmpar(1)

     Parameter adjustments */
    --wa4;
    --fvec;
    --wa3;
    --wa2;
    --wa1;
    --qtf;
    --ipvt;
    --diag;
    --x;
    fjac_dim1 = *ldfjac;
    fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;
    --y;

    /* Function Body */
    temp = 0.;
    nprint = 0;
    *info = 0;
    iflag = 0;
    *nfev = 0;
    *njev = 0;
/*     check the input parameters for errors. */

    if (*n <= 0 || *m < *n || *ldfjac < *m || *ftol < zero || *xtol < zero ||

	    *gtol < zero || *maxfev <= 0 || *factor <= zero) {
	goto L300;
    }
    if (*mode != 2) {
	goto L20;
    }
    for (j = 1; j <= *n; ++j) {
	if (diag[j] <= zero) {
	    goto L300;
	}
    }
L20:

/*     evaluate the function at the starting point
     and calculate its norm. */

    iflag = 1;
    (*fcn)(&x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag, &y[1]);
    *nfev = 1;
    if (iflag < 0) {
	goto L300;
    }

    tolsfd_1.fnorm = fmin2(enorm_(m, &fvec[1]), mauxfd_1.bignum);

/*     initialize levenberg-marquardt parameter and iteration counter. */

    par = zero;
    iter = 1;

/*     beginning of the outer loop. */

L30:

/*        calculate the jacobian matrix. */

    iflag = 2;
    (*fcn)(&x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag, &y[1]);
    ++(*njev);
    if (iflag < 0) {
	goto L300;
    }

/*        if requested, call fcn to enable printing of iterates. */

    if (nprint <= 0) {
	goto L40;
    }
    iflag = 0;
    if ((iter - 1) % nprint == 0) {
	(*fcn)(&x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag, &y[1]);
    }
    if (iflag < 0) {
	goto L300;
    }
L40:

/*        compute the qr factorization of the jacobian. */

    qrfac_(m, n, &fjac[fjac_offset], ldfjac, &c_true, &ipvt[1], n, &wa1[1], &
	   wa2[1], &wa3[1]);

/*        on the first iteration and if mode is 1, scale according
        to the norms of the columns of the initial jacobian. */

    if (iter != 1) {
	goto L80;
    }
    if (*mode == 2) {
	goto L60;
    }
    for (j = 1; j <= *n; ++j) {
	diag[j] = wa2[j];
	if (wa2[j] == zero) {
	    diag[j] = one;
	}
    }
L60:

/*        on the first iteration, calculate the norm of the scaled x
        and initialize the step bound delta. */

    for (j = 1; j <= *n; ++j) {
	wa3[j] = diag[j] * x[j];
    }
    xnorm = enorm_(n, &wa3[1]);
    tolsfd_1.delta = *factor * xnorm;
    if (tolsfd_1.delta == zero) {
	tolsfd_1.delta = *factor;
    }
L80:

/*        form (q transpose)*fvec and store the first n components in
        qtf. */

    for (i__ = 1; i__ <= *m; ++i__) {
	wa4[i__] = fvec[i__];
    }
    for (j = 1; j <= *n; ++j) {
	if (fjac[j + j * fjac_dim1] == zero) {
	    goto L120;
	}
	sum = zero;
	for (i__ = j; i__ <= *m; ++i__) {
	    sum += fjac[i__ + j * fjac_dim1] * wa4[i__];
	}
	temp = -sum / fjac[j + j * fjac_dim1];
	for (i__ = j; i__ <= *m; ++i__) {
	    wa4[i__] += fjac[i__ + j * fjac_dim1] * temp;
	}
L120:
	fjac[j + j * fjac_dim1] = wa1[j];
	qtf[j] = wa4[j];
    }

/*        compute the norm of the scaled gradient. */

    tolsfd_1.gnorm = zero;
    if (tolsfd_1.fnorm == zero) {
	goto L170;
    }
    for (j = 1; j <= *n; ++j) {
	l = ipvt[j];
	if (wa2[l] != zero) {
	    sum = zero;
	    for (i__ = 1; i__ <= j; ++i__) {
		sum += fjac[i__ + j * fjac_dim1] * (qtf[i__] / tolsfd_1.fnorm);
	    }
	    tolsfd_1.gnorm = fmax2(tolsfd_1.gnorm, fabs(sum / wa2[l]));
	}
    }
L170:

/*        test for convergence of the gradient norm. */

    if (tolsfd_1.gnorm <= *gtol) {
	*info = 4;
    }
    if (*info != 0) {
	goto L300;
    }

/*        rescale if necessary. */

    if (*mode == 2) {
	goto L190;
    }
    for (j = 1; j <= *n; ++j)
	diag[j] = fmax2(diag[j], wa2[j]);

L190:

/*        beginning of the inner loop. */

L200:
/*           determine the levenberg-marquardt parameter. */

    lmpar_(n, &fjac[fjac_offset], ldfjac, &ipvt[1], &diag[1], &qtf[1], &
	    tolsfd_1.delta, &par, &wa1[1], &wa2[1], &wa3[1], &wa4[1]);

/*           store the direction p and x + p. calculate the norm of p. */

    for (j = 1; j <= *n; ++j) {
	wa1[j] = -wa1[j];
	wa2[j] = x[j] + wa1[j];
	wa3[j] = diag[j] * wa1[j];
    }
    enorm_n = enorm_(n, &wa3[1]);

/*           on the first iteration, adjust the initial step bound. */

    if (iter == 1) {
	tolsfd_1.delta = min(tolsfd_1.delta,enorm_n);
    }

/*           evaluate the function at x + p and calculate its norm. */

    iflag = 1;
    (*fcn)(&wa2[1], &wa4[1], &fjac[fjac_offset], ldfjac, &iflag, &y[1]);
    ++(*nfev);
    if (iflag < 0) {
	goto L300;
    }
    fnorm1 = fmin2(enorm_(m, &wa4[1]), mauxfd_1.bignum);

/*           compute the scaled actual reduction. */

    actred = -one;
    if (p1 * fnorm1 < tolsfd_1.fnorm) {
/* Computing 2nd power */
	d__1 = fnorm1 / tolsfd_1.fnorm;
	actred = one - d__1 * d__1;
    }
/*          actred = (fnorm*fnorm - fnorm1*fnorm1)

           compute the scaled predicted reduction and
           the scaled directional derivative. */

    for (j = 1; j <= *n; ++j) {
	wa3[j] = zero;
	l = ipvt[j];
	temp = wa1[l];
	for (i__ = 1; i__ <= j; ++i__) {
	    wa3[i__] += fjac[i__ + j * fjac_dim1] * temp;
	}
    }
    temp1 = enorm_(n, &wa3[1]) / tolsfd_1.fnorm;
    temp2 = sqrt(par) * enorm_n / tolsfd_1.fnorm;
/* Computing 2nd power */
    d__1 = temp1;
/* Computing 2nd power */
    d__2 = temp2;
    prered = d__1 * d__1 + d__2 * d__2 / p5;
/*           temp1  = enorm(n,wa3)
           temp2  = (dsqrt(par)*enorm_n)
           prered = (temp1**2 + 2.d0*temp2**2)
 Computing 2nd power */
    d__1 = temp1;
/* Computing 2nd power */
    d__2 = temp2;
    dirder = -(d__1 * d__1 + d__2 * d__2);

/*           compute the ratio of the actual to the predicted
           reduction. */

    ratio = zero;
    if (prered != zero) {
	ratio = actred / prered;
    }

/*           update the step bound. */

    if (ratio > p25) {
	goto L240;
    }
    if (actred >= zero) {
	temp = p5;
    }
    if (actred < zero) {
	temp = p5 * dirder / (dirder + p5 * actred);
    }
    if (p1 * fnorm1 >= tolsfd_1.fnorm || temp < p1) {
	temp = p1;
    }
    tolsfd_1.delta = temp * fmin2(tolsfd_1.delta, enorm_n / p1);
    par /= temp;
    goto L260;
L240:
    if (par != zero && ratio < p75) {
	goto L250;
    }
    tolsfd_1.delta = enorm_n / p5;
    par = p5 * par;
L250:
L260:

/*           test for successful iteration. */

    if (ratio < p0001) {
	goto L290;
    }

/*           successful iteration. update x, fvec, and their norms. */


    for (j = 1; j <= *n; ++j) {
	x[j] = wa2[j];
	wa2[j] = diag[j] * x[j];
    }
    for (i__ = 1; i__ <= *m; ++i__) {
	fvec[i__] = wa4[i__];
    }
    xnorm = enorm_(n, &wa2[1]);
    tolsfd_1.fnorm = fnorm1;
    ++iter;

L290:
/*           tests for convergence. */

    if (abs(actred) <= *ftol && prered <= *ftol && p5 * ratio <= one) {
	*info = 1;
    }
    if (tolsfd_1.fnorm <= *ftol) {
	*info = 1;
    }
    if (tolsfd_1.delta <= *xtol) {
	*info = 2;
    }
    if (abs(actred) <= *ftol && prered <= *ftol && p5 * ratio <= one &&
	*info == 2) {
	*info = 3;
    }
    if (*info != 0) {
	goto L300;
    }

/*           tests for termination and stringent tolerances. */

    if (*nfev >= *maxfev) {
	*info = 5;
    }
    if (fabs(actred) <= machfd_.epsmax && prered <= machfd_.epsmax &&
	p5 * ratio <= one) {
	*info = 6;
    }
    if (tolsfd_1.delta <= machfd_.epsmax) {
	*info = 7;
    }
    if (tolsfd_1.gnorm <= machfd_.epsmax) {
	*info = 8;
    }
    if (*info != 0) {
	goto L300;
    }

/*           end of the inner loop. repeat if iteration unsuccessful. */

    if (ratio < p0001) {
	goto L200;
    }

/*        end of the outer loop. */

    goto L30;
L300:

/*     termination, either normal or user imposed. */

    if (iflag < 0) {
	*info = iflag;
    }
    iflag = 0;
    if (nprint > 0) {
	(*fcn)(&x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag, &y[1]);
    }
    return 0;
} /* lmder1_ */


double enorm_(int *n, double *x)
{
    /* Initialized data */

    static double one = 1.;
    static double zero = 0.;
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
    s1 = zero;
    s2 = zero;
    s3 = zero;
    x1max = zero;
    x3max = zero;
    floatn = (double) (*n);
    agiant = rgiant / floatn;

    for (i__ = 1; i__ <= *n; ++i__) {
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
		s1 = one + s1 * (d__1 * d__1);
		x1max = xabs;
	    }

	} else { /*				sum for small components. */

	    if (xabs <= x3max) {

		if (xabs != zero) {
		    /* Computing 2nd power */
		    d__1 = xabs / x3max;
		    s3 += d__1 * d__1;
		}

	    } else {
		/* Computing 2nd power */
		d__1 = x3max / xabs;
		s3 = one + s3 * (d__1 * d__1);
		x3max = xabs;
	    }
	}

    } /* for(i ) */

/*     calculation of norm. */

    if (s1 == zero) {

	if (s2 == zero) {
	    ret_val = x3max * sqrt(s3);
	}
	else {
	    if (s2 >= x3max)
		ret_val = sqrt(s2 * (one + x3max / s2 * (x3max * s3)));
	    else /* (s2 < x3max) */
		ret_val = sqrt(x3max * (s2 / x3max + x3max * s3));
	}
    }
    else {
	ret_val = x1max * sqrt(s1 + s2 / x1max / x1max);
    }

    return ret_val;
} /* enorm_ */


static
void qrfac_(int *m, int *n, double *a, int *lda,
	    /*logical*/int *pivot, int *ipvt, int *lipvt, double *rdiag,
	    double *acnorm, double *wa)
{
    /* Initialized data */

    static double one = 1.;
    static double p05 = .05;
    static double zero = 0.;

    /* System generated locals */
    int a_dim1, a_offset, i__2, i__3;
    double d__1, d__2, d__3;

    /* Local variables */
    int i__, j, k, kmax, jp1, minmn;
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

    /* Function Body
     epsmch is the machine precision.

     epsmch = dpmpar(1)

     compute the initial column norms and initialize several arrays. */

    for (j = 1; j <= *n; ++j) {
	acnorm[j] = enorm_(m, &a[j * a_dim1 + 1]);
	rdiag[j] = acnorm[j];
	wa[j] = rdiag[j];
	if (*pivot) {
	    ipvt[j] = j;
	}
    }

/*     reduce a to r with householder transformations. */

    minmn = min(*m,*n);
    for (j = 1; j <= minmn; ++j) {
	if (! (*pivot)) {
	    goto L40;
	}

/*        bring the column of largest norm into the pivot position. */

	kmax = j;
	for (k = j; k <= *n; ++k) {
	    if (rdiag[k] > rdiag[kmax]) {
		kmax = k;
	    }
	}
	if (kmax == j) {
	    goto L40;
	}
	for (i__ = 1; i__ <= *m; ++i__) {
	    temp = a[i__ + j * a_dim1];
	    a[i__ + j * a_dim1] = a[i__ + kmax * a_dim1];
	    a[i__ + kmax * a_dim1] = temp;
	}
	rdiag[kmax] = rdiag[j];
	wa[kmax] = wa[j];
	k = ipvt[j];
	ipvt[j] = ipvt[kmax];
	ipvt[kmax] = k;
L40:

/*        compute the householder transformation to reduce the
        j-th column of a to a multiple of the j-th unit vector. */

	i__2 = *m - j + 1;
	ajnorm = enorm_(&i__2, &a[j + j * a_dim1]);

	if (ajnorm == zero) {
	    goto L100;
	}
	if (a[j + j * a_dim1] < zero) {
	    ajnorm = -ajnorm;
	}
	for (i__ = j; i__ <= *m; ++i__) {
	    a[i__ + j * a_dim1] /= ajnorm;
	}
	a[j + j * a_dim1] += one;

/*        apply the transformation to the remaining columns
        and update the norms. */

	jp1 = j + 1;
	if (*n < jp1) {
	    goto L100;
	}
	for (k = jp1; k <= *n; ++k) {
	    sum = zero;
	    for (i__ = j; i__ <= *m; ++i__) {
		sum += a[i__ + j * a_dim1] * a[i__ + k * a_dim1];
	    }
	    temp = sum / a[j + j * a_dim1];
	    for (i__ = j; i__ <= *m; ++i__) {
		a[i__ + k * a_dim1] -= temp * a[i__ + j * a_dim1];
	    }
	    if (! (*pivot) || rdiag[k] == zero) {
		goto L80;
	    }
	    temp = a[j + k * a_dim1] / rdiag[k];
/* Computing MAX
 Computing 2nd power */
	    d__3 = temp;
	    d__1 = zero, d__2 = one - d__3 * d__3;
	    rdiag[k] *= sqrt((max(d__1,d__2)));
/* Computing 2nd power */
	    d__1 = rdiag[k] / wa[k];
	    if (p05 * (d__1 * d__1) > machfd_.epsmax) {
		goto L80;
	    }
	    i__3 = *m - j;
	    rdiag[k] = enorm_(&i__3, &a[jp1 + k * a_dim1]);
	    wa[k] = rdiag[k];
L80:
	    ;
	}
L100:
	rdiag[j] = -ajnorm;
    }
    return;
} /* qrfac_ */

/* Subroutine */
int lmpar_(int *n, double *r__, int *ldr,
	   int *ipvt, double *diag, double *qtb, double *delta,
	   double *par, double *x, double *sdiag,
	   double *wa1, double *wa2)
{
    /* Initialized data */

    static double p1 = .1;
    static double p001 = .001;
    static double zero = 0.;

    /* System generated locals */
    int r_dim1, r_offset;
    double d__1, d__2;

    /* Local variables */
    static int i__, j, k, l;
    static double fp;
    static int jm1, jp1;
    static double sum, parc, parl;
    static int iter;
    static double temp, paru, dwarf;
    static int nsing;
    static double gnorm, dxnorm;

/*     **********

     subroutine lmpar

     given an m by n matrix a, an n by n nonsingular diagonal
     matrix d, an m-vector b, and a positive number delta,
     the problem is to determine a value for the parameter
     par such that if x solves the system

           a*x = b ,     sqrt(par)*d*x = 0 ,

     in the least squares sense, and dxnorm is the euclidean
     norm of d*x, then either par is zero and

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

     **********
     double precision dpmpar,enorm
     Parameter adjustments */
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

    /* Function Body
     dwarf is the smallest positive magnitude.

     dwarf = dpmpar(2) */
    dwarf = machfd_.fltmin;

/*     compute and store in x the gauss-newton direction. if the
     jacobian is rank-deficient, obtain a least squares solution. */

    nsing = *n;
    for (j = 1; j <= *n; ++j) {
	wa1[j] = qtb[j];
	if (r__[j + j * r_dim1] == zero && nsing == *n) {
	    nsing = j - 1;
	}
	if (nsing < *n) {
	    wa1[j] = zero;
	}
    }
    if (nsing < 1) {
	goto L50;
    }
    for (k = 1; k <= nsing; ++k) {
	j = nsing - k + 1;
	wa1[j] /= r__[j + j * r_dim1];
	temp = wa1[j];
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L30;
	}
	for (i__ = 1; i__ <= jm1; ++i__) {
	    wa1[i__] -= r__[i__ + j * r_dim1] * temp;
	}
L30:
	;
    }
L50:
    for (j = 1; j <= *n; ++j) {
	l = ipvt[j];
	x[l] = wa1[j];
    }

/*     initialize the iteration counter.
     evaluate the function at the origin, and test
     for acceptance of the gauss-newton direction. */

    iter = 0;
    for (j = 1; j <= *n; ++j) {
	wa2[j] = diag[j] * x[j];
    }
    dxnorm = enorm_(n, &wa2[1]);
    fp = dxnorm - *delta;
    if (fp <= p1 * *delta) {
	goto L220;
    }

/*     if the jacobian is not rank deficient, the newton
     step provides a lower bound, parl, for the zero of
     the function. otherwise set this bound to zero. */

    parl = zero;
    if (nsing < *n) {
	goto L120;
    }
    for (j = 1; j <= *n; ++j) {
	l = ipvt[j];
	wa1[j] = diag[l] * (wa2[l] / dxnorm);
    }
    for (j = 1; j <= *n; ++j) {
	sum = zero;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L100;
	}
	for (i__ = 1; i__ <= jm1; ++i__) {
	    sum += r__[i__ + j * r_dim1] * wa1[i__];
	}
L100:
	wa1[j] = (wa1[j] - sum) / r__[j + j * r_dim1];
    }
    temp = enorm_(n, &wa1[1]);
    parl = fp / *delta / temp / temp;
L120:

/*     calculate an upper bound, paru, for the zero of the function. */

    for (j = 1; j <= *n; ++j) {
	sum = zero;
	for (i__ = 1; i__ <= j; ++i__) {
	    sum += r__[i__ + j * r_dim1] * qtb[i__];
	}
	l = ipvt[j];
	wa1[j] = sum / diag[l];
    }
    gnorm = enorm_(n, &wa1[1]);
    paru = gnorm / *delta;
    if (paru == zero) {
	paru = dwarf / min(*delta,p1);
    }

/*     if the input par lies outside of the interval (parl,paru),
     set par to the closer endpoint. */

    *par = max(*par,parl);
    *par = min(*par,paru);
    if (*par == zero) {
	*par = gnorm / dxnorm;
    }

/*     beginning of an iteration. */

L150:
    ++iter;

/*        evaluate the function at the current value of par. */

    if (*par == zero) {
/* Computing MAX */
	d__1 = dwarf, d__2 = p001 * paru;
	*par = max(d__1,d__2);
    }
    temp = sqrt(*par);
    for (j = 1; j <= *n; ++j) {
	wa1[j] = temp * diag[j];
    }
    qrsolv_(n, &r__[r_offset], ldr, &ipvt[1], &wa1[1], &qtb[1], &x[1], &sdiag[
	    1], &wa2[1]);
    for (j = 1; j <= *n; ++j) {
	wa2[j] = diag[j] * x[j];
    }
    dxnorm = enorm_(n, &wa2[1]);
    temp = fp;
    fp = dxnorm - *delta;

/*        if the function is small enough, accept the current value
        of par. also test for the exceptional cases where parl
        is zero or the number of iterations has reached 10. */

    if (abs(fp) <= p1 * *delta || (parl == zero && fp <= temp && temp < zero) ||
	iter == 10) {
	goto L220;
    }

/*        compute the newton correction. */

    for (j = 1; j <= *n; ++j) {
	l = ipvt[j];
	wa1[j] = diag[l] * (wa2[l] / dxnorm);
    }
    for (j = 1; j <= *n; ++j) {
	wa1[j] /= sdiag[j];
	temp = wa1[j];
	jp1 = j + 1;
	if (*n < jp1) {
	    goto L200;
	}
	for (i__ = jp1; i__ <= *n; ++i__) {
	    wa1[i__] -= r__[i__ + j * r_dim1] * temp;
	}
L200:
	;
    }
    temp = enorm_(n, &wa1[1]);
    parc = fp / *delta / temp / temp;

/*        depending on the sign of the function, update parl or paru. */

    if (fp > zero) {
	parl = max(parl,*par);
    }
    if (fp < zero) {
	paru = min(paru,*par);
    }

/*        compute an improved estimate for par.

 Computing MAX */
    d__1 = parl, d__2 = *par + parc;
    *par = max(d__1,d__2);

/*        end of an iteration. */

    goto L150;
L220:

/*     termination. */

    if (iter == 0) {
	*par = zero;
    }
    return 0;
} /* lmpar_ */

/* Subroutine */
static
void qrsolv_(int *n, double *r__, int *ldr,
	     int *ipvt, double *diag, double *qtb, double *x,
	     double *sdiag, double *wa)
{
    /* Initialized data */

    static double p5 = .5;
    static double p25 = .25;
    static double zero = 0.;

    /* System generated locals */
    int r_dim1, r_offset;
    double d__1, d__2;

    /* Local variables */
    int i__, j, k, l, jp1, kp1, nsing;
    double tan__, cos__, sin__, sum, temp, cotan, qtbpj;

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

    for (j = 1; j <= *n; ++j) {
	for (i__ = j; i__ <= *n; ++i__) {
	    r__[i__ + j * r_dim1] = r__[j + i__ * r_dim1];
	}
	x[j] = r__[j + j * r_dim1];
	wa[j] = qtb[j];
    }

/*     eliminate the diagonal matrix d using a givens rotation. */

    for (j = 1; j <= *n; ++j) {

        /* prepare the row of d to be eliminated, locating the
	   diagonal element using p from the qr factorization. */

	l = ipvt[j];
	if (diag[l] == zero) {
	    goto L90;
	}
	for (k = j; k <= *n; ++k) {
	    sdiag[k] = zero;
	}
	sdiag[j] = diag[l];

/*        the transformations to eliminate the row of d
        modify only a single element of (q transpose)*b
        beyond the first n, which is initially zero. */

	qtbpj = zero;
	for (k = j; k <= *n; ++k) {

	/* determine a givens rotation which eliminates the
	   appropriate element in the current row of d. */

	    if (sdiag[k] == zero) {
		goto L70;
	    }
	    if (fabs(r__[k + k * r_dim1]) >= (d__2 = sdiag[k],

		    abs(d__2))) {
		goto L40;
	    }
	    cotan = r__[k + k * r_dim1] / sdiag[k];
/* Computing 2nd power */
	    d__1 = cotan;
	    sin__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	    cos__ = sin__ * cotan;
	    goto L50;
L40:
	    tan__ = sdiag[k] / r__[k + k * r_dim1];
/* Computing 2nd power */
	    d__1 = tan__;
	    cos__ = p5 / sqrt(p25 + p25 * (d__1 * d__1));
	    sin__ = cos__ * tan__;
L50:

/*           compute the modified diagonal element of r and
           the modified element of ((q transpose)*b,0). */

	    r__[k + k * r_dim1] = cos__ * r__[k + k * r_dim1] + sin__ * sdiag[
		    k];
	    temp = cos__ * wa[k] + sin__ * qtbpj;
	    qtbpj = -sin__ * wa[k] + cos__ * qtbpj;
	    wa[k] = temp;

/*           accumulate the tranformation in the row of s. */

	    kp1 = k + 1;
	    if (*n < kp1) {
		goto L70;
	    }
	    for (i__ = kp1; i__ <= *n; ++i__) {
		temp = cos__ * r__[i__ + k * r_dim1] + sin__ * sdiag[i__];
		sdiag[i__] = -sin__ * r__[i__ + k * r_dim1] + cos__ * sdiag[
			i__];
		r__[i__ + k * r_dim1] = temp;
	    }
L70:
	    ;
	}
L90:

/*        store the diagonal element of s and restore
        the corresponding diagonal element of r. */

	sdiag[j] = r__[j + j * r_dim1];
	r__[j + j * r_dim1] = x[j];
    }

/*     solve the triangular system for z. if the system is
     singular, then obtain a least squares solution. */

    nsing = *n;
    for (j = 1; j <= *n; ++j) {
	if (sdiag[j] == zero && nsing == *n) {
	    nsing = j - 1;
	}
	if (nsing < *n) {
	    wa[j] = zero;
	}
    }
    if (nsing < 1) {
	goto L150;
    }
    for (k = 1; k <= nsing; ++k) {
	j = nsing - k + 1;
	sum = zero;
	jp1 = j + 1;
	if (nsing < jp1) {
	    goto L130;
	}
	for (i__ = jp1; i__ <= nsing; ++i__) {
	    sum += r__[i__ + j * r_dim1] * wa[i__];
	}
L130:
	wa[j] = (wa[j] - sum) / sdiag[j];
    }
L150:

/*     permute the components of z back to components of x. */

    for (j = 1; j <= *n; ++j) {
	l = ipvt[j];
	x[l] = wa[j];
    }
    return;
} /* qrsolv_ */

