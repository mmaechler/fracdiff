\name{fracdiff}
\alias{fracdiff}
\title{ML Estimates for Fractionally-Differenced ARIMA (p,d,q) models}
\description{
  Calculates the maximum likelihood estimators of the parameters
  of a fractionally-differenced ARIMA (p,d,q) model, together (if possible)
  with their estimated covariance and correlation matrices and standard
  errors, as well as the value of the maximized likelihood.  The
  likelihood is approximated using the fast and accurate method of
  Haslett and Raftery (1989).
}
\usage{
fracdiff(x, nar = 0, nma = 0,
         ar = rep(NA, max(nar, 1)), ma = rep(NA, max(nma, 1)),
         dtol = NULL, drange = c(0, 0.5), h, M = 100, trace = 0)
}
\arguments{
  \item{x}{time series (numeric vector) for the ARIMA model}
  \item{nar}{number of autoregressive parameters \eqn{p}.}
  \item{nma}{number of moving average parameters \eqn{q}.}
  \item{ar}{initial autoregressive parameters.}
  \item{ma}{initial moving average parameters.}
  \item{dtol}{interval of uncertainty for \eqn{d}.  If \code{dtol} is
    negative or NULL, the fourth root of machine precision will be used.
    \code{dtol} will be altered if necessary by the program.}
  \item{drange}{interval over which the likelihood function is to be
    maximized as a function of \eqn{d}.}
  \item{h}{size of finite difference interval for numerical derivatives.
    By default (or if negative),
%% Only found the following by reading ../src/fdhess.f :
    \code{h = min(0.1, eps.5 * (1+ abs(cllf)))}, where
    \code{clff := log. max.likelihood} (as returned) and
    \code{eps.5 := sqrt(.Machine$double.neg.eps)} (typically 1.05e-8).

    This is used to compute a finite difference approximation to the
    Hessian, and hence only influences the cov, cor, and std.error
    computations; use \code{\link{fracdiff.var}()} to change this
    \emph{after} the estimation process.
  }
  \item{M}{number of terms in the likelihood approximation (see Haslett
    and Raftery 1989).}
  \item{trace}{optional integer, specifying a trace level.  If positive,
    currently the \dQuote{outer loop} iterations produce one line of
    diagnostic output.}
}
\details{
  The \pkg{fracdiff} package has --- for historical reason, namely,
  S-plus \code{arima()} compatibility --- used an unusual
  parametrization for the MA part, see also the \sQuote{Details} section
  in \code{\link[stats]{arima}} (in standard \R's \pkg{stats} package).
  The ARMA (i.e., \eqn{d = 0}) model in \code{fracdiff()} and
  \code{\link{fracdiff.sim}()} is

  \deqn{X_t - a_1X_{t-1} - \cdots - a_pX_{t-p} = e_t - b_1e_{t-1} - \dots - b_qe_{t-q},}{%
    X[t] - a[1]X[t-1] - \dots - a[p]X[t-p] = e[t] - b[1]e[t-1] - \dots - b[q]e[t-q],}

  where \eqn{e_i}{e[i]} are mean zero i.i.d., for \code{fracdiff()}'s
  estimation,  \eqn{e_i \sim \mathcal{N}(0,\sigma^2)}{e[i] ~ N(0, s^2)}.
  This model indeed has the signs of the MA coefficients \eqn{b_j}{b[j]}
  \emph{inverted}, compared to other parametrizations, including
  Wikipedia's
  \url{https://en.wikipedia.org/wiki/Autoregressive_moving-average_model}
  and the one of \code{\link[stats]{arima}}.

  Note that \code{NA}'s in the initial values for \code{ar} or \code{ma}
  are replaced by \eqn{0}'s.
}
\value{
  an object of S3 \code{\link{class}} \code{"fracdiff"}, which is
  a list with components:
  \item{log.likelihood}{logarithm of the maximum likelihood}
  \item{d}{optimal fractional-differencing parameter}
  \item{ar}{vector of optimal autoregressive parameters}
  \item{ma}{vector of optimal moving average parameters}
  \item{covariance.dpq}{covariance matrix of the parameter estimates
    (order : d, ar, ma).}
  \item{stderror.dpq}{standard errors of the parameter estimates
  \code{c(d, ar, ma)}.}
  \item{correlation.dpq}{correlation matrix of the parameter estimates
    (order : d, ar, ma).}
  \item{h}{interval used for numerical derivatives, see \code{h} argument.}
  \item{dtol}{interval of uncertainty for d; possibly altered from input
    \code{dtol}.}
  \item{M}{as input.}
  \item{hessian.dpq}{the approximate Hessian matrix \eqn{H} of 2nd order
    partial derivatives of the likelihood with respect to the
    parameters; this is (internally) used to compute
    \code{covariance.dpq}, the approximate asymptotic covariance matrix as
    \eqn{C = (-H)^{-1}}.}
}
\note{
  Ordinarily, \code{nar} and \code{nma} should not be too large (say < 10)
  to avoid degeneracy in the model.  The function
  \code{\link{fracdiff.sim}} is available for generating test problems.
}
\section{Method}{
  The optimization is carried out in two levels:\cr
  an outer univariate unimodal
  optimization in d over the interval \code{drange} (typically [0,.5]),
  using Brent's \code{fmin} algorithm), and\cr
  an inner nonlinear least-squares optimization in the AR and MA parameters to
  minimize white noise variance (uses the MINPACK subroutine \code{lm}DER).
  written by Chris Fraley (March 1991).
}
\section{Warning}{
  The variance-covariance matrix and consequently the standard errors
  may be quite inaccurate, see the example in \code{\link{fracdiff.var}}.
}
\references{
  J. Haslett and A. E. Raftery (1989)
  Space-time Modelling with Long-memory Dependence: Assessing Ireland's
  Wind Power Resource (with Discussion);
  \emph{Applied Statistics} \bold{38}, 1--50.

  R. Brent (1973)
  \emph{Algorithms for Minimization without Derivatives}, Prentice-Hall

  J. J. More, B. S. Garbow, and K. E. Hillstrom (1980)
  \emph{Users Guide for MINPACK-1}, Technical Report ANL-80-74,
  Applied Mathematics Division, Argonne National Laboratory.
}
\seealso{
  \code{\link{coef.fracdiff}} and other methods for \code{"fracdiff"}
  objects;
  \code{\link{fracdiff.var}()} for re-estimation of variances or
  standard errors;
  \code{\link{fracdiff.sim}}
}
\examples{
ts.test <- fracdiff.sim( 5000, ar = .2, ma = -.4, d = .3)
fd. <- fracdiff( ts.test$series,
                 nar = length(ts.test$ar), nma = length(ts.test$ma))
fd.
## Confidence intervals
confint(fd.)

## with iteration output
fd2 <- fracdiff(ts.test$series, nar = 1, nma = 1, trace = 1)
all.equal(fd., fd2)
}
\keyword{ts}
