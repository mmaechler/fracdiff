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
         dtol,
         drange = c(0, 0.5), h, M = 100)
}
\arguments{
  \item{x}{time series (numeric vector) for the ARIMA model}
  \item{nar}{number of autoregressive parameters \eqn{p}.}
  \item{nma}{number of moving average parameters \eqn{q}.}
  \item{ar}{initial autoregressive parameters}
  \item{ma}{initial moving average parameters}
  \item{dtol}{interval of uncertainty for \eqn{d}. If \code{dtol < 0},
    the fourth root of machine precision will be used.  \code{dtol} will
    be altered if necessary by the program.}
  \item{drange}{interval over which the likelihood function is to be
    maximized as a function of \eqn{d}.}
  \item{h}{finite difference interval.}
  \item{M}{number of terms in the likelihood approximation (see Haslett
    and Raftery 1989).}
}
\value{
  a list containing the following elements :
  \item{log.likelihood}{logarithm of the maximum likelihood}
  \item{d}{optimal fractional-differencing parameter}
  \item{ar}{vector of optimal autoregressive parameters}
  \item{ma}{vector of optimal moving average parameters}
  \item{covariance.dpq}{covariance matrix of the parameter estimates
    (order : d, ar, ma).}
  \item{stderror.dpq}{standard errors of the parameter estimates  c(d, ar, ma).}
  \item{correlation.dpq}{correlation matrix of the parameter estimates
    (order : d, ar, ma).}
  \item{dtol}{interval of uncertainty for d.}
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
  written by Chris Fraley (March 1991)
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
  \code{\link{fracdiff.sim}}
}
\examples{
ts.test <- fracdiff.sim( 5000, ar = .2, ma = -.4, d = .3)
fracdiff( ts.test$series, nar = length(ts.test$ar), nma = length(ts.test$ma))
}
\keyword{ts}
