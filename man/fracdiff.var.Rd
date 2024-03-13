\name{fracdiff.var}
\alias{fracdiff.var}
\title{Recompute Covariance Estimate for fracdiff}
\usage{
fracdiff.var(x, fracdiff.out, h)
}
\arguments{
  \item{x}{a univariate time series or a vector.  Missing values (NAs)
    are not allowed.}
  \item{fracdiff.out}{output from \code{fracdiff} for time series \code{x}.}
  \item{h}{finite-difference interval length (\eqn{ > 0}) for approximating partial
    derivatives with respect to the \code{d} parameter.  Typically smaller
    than the one in \code{fracdiff.out}}
}
\description{
  Allows the finite-difference interval to be altered for recomputation of the
  covariance estimate for \code{fracdiff}.
}
\value{
  an object of S3 \code{\link{class}} \code{"fracdiff"}, i.e., basically
  a list with the same elements as the result from
  \code{\link{fracdiff}}, but with possibly different values for the
  hessian, covariance, and correlation matrices and for standard error,
  as well as for \code{h}.
}
\seealso{
  \code{fracdiff}, also for references.
}
\examples{
## Generate a fractionally-differenced ARIMA(1,d,1) model :
set.seed(5) # reproducibility; x86_64 Lnx: get warning
tst <- fracdiff.sim(500, ar = .2, ma = .4, d = .3)$series
## estimate the parameters in an ARIMA(1,d,1) model for the simulated series
fd.out <- fracdiff(tst, nar= 1, nma = 1) # warning ... maybe change 'h'
summary(fd.out)## *** Warning ... {has been stored} --> h = 7.512e-6

## Modify the covariance estimate by changing the finite-difference interval
(fd.o2 <- fracdiff.var(tst, fd.out, h = 1e-3))
## looks identical as  print(fd.out),
## however these (e.g.) differ :
vcov(fd.out)
vcov(fd.o2)

## A case, were the default variance is *clearly* way too small:
set.seed(1); fdc <- fracdiff(X <- fracdiff.sim(n=100, d=0.25)$series)
fdc
# Confidence intervals just based on asymp.normal approx. and std.errors:
confint(fdc) # ridiculously too narrow
}
\keyword{ts}
