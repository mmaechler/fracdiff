\name{fracdiff.var}
\alias{fracdiff.var}
\title{Recompute Covariance Estimate for fracdiff}
\usage{
fracdiff.var(x, fracdiff.out, h)
}
\arguments{
  \item{x}{a univariate time series or a vector. Missing values (NAs)
    are not allowed.} 
\item{fracdiff.out}{output from \code{fracdiff} for time series
  \code{x}.}
\item{h}{finite-difference interval for approximating partial
  derivatives with respect to the \code{d} parameter.}
}
\description{
  Allows the finite-difference interval to be altered for recomputation of the
  covariance estimate for \code{fracdiff}.
}
\value{
  a list with the same elements as the output to fracdiff, but with 
  possibly different values for the hessian, covariance, and correlation
  matrices and for standard error, as well as for \code{h}.
}
\seealso{
  \code{fracdiff}
}

\examples{
# generate a fractionally-differenced ARIMA(1,d,1)
# model given initial values
ts.test <- fracdiff.sim( 10000, ar = .2, ma = .4, d = .3)
# estimate the parameters in an ARIMA(1,d,1)
# model for the simulated series
fd.out <- fracdiff(ts.test$ser, nar=length(ts.test$ar),
                   nma=length(ts.test$ma))
# modify the covariance estimate by changing the
# finite-difference interval
fracdiff.var(ts.test$series, fd.out, h = .0001) 
}
\keyword{ts}
