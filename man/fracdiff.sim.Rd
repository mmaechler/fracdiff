\name{fracdiff.sim}
\alias{fracdiff.sim}
\title{Simulate fractional ARIMA Time Series}
\description{
  Generates simulated long-memory time series data from the
  fractional ARIMA(p,d,q) model.  This is a test problem generator for
  \code{\link{fracdiff}}.
}
\usage{
fracdiff.sim(n, ar, ma, d, mu = 0.)
}
\arguments{
  \item{n}{length of the time series.}
  \item{ar}{vector of autoregressive parameters.}
  \item{ma}{vector of moving average parameters.}
  \item{d}{fractional differencing parameter.}
  \item{mu}{time series mean}
}
\value{
  a list containing the following elements :
  \item{x}{time series}
  \item{ar}{same as input}
  \item{ma}{same as input}
  \item{d}{same as input}
  \item{mu}{same as input}
  \item{seed}{same as input}
}
\seealso{
  \code{\link{fracdiff}}, also for references.
}
\examples{
## Pretty (too) short to "see" the long memory
fracdiff.sim(100, ar = .2, ma = .4, d = .3)
}
\keyword{ts}
