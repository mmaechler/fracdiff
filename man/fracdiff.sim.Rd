\name{fracdiff.sim}
\alias{fracdiff.sim}
\title{fracdiff.sim: test problem generator for fracdiff}

\usage{
fracdiff.sim( n, ar, ma, d, mu = 0.)

}
\arguments{
  \item{n}{length of the time series }
  \item{ar}{vector of autoregressive parameters }
  \item{ma}{vector of moving average parameters}
  \item{d}{fractional differencing parameter}
  \item{mu}{time series mean}
}

\description{
  Generates simulated long-memory time series data from the
  ARIMAS(p,d,q) model.
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

\references{
  J. Haslett and A. E. Raftery, "Space-time Modelling with Long-memory
  Dependence: Assessing Ireland's Wind Power Resource (with Discussion)",
  Applied Statistics, 38, 1-50.
}

\seealso{
  \code{fracdiff}
}

\examples{
fracdiff.sim( 10000, ar = .2, ma = .4, d = .3)
}

\keyword{ts}
