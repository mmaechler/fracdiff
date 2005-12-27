\name{GPH}
\alias{GPH}
\title{Geweke and Porter-Hudak Estimator for ARFIMA(p,d,q)}

\description{
  Estimate the fractional (or \dQuote{memory}) parameter \eqn{d} in the
  ARFIMA(p,d,q) model by the method of Geweke and Porter-Hudak (GPH).
  The GPH estimator is based on the regression equation using the
  peridogram function as an estimate of the spectral density.
}
\usage{GPH(x, band = NULL)}
\arguments{
  \item{x}{univariate time series}
  \item{band}{the bandwidth used in the regression equation}
}
\details{
  The function also provides the asymptotic standard deviation and the standard
  error deviation of the fractional estimator.

  \code{band = trunc(n^(alpha))}, where \code{0 < alpha < 1} and \eqn{n}
  is the sample size.
  Default \code{alpha = 0.5}.
%% MM: FIXME: I don't understand the above!
}

\value{
  \item{d}{GPH estimate}
  \item{sd.as}{asymptotic standard deviation}
  \item{sd.reg}{standard error deviation}
}

\references{see those in \code{\link{Sperio}}.
}

\author{Valderio A. Reisen and Artur J. Lemonte}

\seealso{\code{\link{Sperio}}, \code{\link{fracdiff}}}

\examples{
memory.long <- fracdiff.sim(1500, d = 0.3)
GPH(memory.long$series)
}
\keyword{ts}
