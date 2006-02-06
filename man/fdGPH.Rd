\name{fdGPH}
\alias{fdGPH}
\title{Geweke and Porter-Hudak Estimator for ARFIMA(p,d,q)}

\description{
  Estimate the fractional (or \dQuote{memory}) parameter \eqn{d} in the
  ARFIMA(p,d,q) model by the method of Geweke and Porter-Hudak (GPH).
  The GPH estimator is based on the regression equation using the
  peridogram function as an estimate of the spectral density.
}
\usage{
fdGPH(x, bandw.exp = 0.5)
}
\arguments{
  \item{x}{univariate time series}
  \item{bandw.exp}{the bandwidth used in the regression equation}
}
\details{
  The function also provides the asymptotic standard deviation and the standard
  error deviation of the fractional estimator.

  The bandwidth is
  \code{bw = trunc(n ^ bandw.exp)}, where 0 < bandw.exp < 1 and n is the sample size.
  Default \code{bandw.exp = 0.5}.
}

\value{
  \item{d}{GPH estimate}
  \item{sd.as}{asymptotic standard deviation}
  \item{sd.reg}{standard error deviation}
}

\references{see those in \code{\link{fdSperio}}.
}

\author{Valderio A. Reisen and Artur J. Lemonte}

\seealso{\code{\link{fdSperio}}, \code{\link{fracdiff}}}

\examples{
memory.long <- fracdiff.sim(1500, d = 0.3)
fdGPH(memory.long$series)
}
\keyword{ts}
