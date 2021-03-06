\name{fdSperio}
\alias{fdSperio}

\title{Sperio Estimate for 'd' in ARFIMA(p,d,q)}
\description{
  This function makes use Reisen (1994) estimator to estimate the memory
  parameter d in the ARFIMA(p,d,q) model.  It is based on the regression
  equation using the smoothed periodogram function as an estimate of the
  spectral density.
}
\usage{
fdSperio(x, bandw.exp = 0.5, beta = 0.9)
}
\arguments{
  \item{x}{univariate time series data.}
  \item{bandw.exp}{numeric: exponent of the bandwidth used in the regression equation.}
  \item{beta}{numeric: exponent of the bandwidth used in the lag Parzen window.}
}
\details{
  The function also provides the asymptotic standard deviation and the
  standard error deviation of the fractional estimator.

  The bandwidths are \code{bw = trunc(n ^ bandw.exp)}, where 0 < bandw.exp < 1
  and n is the sample size.  Default \code{bandw.exp= 0.5};
  \cr
  and \code{bw2 = trunc(n ^ beta)}, where 0 < beta < 1 and n is the
  sample size.  Default \code{beta = 0.9}.
}

\value{
  a list with components
  \item{d}{Sperio estimate}
  \item{sd.as}{asymptotic standard deviation}
  \item{sd.reg}{standard error deviation}
}
\references{
  Geweke, J. and Porter-Hudak, S. (1983)
  The estimation and application of long memory time series models.
  \emph{Journal of Time Series Analysis} \bold{4}(4), 221--238.

  Reisen, V. A. (1994)
  Estimation of the fractional difference parameter in the ARFIMA(p,d,q)
  model using the smoothed periodogram.
  \emph{Journal Time Series Analysis}, \bold{15}(1), 335--350.

  Reisen, V. A., B. Abraham, and E. M. M. Toscano (2001)
  Parametric and semiparametric estimations of stationary univariate
  ARFIMA model.
  \emph{Brazilian Journal of Probability and Statistics} \bold{14}, 185--206.
}
\author{Valderio A. Reisen \email{valderio@cce.ufes.br} and Artur J. Lemonte}

\seealso{\code{\link{fdGPH}}, \code{\link{fracdiff}}
}

\examples{
memory.long <- fracdiff.sim(1500, d = 0.3)
spm <- fdSperio(memory.long$series)
str(spm, digits=6)
}
\keyword{ts}
