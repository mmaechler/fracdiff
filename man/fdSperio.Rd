\name{Sperio}
\alias{Sperio}

\title{Estimate the fractional parameter in the ARFIMA(p,d,q) model using
  the smoothed periodogram function
}

\description{
  This function makes use Reisen (1994) estimator to estimate the memory
  parameter d in the ARFIMA model. It is based on the regression
  equation using the smoothed peridogram function as an estimate of the
  spectral density.
}

\usage{Sperio(x, band = NULL, beta = NULL)}

\arguments{
  \item{x}{The univariate time series data}
  \item{band}{The bandwidth used in the regression equation}
  \item{beta}{The bandwidth used in the lag Parzen window}
}

\details{
  \item{The function also provides the asymptotic standard deviation and the standard
        error deviation of the fractional estimator}
  \item{band = trunc(n^(alpha)), where 0 < alpha < 1 and n is the sample size.
        Default alpha = 0.5}
  \item{beta = trunc(n^(beta)), where 0 < beta < 1 and n is the sample size.
        Default beta = 0.9}
        }

\value{
  \item{d}{Sperio estimate}
  \item{sd.as}{asymptotic standard deviation}
  \item{sd.reg}{standard error deviation}
}
\references{
  Geweke, J. and Porter-Hudak, S. (1983)
  The estimation and application of long memory time series models.
  \emph{Journal of Time Series Analysis} \bold{4}(4), 221--238.

  Reisen, V. A. (1994).
  Estimation of the fractional difference parameter in the ARFIMA(p,d,q)
  model using the smoothed periodogram.
  \emph{Journal Time Series Analysis}, \bold{15}(1), 335--350.

  Reisen, V. A., B. Abraham, and E. M. M. Toscano (2001).
  Parametric and semiparametric estimations of stationary univariate
  ARFIMA model.
  \emph{Brazilian Journal of Probability and Statistics} \bold{14}, 185--206.
}
\author{Valderio A. Reisen \email{valderio@cce.ufes.br} and Artur J. Lemonte}

\seealso{GPH, fracdiff}

\examples{
memory.long <- fracdiff.sim(1500, d = 0.3)
Sperio(memory.long$series)
}
\keyword{ts}
