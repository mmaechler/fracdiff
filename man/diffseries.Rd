\name{diffseries}
\alias{diffseries}
\title{Fractionally Differenciate Data}

\description{
  Differenciates the time series data using
  the approximated binomial expression of the long-memory filter and an estimate of
  the memory parameter in the ARFIMA(p,d,q) model.
}
\usage{
diffseries(x, d)
}
\arguments{
  \item{x}{numeric vector or univariate time series.}
  \item{d}{number specifiying the fractional difference order.}
}

\value{the fractionally differenced series \code{x}.}
\details{
  Since 2018, we are using (an important correction of) the fast
  algorithm based on the discrete Fourier transform (\code{\link{fft}})
  by Jensen and Nielsen which is significantly faster for large
  \code{n = length(x)}.
}
\references{
  See those in \code{\link{fdSperio}}; additionally

  Reisen, V. A. and Lopes, S. (1999)
  Some simulations and applications
  of forecasting long-memory time series  models;
  \emph{Journal of Statistical Planning and Inference} \bold{80}, 269--287.

  Reisen, V. A. Cribari-Neto, F. and Jensen, M.J. (2003)
  Long Memory Inflationary Dynamics. The case of Brazil.
  \emph{Studies in Nonlinear Dynamics and Econometrics} \bold{7}(3), 1--16.

  Jensen, Andreas Noack and Nielsen, Morten \enc{Ørregaard}{Oerregaard} (2014)
  A Fast Fractional Difference Algorithm.
  \emph{Journal of Time Series Analysis} \bold{35}(5), 428--436;
  \doi{10.1111/jtsa.12074}.
}

\author{Valderio A. Reisen \email{valderio@cce.ufes.br} and Artur
  J. Lemonte (first slow version), now hidden as \code{diffseries.0()}.

  Current version: Jensen and Nielsen (2014); tweaks by Martin Maechler, 2018.
}

\seealso{\code{\link{fracdiff.sim}}}

\examples{
memory.long <- fracdiff.sim(80, d = 0.3)
str(mGPH <- fdGPH(memory.long$series))
r <- diffseries(memory.long$series, d = mGPH$d)
#acf(r) # shouldn't show structure - ideally
}
\keyword{ts}
