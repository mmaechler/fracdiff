\name{diffseries}
\alias{diffseries}

\title{Differenciate the data using an estimate of the fractional parameter}

\description{
  This function differenciates the time series data using
  the approximated binomial expression of the long-memory filter and an estimate of
  the memory parameter in the ARFIMA(p,d,q) model.
}

\usage{diffseries(x, d)}

\arguments{
  \item{x}{numeric vector or univariate time series}.
  \item{d}{number specifiying the fractional difference order.}
}

\value{the fractionally differenced series \code{x}.}

\references{
  See those in \code{\link{Sperio}}; additionally

  Reisen, V. A. and Lopes, S. (1999)
  Some simulations and applications
  of forecasting long-memory time series  models;
  \emph{Journal of Statistical Planning and Inference} \bold{80}, 269--287.

  Reisen, V. A. Cribari-Neto, F. and Jensen, M.J. (2003)
  Long Memory Inflationary Dynamics. The case of Brazil.
  \emph{Studies in Nonlinear Dynamics and Econometrics} \bold{7}(3), 1--16.
}

\author{Valderio A. Reisen (valderio@cce.ufes.br) and Artur J. Lemonte}

\seealso{\code{\link{fracdiff.sim}}}

\examples{
memory.long <- fracdiff.sim(80, d = 0.3)
mGPH <- GPH(memory.long$series)
r <- diffseries(memory.long$series, d = mGPH$d)
acf(r) # shouldn't show structure - ideally
}
\keyword{ts}
