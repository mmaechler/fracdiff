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

\value{\item{sdiff}{The differenciated time series data}}

\references{
  \item{Geweke, J. and Porter-Hudak, S. (1983). The estimation and application of
        long memory time series models. Journal of Time Series Analysis, 4(4), 221-238}
  \item{Reisen, V. A. (1994). Estimation of the fractional difference parameter in the
        ARFIMA(p,d,q) model using the smoothed periodogram. Journal Time Series Analysis,
	    15(1), 335-350}
  \item{Reisen, V. A., B. Abraham, and E. M. M. Toscano, (2001). Parametric and semiparametric
        estimations of stationary univariate ARFIMA model. Brazilian Journal of Probability
        and Statistics, 14, 185-206.}
  \item{Reisen, V. A. and Lopes, S. (1999). Some simulations and applications of forecasting long-memory time series models,
       Journal of Statistical Planning and Inference, 80, 269-287.}
  \item{Reisen, V. A. Cribari-Neto, F. and Jensen, M.J. (2003). Long Memory Inflationary Dynamics. The case of Brazil,
       Studies in Nonlinear Dynamics and Econometrics, 7(3), 1-16.}
           }

\author{Valderio A. Reisen (valderio@cce.ufes.br) and Artur J. Lemonte}

\seealso{fracdiff.sim}

\examples{
memory.long <- fracdiff.sim(80, d = 0.3)
mGPH <- GPH(memory.long$series)$d
r <- diffseries(memory.long$series, d = mGPH$d)
}
\keyword{ts}
