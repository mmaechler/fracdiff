\name{fracdiff-methods}
\alias{coef.fracdiff}
\alias{logLik.fracdiff}
\alias{print.fracdiff}
\alias{fitted.fracdiff}
\alias{residuals.fracdiff}
\alias{vcov.fracdiff}
\alias{summary.fracdiff}
\alias{print.summary.fracdiff}
%
\title{Many Methods for "fracdiff" Objects}
\description{
  Many \dQuote{accessor} methods for \code{\link{fracdiff}} objects,
  notably \code{\link{summary}}, \code{\link{coef}}, \code{\link{vcov}}, and
  \code{\link{logLik}}; further \code{\link{print}()} methods were needed.
}
\usage{
\method{coef}{fracdiff}(object, \dots)
\method{logLik}{fracdiff}(object, \dots)
\method{print}{fracdiff}(x, digits = getOption("digits"), \dots)
\method{summary}{fracdiff}(object, symbolic.cor = FALSE, \dots)
\method{print}{summary.fracdiff}(x, digits = max(3, getOption("digits") - 3),
        correlation = FALSE, symbolic.cor = x$symbolic.cor,
        signif.stars = getOption("show.signif.stars"), \dots)
\method{fitted}{fracdiff}(object, \dots)
\method{residuals}{fracdiff}(object, \dots)
\method{vcov}{fracdiff}(object, \dots)
}
\arguments{
  \item{x, object}{object of class \code{fracdiff}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{\dots}{further arguments passed from and to methods.}
  \item{correlation}{logical; if \code{TRUE}, the correlation matrix of
    the estimated parameters is returned and printed.}
  \item{symbolic.cor}{logical. If \code{TRUE}, print the correlations in
    a symbolic form (see \code{\link{symnum}}) rather than as numbers.}
  \item{signif.stars}{logical. If \code{TRUE}, \dQuote{significance stars}
    are printed for each coefficient.}
}
\author{Martin Maechler; Rob Hyndman contributed the
  \code{\link{residuals}()} and \code{\link{fitted}()} methods.}
\seealso{\code{\link{fracdiff}} to get \code{"fracdiff"} objects,
  \code{\link{confint.fracdiff}} for the \code{\link{confint}} method;
  further, \code{\link{fracdiff.var}}.
}
\examples{
set.seed(7)
ts4 <- fracdiff.sim(10000, ar = c(0.6, -.05, -0.2), ma = -0.4, d = 0.2)
modFD <- fracdiff( ts4$series, nar = length(ts4$ar), nma = length(ts4$ma))
## -> warning (singular Hessian) %% FIXME ???
coef(modFD) # the estimated parameters
vcov(modFD)
smFD <- summary(modFD)
smFD
coef(smFD) # gives the whole table
AIC(modFD) # AIC works because of the logLik() method
stopifnot(exprs = {

})
}
\keyword{print}
\keyword{models}
