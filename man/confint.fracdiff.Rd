\name{confint.fracdiff}
\alias{confint.fracdiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
confint.fracdiff(object, parm, level = 0.95, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ ~~Describe \code{object} here~~ }
  \item{parm}{ ~~Describe \code{parm} here~~ }
  \item{level}{ ~~Describe \code{level} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (object, parm, level = 0.95, ...) 
{
    p <- length(cf <- coef(object))
    stopifnot(p >= 1, length(level) == 1, 0 < level, level < 
        1)
    se <- object$stderror.dpq
    pnames <- names(cf)
    names(se) <- pnames
    if (missing(parm)) 
        parm <- 1:p
    else if (is.character(parm)) 
        parm <- match(parm, pnames, nomatch = 0)
    cf <- cf[parm]
    se <- se[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    CI <- cf + outer(se, qnorm(a))
    dimnames(CI)[[2]] <- paste(format(100 * a), "\%")
    CI
  }
}
\keyword{ ~kwd1 }% at least one, from doc/KEYWORDS
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
