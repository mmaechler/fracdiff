#### Methods for "fracdiff"  objects
#### -------------------------------

coef.fracdiff <- function(object, ...) unlist(object[c("d", "ar", "ma")])

vcov.fracdiff      <- function(object, ...) object$covariance.dpq

residuals.fracdiff <- function(object, ...) .NotYetImplemented()
fitted.fracdiff    <- function(object, ...) .NotYetImplemented()

logLik.summary.fracdiff <-
logLik.fracdiff <- function(object, ...)
{
    r <- object$log.likelihood
    attr(r, "df") <- length(coef(object)) + 1 # "+ 1" : sigma^2
    attr(r, "nobs") <- attr(r, "nall") <- object$n
    class(r) <- "logLik"
    r
}

print.fracdiff <- function(x, digits = getOption("digits"), ...)
{
    cat("\nCall:\n ", deparse(x$call), "\n\nCoefficients:\n")
    print(coef(x), digits = digits, ...)
    ## print.default(x, digits = digits, ...)too cheap to be true
    if(x$msg != "ok") cat("*** Warning during fit:", x$msg,"\n")
    cat("a list with components:\n")
    print(names(x), ...)
    invisible(x)
}

summary.fracdiff <- function(object, symbolic.cor = FALSE, ...)
{
    ## add a 'coef' matrix (and not much more):
    cf <- coef(object)
    se <- object$stderror.dpq
    cf <- cbind("Estimate" = cf,
                "Std. Error"= se, "z value" = cf / se,
                "Pr(>|z|)" = 2 * pnorm(-abs(cf / se)))
    object$coef <- cf
    object$symbolic.cor <- symbolic.cor
    ## remove those components we have in 'coef' anyway
    object$d <- object$ar <- object$ma <- object$stderror.dpq <- NULL
    class(object) <- "summary.fracdiff"
    object
}

print.summary.fracdiff <-
    function(x, digits = max(3, getOption("digits") - 3),
             correlation = FALSE,
             symbolic.cor = x$symbolic.cor,
             signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n ", deparse(x$call), "\n\nCoefficients:\n")
    printCoefmat(x$coef, digits = digits, signif.stars = signif.stars, ...)
    if(x$msg != "ok")
	cat("*** Warning during fit:", x$msg,"\n")
    cat("[d.tol = ", formatC(x$d.tol),", M = ", x$M,", h = ",formatC(x$h),
	## really not much informative: "length.w = ", x$length.w,
	"]\n", sep='')
    print(logLik(x), digits = digits, ...)
    if (correlation && !is.null(correl <- x$correlation.dpq)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (is.logical(symbolic.cor) && symbolic.cor) {
                print(symnum(correl, abbr.col = NULL))
            }
            else {
                correl <- format(round(correl, 2), nsmall = 2, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop = FALSE], quote = FALSE)
            }
        }
    }
    invisible(x)
}


### This and coef.fracdiff()  were supplied

## From: Spencer Graves <spencer.graves@pdf.com>
## To: Melissa Ann Haltuch <mhaltuch@u.washington.edu>
## CC: r-help@stat.math.ethz.ch, Martin Maechler <maechler@stat.math.ethz.ch>
## Subject: Re: [R] fracdiff
## Date: Sun, 23 Jul 2006 03:40:08 +0800

confint.fracdiff <- function(object, parm, level = 0.95, ...)
{
    p <- length(cf <- coef(object))
    stopifnot(p >= 1, length(level) == 1, 0 < level, level < 1)
    se <- object$stderror.dpq
    pnames <- names(cf)
    names(se) <- pnames
    if (missing(parm))
        parm <- 1:p
    else if (is.character(parm))
        parm <- match(parm, pnames, nomatch = 0)
    cf <- cf[parm]
    se <- se[parm]

    a <- (1-level)/2
    a <- c(a, 1 - a)
    CI <- cf + outer(se, qnorm(a))
    dimnames(CI)[[2]] <- paste(format(100*a), "%")
    CI
}
