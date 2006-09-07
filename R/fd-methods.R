
### These two were supplied
## From: Spencer Graves <spencer.graves@pdf.com>
## To: Melissa Ann Haltuch <mhaltuch@u.washington.edu>
## CC: r-help@stat.math.ethz.ch, Martin Maechler <maechler@stat.math.ethz.ch>
## Subject: Re: [R] fracdiff
## Date: Sun, 23 Jul 2006 03:40:08 +0800

coef.fracdiff <- function(object) {
    unlist(object[c("d", "ar", "ma")])
}

confint.fracdiff <- function(object, parm, level=0.95, ...)
{
    p <- length(b <- coef.fracdiff(object))
    stopifnot(p >= 1)
    se <- object$stderror.dpq
    names(se) <- names(b)

    if(!missing(parm)) {
        b <- b[parm]
        se <- se[parm]
    }

    conf.c <- (1-level)/2
    conf2 <- c(conf.c, 1-conf.c)
    CI <- b + outer(se, qnorm(conf2))
    dimnames(CI)[[2]] <- paste(format(100*conf2), "%")
    CI
}
