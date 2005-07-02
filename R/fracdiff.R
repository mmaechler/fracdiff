####-*- mode: R; kept-old-versions: 12;  kept-new-versions: 20; -*-

### Original file:
### copyright 1991 Department of Statistics, Univeristy of Washington

### Patched by Friedrich.Leisch, for use with R, 22.1.1997
### fixed & changed by Martin Maechler, since Dec 2003

".First.lib" <- function(lib, pkg) library.dynam("fracdiff", pkg, lib)


fracdiff <- function(x, nar = 0, nma = 0,
                     ar = rep(NA, max(nar, 1)), ma = rep(NA, max(nma, 1)),
                     dtol = NULL, drange = c(0, 0.5), h, M = 100)
{
    ## #########################################################################
    ##
    ##   x      - time series for the ARIMA model
    ##   nar    - number of autoregressive parameters
    ##   nma    - number of moving average parameters
    ##   ar     - initial autoregressive parameters
    ##   ma     - initial moving average parameters
    ##   dtol   - desired accurcay for d
    ##            by default (and if negative), (4th root of machine precision)
    ##		  is used.  dtol will be changed internally if necessary
    ##   drange - interval over which the likelihood function is to be maximized
    ##            as a function of d
    ##   h      - finite difference interval
    ##   M      - number of terms in the likelihood approximation
    ##
    ##           (see Haslett and Raftery 1989)
    ##
    ## ########################################################################

    if(any(is.na(x)))
        stop("missing values not allowed in time series")
    if(is.matrix(x) && ncol(x) > 2)
        stop("multivariate time series not allowed")
    n <- length(x)
    if(round(nar) != nar || nar < 0 || round(nma) != nma || nma < 0)
        stop("'nar' and 'nma' must be non-negative integers")
    npq <- as.integer(nar + nma)
    npq1 <- as.integer(npq + 1:1)
    lenw <- max(npq + 2*(n + M),
                3*n + (n+6)*npq + npq %/% 2 + 1,
                31 * 12, ## << added because checked in ../src/fdcore.f
                (3 + 2*npq1) * npq1 + 1)## << this is *not* checked (there)
    lenw <- as.integer(lenw)
    ar[is.na(ar)] <- 0
    ma[is.na(ma)] <- 0
    if(is.null(dtol))
        dtol <- .Machine$double.eps^0.25 # ~ 1.22e-4
    ## if dtol < 0: the fortran code will choose defaults
    x <- as.double(x)

    ## this also initializes "common blocks" that are used
    ## in .Fortran(.) calls :
    result <- .Fortran("fracdf",
                       x,
                       n,
                       as.integer(M),
                       as.integer(nar),
                       as.integer(nma),
                       dtol = as.double(dtol),
                       drange = as.double(drange),
                       hood = double(1),
                       d = double(1),
                       ar = as.double(ar),
                       ma = as.double(ma),
                       w = double(lenw),
                       lenw = lenw,
                       iw = integer(npq),## <<< new int-work array
                       info = integer(1),
                       .Machine$double.xmin,
                       .Machine$double.xmax,
                       .Machine$double.neg.eps,
                       .Machine$double.eps,
                       PACKAGE = "fracdiff")
    if(result$info)
        switch(result$info,
               stop("insufficient workspace; need ", result$lenw,
                    " instead of just ", lenw),
               stop("error in gamma function"),
               stop("invalid MINPACK input"),
               warning("warning in gamma function"),
               warning("optimization failure"),
               warning("optimization limit reached"))

    hess <- .C("fdhpq",
               ## x,
               hess = double(npq1 * npq1),
               npq1,
               result$w,
               PACKAGE = "fracdiff")$hess

    temp <- .C("fdcov",
               x,
               result$d,
               h = as.double(if(missing(h)) -1 else h),
               hd = double(npq1),
               cov = hess, npq1,
               cor = hess, npq1,
               se = double(npq1),
               result$w,
               info = integer(1),
               PACKAGE = "fracdiff")
    if(temp$info)
        switch(temp$info,
               warning("warning in gamma function"),
               warning("singular Hessian"),
               warning("unable to compute correlation matrix"),
               stop("error in gamma function"))
    if(npq == 0) {
        result$ar <- NULL
        result$ma <- NULL
    }
    nam <- "d"
    if(nar) nam <- c(nam, paste("ar", 1:nar, sep = ""))
    if(nma) nam <- c(nam, paste("ma", 1:nma, sep = ""))
    hess <- matrix(hess, nrow = npq1, ncol = npq1, dimnames = list(nam, nam))
    hess[1, ] <- temp$hd
    hess[row(hess) > col(hess)] <- hess[row(hess) < col(hess)]
    se.ok <- temp$info != 0 || temp$info < 3
    list(log.likelihood = result$hood,
         d = result$d, ar = result$ar, ma = result$ma,
         covariance.dpq  = array(temp$cov, c(npq1, npq1), list(nam, nam)),
         stderror.dpq    = if(se.ok) temp$se, # else NULL
         correlation.dpq = if(se.ok) array(temp$cor, c(npq1, npq1)),
         h = temp$h, d.tol = result$dtol, M = M, hessian.dpq = hess,
         length.w = lenw)
}

fracdiff.var <- function(x, fracdiff.out, h)
{
    if(!is.numeric(h))
        stop("h must be numeric")
    if(!is.list(fracdiff.out) || !is.numeric(M <- fracdiff.out$M))
        stop("invalid ", sQuote("fracdiff.out"))
    p <- length(fracdiff.out$ar)
    q <- length(fracdiff.out$ma)
    n <- length(x)
    npq <- p + q
    npq1 <- npq + 1
    lwork <- max(npq + 2 * (n + M),
                 3 * n + (n + 6) * npq + npq %/% 2 + 1,
                 (3 + 2 * npq1) * npq1 + 1)
    ## Initialize
    .C("fdcom",
       n,
       as.integer(M),
       (p),
       (q),
       as.double(fracdiff.out$log.likelihood),
       .Machine$double.xmin,
       .Machine$double.xmax,
       .Machine$double.neg.eps,
       .Machine$double.eps,
       PACKAGE = "fracdiff")
    ## Re compute Covariance Matrix:
    temp <- .C("fdcov",
               as.double(x),
               as.double(fracdiff.out$d),
               h = as.double(h),
               hd = double(npq1),
               cov = as.double(fracdiff.out$hessian.dpq),
               as.integer(npq1),
               cor = as.double(fracdiff.out$hessian.dpq),
               as.integer(npq1),
               se = double(npq1),
               as.double(c(fracdiff.out$ma,
                           fracdiff.out$ar,
                           rep(0, lwork))),
               info = integer(1),
               PACKAGE = "fracdiff")
    if(temp$info)
        switch(temp$info,
               warning("warning in gamma function"),
               warning("singular Hessian"),
               warning("unable to compute correlation matrix"),
               stop("error in gamma function"))
    nam <- "d"
    if(p) nam <- c(nam, paste("ar", 1:p, sep = ""))
    if(q) nam <- c(nam, paste("ma", 1:q, sep = ""))
    fracdiff.out$h <- temp$h
    fracdiff.out$covariance.dpq <-
        matrix(temp$cov, nrow = npq1, ncol = npq1, dimnames = list(nam, nam))
    se.ok <- temp$info != 0 || temp$info < 3 ## << FIXME -- illogical!!
    fracdiff.out$stderror.dpq    <- if(se.ok) temp$se # else NULL
    fracdiff.out$correlation.dpq <- if(se.ok) array(temp$cor, c(npq1, npq1))
    fracdiff.out$hessian.dpq[1, ] <- temp$hd
    fracdiff.out$hessian.dpq[, 1] <- temp$hd
    fracdiff.out
}

## MM:  Added things for more  arima.sim() compatibility.
##      really, 'mu' is nonsense since can be done separately (or via 'innov').
fracdiff.sim <- function(n, ar = NULL, ma = NULL, d, rand.gen = rnorm,
                         innov = rand.gen(n+q, ...), n.start = NA,
                         allow.0.nstart = FALSE, # <- for back-compatibility
                         ..., mu = 0)
{
    p <- length(ar)
    q <- length(ma)
    if(p) {
        minroots <- min(Mod(polyroot(c(1, -ar))))
        if(minroots <= 1) {
            warning("'ar' part of fracdiff model is not stationary!!")
            minroots <- 1.01 # -> n.start= 603 by default
        }
    }
    if(is.na(n.start))
        n.start <- p + q + ifelse(p > 0, ceiling(6/log(minroots)), 0)
    if(n.start < p + q && !allow.0.nstart)
        stop("burn-in 'n.start' must be as long as 'ar + ma'")
    if(length(innov) < n+q) stop("'innov' must have length >= n + q")
    y <- c(rand.gen(n.start, ...), innov[1:(n+q)])
    stopifnot(is.double(y), length(y) == n + q + n.start)
    y <- .C("fdsim",
            as.integer(n + n.start),
            (p),
            (q),
            as.double(ar),
            as.double(ma),
            as.double(d),
            as.double(mu),
            y = y,
            s = double(length(y)),
            .Machine$double.xmin,
            .Machine$double.xmax,
            .Machine$double.neg.eps,
            .Machine$double.eps,
            PACKAGE = "fracdiff")[["s"]][n.start + 1:n]
    list(series = y, ar = ar, ma = ma, d = d, mu = mu, n.start = n.start)
}
