### Original file: 
### copyright 1991 Department of Statistics, Univeristy of Washington

### Patched by Friedrich.Leisch
### for use with R, 22.1.1997


".First.lib" <- function(lib, pkg) library.dynam("fracdiff", pkg, lib)

"fracdiff"<-
function(x, nar = 0, nma = 0, ar = rep(NA, max(nar, 1)), ma = rep(NA, max(nma, 
	1)), dtol, drange = c(0, 0.5), h, M = 100)
{
#############################################################################
# 
#   x      - time series for the ARIMA model
#   nar    - number of autoregressive parameters
#   nma    - number of moving average parameters
#   ar     - initial autoregressive parameters
#   ma     - initial moving average parameters
#   dtol   - desired accurcay for d
#            if dtol is negative, then the default 
#              (4th root of machine precision) is used
#               dtol will be changed internally if necessary
#   drange - interval over which the likelihood function is to be maximized
#            as a function of d
#   h      - finite difference interval 
#   M      - number of terms in the likelihood approximation 
#
#           (see Haslett and Raftery 1989)
#
#############################################################################
	n <- length(x)
	npq <- nar + nma
	npq1 <- npq + 1
	lwork <- max(npq + 2 * (n + M), 3 * n + (n + 6) * npq + npq %/% 2 + 1, (
		3 + 2 * npq1) * npq1 + 1)
	ar[is.na(ar)] <- 0
	ma[is.na(ma)] <- 0
	if(is.matrix(x) && ncol(x) > 2)
		stop("multivariate time series not allowed")
	if(any(is.na(x)))
		stop("missing values not allowed in time series")
	if(missing(dtol))
		dtol <- .Machine$double.eps^0.25
	result <- .Fortran("fracdf",
                           as.double(x),
                           as.integer(n),
                           as.integer(M),
                           as.integer(nar),
                           as.integer(nma),
                           dtol = as.double(dtol),
                           drange = as.double(drange),
                           hood = double(1),
                           d = double(1),
                           ar = as.double(ar),
                           ma = as.double(ma),
                           w = double(lwork),
                           as.integer(lwork),
                           info = integer(1),
                           as.double(.Machine$double.xmin),
                           as.double(.Machine$double.xmax),
                           as.double(.Machine$double.neg.eps),
                           as.double(.Machine$double.eps),
                           PACKAGE = "fracdiff")
	if(result$info)
		switch(result$info,
			stop("insufficient workspace"),
			stop("error in gamma function"),
			stop("invalid MINPACK input"),
			warning("warning in gamma function"),
			warning("optimization failure"),
			warning("optimization limit reached"))
	temp <- .Fortran("fdhpq",
                         as.double(x),
                         hess = double(npq1 * npq1),
                         as.integer(npq1),
                         double(lwork),
                         PACKAGE = "fracdiff")
	hess <- temp$hess
	temp <- .Fortran("fdcov",
                         as.double(x),
                         as.double(result$d),
                         h = as.double(if(missing(h)) -1 else h),
                         hd = double(npq1),
                         cov = as.double(hess),
                         as.integer(npq1),
                         cor = as.double(hess),
                         as.integer(npq1),
                         se = double(npq1),
                         double(lwork),
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
	if(nar)
		nam <- c(nam, paste("ar", 1:nar, sep = ""))
	if(nma)
		nam <- c(nam, paste("ma", 1:nma, sep = ""))
	hess <- matrix(hess, nrow = npq1, ncol = npq1, dimnames = list(nam, nam
		))
	hess[1,  ] <- temp$hd
	hess[row(hess) > col(hess)] <- hess[row(hess) < col(hess)]
	list(log.likelihood = result$hood, d = result$d, ar = result$ar, ma = 
		result$ma, covariance.dpq = array(temp$cov, c(npq1, npq1), list(
		nam, nam)), stderror.dpq = if(temp$info !=0 || temp$info < 3) temp$
			se else NULL, correlation.dpq = if(temp$info !=0 || temp$
		info < 3) array(temp$cor, c(npq1, npq1)) else NULL, h = temp$h, 
		d.tol = result$dtol, M = M, hessian.dpq = hess)
}
"fracdiff.var"<-
function(x, fracdiff.out, h)
{
	if(!is.numeric(h))
		stop("h must be numeric")
	p <- length(fracdiff.out$ar)
	q <- length(fracdiff.out$ma)
	n <- length(x)
	M <- fracdiff.out$M
	npq <- p + q
	npq1 <- npq + 1
	lwork <- max(npq + 2 * (n + M), 3 * n + (n + 6) * npq + npq %/% 2 + 1, (
		3 + 2 * npq1) * npq1 + 1)
	.Fortran("fdcom",
                 as.integer(n),
                 as.integer(M),
                 as.integer(p),
                 as.integer(q),
                 as.double(fracdiff.out$log.likelihood),
                 as.double(.Machine$double.xmin),
                 as.double(.Machine$double.xmax),
                 as.double(.Machine$double.neg.eps),
                 as.double(.Machine$double.eps),
                 PACKAGE = "fracdiff")
	temp <- .Fortran("fdcov",
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
	if(p)
		nam <- c(nam, paste("ar", 1:p, sep = ""))
	if(q)
		nam <- c(nam, paste("ma", 1:q, sep = ""))
	fracdiff.out$h <- temp$h
	fracdiff.out$covariance.dpq <- matrix(temp$cov, nrow = npq1, ncol = 
		npq1, dimnames = list(nam, nam))
	fracdiff.out$stderror.dpq <- if(!temp$info || temp$info < 3) temp$se
		 else NULL
	fracdiff.out$correlation.dpq <- if(!temp$info || temp$info < 3) array(
			temp$cor, c(npq1, npq1))
	fracdiff.out$hessian.dpq[1,  ] <- temp$hd
	fracdiff.out$hessian.dpq[, 1] <- temp$hd
	fracdiff.out
}
"fracdiff.sim"<-
function(n, ar, ma, d, mu = 0)
{
	p <- length(ar)
	q <- length(ma)
	temp <- .Fortran("fdsim",
                         as.integer(n),
                         as.integer(p),
                         as.integer(q),
                         as.double(ar),
                         as.double(ma),
                         as.double(d),
                         as.double(mu),
                         as.double(rnorm(n + q)),
                         x = double(n + q),
                         as.double(.Machine$double.xmin),
                         as.double(.Machine$double.xmax),
                         as.double(.Machine$double.neg.eps),
                         as.double(.Machine$double.eps),
                         PACKAGE = "fracdiff")$x[1:n]
	list(series = temp, ar = ar, ma = ma, d = d, mu = mu)
}


