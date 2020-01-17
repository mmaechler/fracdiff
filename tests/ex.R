library(fracdiff)

.proctime00 <- proc.time()

set.seed(107)
options(digits = 5)

## 1)

x1 <- fracdiff.sim(5000, ar = .2, ma = -.4, d = .3, n.start=0, allow.0 = TRUE)
(fd1 <- fracdiff(x1$series, nar = 1, nma = 1, dtol = 1e-10))
vcov(fd1)
logLik(fd1)

fdCOVcomp <-
    c("h", "covariance.dpq", "stderror.dpq", "correlation.dpq", "hessian.dpq")
fd1. <- fracdiff.var(x1$series, fd1, h = fd1$h / 2)
dns <- dimnames(fd1.$covariance.dpq)

## dput(sapply(fd1.[fdCOVcomp], signif, digits = 4)) # edited:
fd1.L <- list(
    h = 3.7155e-05,
    covariance.dpq = matrix(c(0.0005966, -0.0008052, -0.0001897,
                              -0.0008052, 0.001612, 0.000664,
                              -0.0001897, 0.000664, 0.0005485),
                            3L, 3L, dimnames = dns),
    stderror.dpq = c(0.02443, 0.04015, 0.02342),
    correlation.dpq = matrix(c(1, -0.821, -0.3316,
                               -0.821, 1,  0.7061,
                               -0.3316, 0.7061, 1), 3),
    hessian.dpq = matrix(c(-8252, -5875, 4258,
                           -5875, -5420, 4529,
                           4258, 4529, -5834),
                         3L, 3L, dimnames = dns))
stopifnot(all.equal(fd1.[fdCOVcomp], fd1.L, tolerance = 2e-4))

fd1u <- fracdiff.var(x1$series, fd1, h = fd1$h * 8)
## dput(sapply(fd1u[fdCOVcomp], signif, digits = 4)) :
fd1uL <- list(
    h = 0.0005945,
    covariance.dpq = matrix(c(0.0005966, -0.0008052, -0.0001897,
                              -0.0008052, 0.001612, 0.000664,
                              -0.0001897, 0.000664, 0.0005485),
                            3L, 3L, dimnames = dns),
    stderror.dpq = c(0.02443, 0.04015, 0.02342),
    correlation.dpq = matrix(c(1, -0.821, -0.3316,
                               -0.821, 1, 0.7061,
                               -0.3316, 0.7061, 1), 3),
    hessian.dpq = matrix(c(-8252, -5875, 4258,
                           -5875, -5420, 4529,
                            4258,  4529, -5834),
                            3L, 3L, dimnames = dns))
stopifnot( all.equal(fd1u[fdCOVcomp], fd1uL, tolerance = 2e-4) )

## 2)

x2 <-  fracdiff.sim( 2048, ar = .8, ma = -.4, d = .3, n.start=0, allow.0 = TRUE)
## -> NA's and problems
fd2 <- fracdiff(x2$series, nar = length(x2$ar), nma = length(x2$ma))
summary(fd2)

fd2. <- fracdiff.var(x2$series, fd2, h = fd2$h / 2)
sapply(fd2.[fdCOVcomp], signif, digits = 4)
fd2u <- fracdiff.var(x2$series, fd2, h = fd2$h * 8)#-> warning, unable .. corr...
sapply(fd2u[fdCOVcomp], signif, digits = 4)

## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
