library(fracdiff)

doExtras <- fracdiff:::doExtras()
.proctime00 <- proc.time()

set.seed(107)
options(digits = 5)

## 1)

x1 <- fracdiff.sim(5000, ar = .2, ma = -.4, d = .3, n.start=0, allow.0 = TRUE)
summary(fd1 <- fracdiff(x1$series, nar = 1, nma = 1, dtol = 1e-10))
vcov(fd1)
logLik(fd1)
stopifnot(all.equal(structure(-7051.5027, df = 4L, nall = 5000L, nobs = 5000L, class = "logLik"),
                    logLik(fd1)))
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
if(doExtras)
    print(all.equal(fd1u[fdCOVcomp], fd1uL, tolerance = 0))
stopifnot(all.equal(fd1u[fdCOVcomp], fd1uL, tolerance = 2e-4) )

## 2)

x2 <-  fracdiff.sim( 2048, ar = .8, ma = -.4, d = .3, n.start=0, allow.0 = TRUE)
## -> NA's and problems
fd2 <- fracdiff(x2$series, nar = length(x2$ar), nma = length(x2$ma))
sfd2 <- summary(fd2)
ss2 <- lapply(sfd2[setdiff(names(sfd2), c("residuals", "fitted"))],
              function(.) if(is.numeric(.)) signif(., digits = 7) else .)
ss2$coefficients <- ss2$coefficients[, -4] # drop p values
ss2S <- list(
    log.likelihood = -2924.262, n = 2048,
    msg = c(fracdf = "ok", fdcov = "ok"),
    covariance.dpq = matrix(c(0.0004182859, -0.0007078449, -6.753008e-05,
                             -0.0007078449,  0.001712827,  0.0002692938,
                             -6.753008e-05,  0.0002692938, 0.0002572701), 3L,
                            dimnames = dns),
    fnormMin = 45.62935, sigma = 1.008768,
    correlation.dpq = matrix(c(1, -0.8362667, -0.2058572,
                               -0.8362667,  1, 0.405672,
                               -0.2058572, 0.405672,   1), 3, dimnames = dns),
    h = 3.082264e-05, d.tol = 0.0001220703, M = 100,
    hessian.dpq = matrix(c(-8557.83, -3810.55,  1742.32,
                           -3810.55, -2395.564, 1507.303,
                           1742.32, 1507.303, -5007.377), 3, dimnames = dns),
    length.w = 10254,
    call = quote(
        fracdiff(x = x2$series, nar = length(x2$ar), nma = length(x2$ma))),
    coefficients = matrix(c(0.3374173,  0.7709664, -0.3810478,
                            0.02045204, 0.04138631, 0.01603964,
                            16.49798,  18.62854,  -23.75663),
                          3, 3, dimnames = dimnames(ss2$coefficients)),
    df = 4, aic = 5856.524, symbolic.cor = FALSE)
##
if(doExtras)
    print(all.equal(ss2S, ss2, tol = 0)) # 0.0001273 (32b Win); TRUE (64b F30, gcc)
stopifnot(all.equal(ss2S, ss2, tol = 4e-4))

fd2. <- fracdiff.var(x2$series, fd2, h = fd2$h / 2)
sfd2. <-         sapply(fd2.[fdCOVcomp], signif, digits = 4)
sfd2S <- ## dput(sapply(fd2.[fdCOVcomp], signif, digits = 5))
    list(h = 1.5411e-05,
         covariance.dpq = matrix(c( 5.4726e-05,-9.261e-05, -8.8353e-06,
                                   -9.261e-05,  0.0006717,  0.00016997,
                                   -8.8353e-06, 0.00016997, 0.00024779), 3, dimnames=dns),
         stderror.dpq = c(0.0073977, 0.025917, 0.015741),
         correlation.dpq = matrix(c(1, -0.48303, -0.075871,
                                    -0.48303, 1, 0.41661,
                                    -0.075871, 0.41661, 1), 3),
         hessian.dpq = matrix(c(-24440,  -3810.6, 1742.3,
                                -3810.6, -2395.6, 1507.3,
                                 1742.3,  1507.3,-5007.4), 3, dimnames=dns))
##
if(doExtras)
    print(all.equal(sfd2S, sfd2., tol =   0 , countEQ=TRUE)) # 8.7655e-5
stopifnot(all.equal(sfd2S, sfd2., tol = 2e-4, countEQ=TRUE))

fd2u <- fracdiff.var(x2$series, fd2, h = fd2$h * 8)#-> warning, unable .. corr...
##= no se.ok -->
fdCOV.0 <- setdiff(fdCOVcomp, c("stderror.dpq", "correlation.dpq"))
sd2u <- sapply(fd2u[fdCOV.0], signif, digits = 4)
sd2uS <- list(  ## dput(sapply(sd2u[fdCOVcomp], signif, digits = 5))
    h = 0.0002466,
    covariance.dpq = matrix(c(-0.0003545, 6e-04, 5.724e-05,
                              6e-04, -0.0005003, 5.816e-05,
                              5.724e-05, 5.816e-05, 0.0002371), 3, dimnames=dns),
    ## stderror.dpq = c(0, 0, 0.0154),
    ## correlation.dpq = matrix(0, 3,3),
    hessian.dpq = matrix(c(-3347, -3811, 1742,
                           -3811, -2396, 1507,
                            1742,  1507,-5007), 3, dimnames=dns))
##
if(doExtras)
    print(all.equal(sd2uS, sd2u, tol =   0 , countEQ=TRUE))# 0.000103 (32b Win); T.(64b F30)
stopifnot(all.equal(sd2uS, sd2u, tol = 4e-4, countEQ=TRUE))

