### Very similar to ./ex.R  but using *MORE* precision
###  ===> non-portable but better for consistency checking (Development)
library(fracdiff)

.proctime00 <- proc.time()
set.seed(107)
options(digits = 10)

## 1)

x1 <-  fracdiff.sim( 5000, ar = .2, ma = -.4, d = .3, n.start=0, allow.0 = TRUE)
(fd1 <- fracdiff(x1$series, nar = 1, nma = 1, dtol = 1e-10))

fdCOVcomp <-
    c("h", "covariance.dpq", "stderror.dpq", "correlation.dpq", "hessian.dpq")
fd1. <- fracdiff.var(x1$series, fd1, h = fd1$h / 8)
fd1.[fdCOVcomp]
fd1u <- fracdiff.var(x1$series, fd1, h = fd1$h * 8)
sapply(fd1u[fdCOVcomp], signif, digits= 8)

## 2)

x2 <-  fracdiff.sim( 2048, ar = .8, ma = -.4, d = .3, n.start=0, allow.0 = TRUE)# -> NA's and problems
(fd2 <- fracdiff(x2$series, nar = length(x2$ar), nma = length(x2$ma)))

(fd2. <- fracdiff.var(x2$series, fd2, h = fd2$h / 8))[fdCOVcomp]
(fd2u <- fracdiff.var(x2$series, fd2, h = fd2$h * 8))[fdCOVcomp]

## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
