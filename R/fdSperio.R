#### by Valderio Reisen  -- Dec.2005

## MM(FIXME): This is "in parallel"  to GPH() , see ./GPH.R

Sperio <- function(x, bandw.exp = 0.5, beta = 0.9)
{
    x <- as.numeric(na.fail(as.ts(x)))
    if (any(is.na(x))) stop("NAs in x")
    if (NCOL(x) > 1) stop("only implemented for univariate time series")
    n <- length(x)
    ## Compute "smoothed" periodogram -- MM (FIXME): use  spec.pgram() !
    g <- trunc(n^bandw.exp)
    j <- 1:g
    w <- 2*pi*j/n
    mx <- mean(x)
    var.x <- sum((x - mx)^2)/n # not /(n-1)
    cov.x <- matrix(0, n-1, 1)
    for (k in 1:(n-1)) {
        cov.x[k] <- (1/n)*sum((x[1:(n-k)] - mx)*(x[(1+k):n] - mx))
    }
    M <- trunc(n^beta)
    Im <- trunc(M/2)
    pw <- numeric(n-1)
    for (k in 1:(n-1)) {
        A_k <- k/M
        pw[k] <-
            if (k <= Im)  1 - 6*A_k^2 *(1 - A_k)
            else if (k <= M) 2*(1 - A_k)^3 else 0
    }
    periodogram <- numeric(g)
    for (i in 1:g) {
        periodogram[i] <- (1/(2*pi))*(var.x + 2*sum(cov.x*pw*cos(w[i]*(1:(n-1)))))
    }
    pos <- j[periodogram > 0]
    y.reg <- log(periodogram[pos])
    x.reg <- 2*log(2*sin(w[pos]/2)) ## = log( (2*sin(..)) ^ 2)
    fit <- lm(y.reg ~ x.reg)
    d.GPH <- coef(fit)[2]
    names(d.GPH) <- NULL
    x.r2 <- sum((x.reg - mean(x.reg))^2)
    var.d <- (0.539285*M/n)/ x.r2
    var.reg <- sum((y.reg - fitted(fit))^2) / ((g - 1) * x.r2)
    list(d = -d.GPH, sd.as = sqrt(var.d), sd.reg = sqrt(var.reg))
}

