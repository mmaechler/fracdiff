#### by Valderio Reisen  -- Dec.2005--
#### Tweaks by MM

## MM(FIXME): This is "in parallel"  to fdSperio() , see ./fdSperio.R

fdGPH <- function(x, bandw.exp = 0.5)
{
    if(NCOL(x) > 1) stop("only implemented for univariate time series")
    x <- as.numeric(na.fail(as.ts(x)))
    if (any(is.na(x))) stop("NAs in x")
    n <- length(x)
    ## Compute "smoothed" periodogram -- MM (FIXME): use  spec.pgram() !
    g <- trunc(n^bandw.exp)
    j <- 1:g
    kk <- 1:(n-1)
    w <- 2*pi*j/n
    x <- x - mean(x)
    var.x <- sum(x^2)/n # not /(n-1)
    cov.x <- numeric(n-1L)
    for (k in kk)
        cov.x[k] <- sum(x[1:(n-k)] * x[(1+k):n]) / n

    periodogram <- numeric(g)
    for (i in 1:g) # unscaled (will scale below)
        periodogram[i] <- var.x + 2*sum(cov.x * cos(w[i]*kk))

    pos <- j[periodogram > 0]
    y.reg <- log(periodogram[pos] / (2*pi))
    x.reg <- 2*log(2*sin(w[pos]/2)) ## = log( (2*sin(..)) ^ 2)
    fit <- lm.fit(cbind(1, x.reg), y.reg)
    d.GPH <- coef(fit)[["x.reg"]]
    x.r2 <- sum((x.reg - mean(x.reg))^2)
    var.d <- pi^2 / (6*x.r2)
    var.reg <- sum(resid(fit)^2) / ((g - 1) * x.r2)
    ## return
    list(d = -d.GPH, sd.as = sqrt(var.d), sd.reg = sqrt(var.reg))
}

