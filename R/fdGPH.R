GPH <- function(x, bandw.exp = 0.5)
{
    x <- as.data.frame(x)
    names(x) <- "series"
    x <- x$series
    if (NCOL(x) > 1)
        stop("only implemented for univariate time series")
    if (any(is.na(x)))
        stop("NAs in x")
    n <- length(x)
    g <- trunc(n^band)
    j <- 1:g
    w <- 2*pi*j/n
    var.x <- sum((x - mean(x))^2)/n # not /(n-1)
    cov.x <- matrix(0, n-1, 1)
    for (k in 1:(n-1)) {
        cov.x[k] <- (1/n)*sum((x[1:(n-k)] - mean(x))*(x[(1+k):n] - mean(x)))
    }
    periodogram <- matrix(0, g, 1)
    for (i in 1:g) {
        periodogram[i] <- (1/(2*pi))*(var.x + 2*sum(cov.x*cos(w[i]*(1:(n-1)))))
    }
    y.reg <- log(periodogram[periodogram > 0])
    x.reg <- log((2*sin(w/2))^2)
    pos <- j[periodogram > 0]
    x.reg <- x.reg[pos]
    fit <- lm(y.reg ~ x.reg)
    d.GPH <- coef(fit)[2]
    names(d.GPH) <- NULL
    x.r2 <- sum((x.reg - mean(x.reg))^2)
    var.d <- pi^2 / (6*x.r2)
    var.reg <- sum( (y.reg - fitted(fit))^2 ) / ((g - 1) * x.r2)
    ## return
    list(d = -d.GPH, sd.as = sqrt(var.d), sd.reg = sqrt(var.reg))
}

