Sperio <- function(x, band = NULL, beta = NULL)
{
    x <- as.data.frame(x)
    names(x) <- "series"
    x <- x$series
    names(x) <- NULL
    if (NCOL(x) > 1)
        stop("only implemented for univariate time series")
    if (any(is.na(x)))
        stop("NAs in x")
    if (is.null(band))
        band <- 0.5
    if (is.null(beta))
        beta <- 0.9
    n <- length(x)
    g <- trunc(n^band)
    j <- 1:g
    w <- 2*pi*j/n
    var.x <- (1/n)*sum((x - mean(x))^2)
    cov.x <- matrix(0, n-1, 1)
    for (k in 1:(n-1)){
        cov.x[k] <- (1/n)*sum((x[1:(n-k)] - mean(x))*(x[(1+k):n] - mean(x)))
    }
    M <- trunc(n^beta)
    Im <- trunc(M/2)
    pw <- matrix(0, n-1, 1)
    A <- matrix(0, n-1, 1)
    for (k in 1:(n-1)){
        A[k] <- k/M
        if (k <= Im) { pw[k] <- 1.0 - 6*(A[k]^2) + 6*(A[k]^3) }
        else if (k <= M) { pw[k] <- 2*((1 - A[k])^3) }
        else { pw[k] <- 0 }
    }
    periodogram <- matrix(0, g, 1)
    for (i in 1:g) {
        periodogram[i] <- (1/(2*pi))*(var.x + 2*sum(cov.x*pw*cos(w[i]*(1:(n-1)))))
    }
    y.reg <- log(periodogram[periodogram > 0])
    x.reg <- log((2*sin(w/2))^2)
    pos <- j[periodogram > 0]
    x.reg <- x.reg[pos]
    fit <- lm(y.reg ~ x.reg)
    d.GPH <- coef(fit)[2]
    names(d.GPH) <- NULL
    var.d <- (0.539285*M/n)/ sum((x.reg - mean(x.reg))^2)
    var.reg <- (sum((y.reg - fitted(fit))^2)/(g - 1)) /
        sum((x.reg - mean(x.reg))^2)
    list(d = -d.GPH, sd.as = sqrt(var.d), sd.reg = sqrt(var.reg))

}

