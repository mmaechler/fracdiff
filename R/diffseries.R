diffseries <- function(x, d)
{
    x <- as.data.frame(x)
    names(x) <- "series"
    x <- x$series
    if (NCOL(x) > 1)
        stop("only implemented for univariate time series")
    if (any(is.na(x)))
        stop("NAs in x")
    n <- length(x)
    x <- x - mean(x)
    PI <- matrix(0, n, 1)
    PI[1] <- -d
    for (k in 2:n){
        PI[k] <- PI[k-1]*(k - 1 - d)/k
    }
    ydiff <- NULL
    ydiff[1] <- x[1]
    for (i in 2:n) {
        ydiff[i] <- x[i] + sum(PI[1:(i-1)]*x[(i-1):1])
    }
    res <- list(sdiff = ydiff)
    res
}

