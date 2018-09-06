#### Fractional differentiation -- the inverse of fractional integration
#### -------------------------- ----------------------------------------

## by Valderio Reisen  -- Dec.2005--
## MM:  This is 'not optimal' -- and I may have better in ../filters.R ? <<< FIXME >>>
diffseries.0 <- diffseries <- function(x, d)
{
    x <- as.data.frame(x)
    names(x) <- "series"
    x <- x$series
    if (NCOL(x) > 1)
        stop("only implemented for univariate time series")
    if (any(is.na(x)))
        stop("NAs in x")
    n <- length(x)
    stopifnot(n >= 2)
    x <- x - mean(x)
    PI <- numeric(n)
    PI[1] <- -d
    for (k in 2:n) {
        PI[k] <- PI[k-1]*(k - 1 - d)/k
    }
    ydiff <- x
    for (i in 2:n) {
        ydiff[i] <- x[i] + sum(PI[1:(i-1)]*x[(i-1):1])
    }
    ## return numeric!
    ydiff
}



## From: alexios ghalanos <alexios@4dscape.com>
## Date: Mon, 13 Jan 2014 19:58:48 +0000
## To: <maechler@stat.math.ethz.ch>
## Subject: fracdiff

## Dear Martin,

## Just a quick note, should it be of interest, that a very fast algorithm
## for diffseries was recently published:
## qed.econ.queensu.ca/working_papers/papers/qed_wp_1307.pdf‎

## Page 6 contains the R code and page 7 the benchmark timings.

## Quick check (win 7 x64, R 3.02) shows a large performance boost (for 'large' n):
#--------------------------------------
## library(microbenchmark)
## library(fracdiff)
## memory.long <- fracdiff.sim(8000, d = 0.3)
## mGPH <- fdGPH(memory.long$series)
## Jensen and Nielsen code:
## (slightly improved by MM)
diffseries <- function(x, d) {
    stopifnot((iT <- length(x)) >= 2)
    np2 <- nextn (2* iT - 1, 2)
    k <- 1:(iT -1)
    b <- c(1, cumprod ((k - d - 1)/ k))
    pad <- rep.int(0, np2-iT)
    dx <- fft(fft(c(b, pad)) *
              fft(c(x, pad)), inverse =TRUE)[1:iT]/ np2
    Re(dx)
}
## microbenchmark(diffseries(memory.long$series, d = mGPH$d),
## diffseries2(memory.long$series, d = mGPH$d))
# Unit: milliseconds
# diffseries  852.314992 (median)
# diffseries2   3.181065 (median)
#------------------------------------------------


## Best Regards,
## Alexios

