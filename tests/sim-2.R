require(fracdiff)

.ptime <- proc.time()
##>>  *no* *.Rout.save here  ===> can well have if(doExtras) ...
(doExtras <- fracdiff:::doExtras())

## confirm that we guessed right:
## fracdiff.sim(....., d = 0, backComp = FALSE)  <===>  arima.sim(....)

AR <- c(0.7, -0.1, 0.2)
MA <- c(-0.5, 0.4, 0.4)
n <- 512 ; sd <- 0.1
n.st <- 10

set.seed(1)
for(i in 1:200) {
    cat(sprintf("%3d ", i))
    p <- sample(0:length(AR), 1)
    q <- sample(0:length(MA), 1)
    .ar <- AR[seq_len(p)]
    .ma <- MA[seq_len(q)]
    n.st <- p+q+ rpois(1, lambda = 2)
    sid <- round(runif(1)* 1000)
    set.seed(sid)
    y1 <- arima.sim(list(ar = .ar, ma = .ma), n = n, n.start = n.st, sd = sd)
    set.seed(sid)
    y2 <- fracdiff.sim(n = n, ar = .ar, ma = - .ma, d = 0,
                       n.start = n.st, sd = sd, backComp = FALSE)
    if(!isTRUE(aeq <- all.equal(c(y1), y2$series, tol = 1e-15)))
        cat("y1 and y2 are not equal: ", aeq,"\n")
    if(!(i %% 10))    cat("\n")
}

cat('Time elapsed: ', proc.time() - .ptime,'\n'); .ptime <- proc.time()

if(doExtras) withAutoprint({
    isee <- 1:500
    do5c <- lapply(isee, function(seed) { set.seed(seed)
        tst <- fracdiff.sim(1000, ar = .6, ma = -.4, d = .3)$ser
        fracdiff(tst, nar = 1, nma = 1)
    })
    ## There were 56 warnings (use warnings() to see them)
    summary(warnings()) ## 56 x  { unable to compute correlation matrix; maybe change 'h' }
    str( msg5c <- sapply(do5c, `[[`, "msg") ) # matrix with 2 rows
    table(msg5c["fdcov", ])
    str(which(prbl <- msg5c["fdcov", ] != "ok"))
    ## int [1:56] 5 22 26 29 37 39 65 66 69 71 ...

    ## now want to use fracdiff.var() with different  'h' to *fix* the failures
    ## (need to re-simulate as we did not keep 'tst's above)
    L <- lapply(isee[prbl], function(seed) {
        set.seed(seed)
        tst <- fracdiff.sim(1000, ar = .6, ma = -.4, d = .3)$ser # -> 56 warnings
        fd. <- fracdiff(tst, nar = 1, nma = 1)
        hs <- fd.[["h"]] / 2^(3:19)
        fd2L <- lapply(hs, function(h) suppressWarnings(fracdiff.var(tst, fd., h = h)))
        table(msgCov <- sapply(fd2L, `[[`, "msg")["fdcov",])
        if(any(ok <- msgCov == "ok")) list(h.old = fd.[["h"]], fits.ok = fd2L[ok]) else list()
    })
    cat("fracdiff.var() fixed:", sum(okL <- lengths(L) > 0), "(out of", length(L), "problem cases)\n")
    ## MM: see  fracdiff.var() fixed: 43 (out of 56 problem cases) -- sometimes really needing  h/32768
    invisible(lapply(L[okL], function(lst) {
        cat(sprintf("h.old= %12g", lst[["h.old"]]))
        if(FALSE) {## too much:
            ## cat("h.old: ", lst[["h.old"]], " summary(<new ones>):\n")
            lapply(lst[["fits.ok"]], function(fm) print(summary(fm)))
        }
        cat(" -- h ratio(s) old/new = 2 ^ ... : ")
        print(log2(lst[["h.old"]] / sapply(lst[["fits.ok"]], `[[`, "h")))
    }))
    ##
    cat('Time elapsed: ', proc.time() - .ptime,'\n')
})
