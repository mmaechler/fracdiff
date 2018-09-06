library(fracdiff)

set.seed(1)
##  examples(fdSperio)
mem.long <- fracdiff.sim(1500, d = 0.3)
spm <- fdSperio(mem.long$series)
str(spm, digits=6)

set.seed(8)
##  examples(fdGPH)
mem.l2 <- fracdiff.sim(1024, d = 0.25)
fdGPH(mem.l2$series)

diffserie0 <- fracdiff:::diffseries0 # the old slow for()-loop one
stopifnot(exprs = {
    all.equal(diffserie0(1:20, d = 1), c(-9.5, rep(1, 20-1)), tol = 1e-15)
    all.equal(diffseries(1:20, d = 1), c(-9.5, rep(1, 20-1)), tol = 1e-13) # fft
    all.equal(diffserie0(-10:10, d = 0), -10:10, tol = 1e-15)
    all.equal(diffseries(-10:10, d = 0), -10:10, tol = 1e-13)
    all.equal(diffserie0(-10:10, d = 1/2),
              diffseries(-10:10, d = 1/2), tol = 1e-13) # see 4.3e-16 on 64b-Lnx
})

set.seed(123)
## example(diffseries)
mem.l3 <- fracdiff.sim(80, d = 0.3)
mGPH <- fdGPH(mem.l3$series)
r0 <- diffserie0(mem.l3$series, d = mGPH$d)
r. <- diffseries(mem.l3$series, d = mGPH$d)
print(r0, digits = 4)
all.equal(r0, r., tol = 0) # average rel.error
stopifnot(all.equal(r0, r., tol = 1e-13))
print(acf(r0)) #
mtext("(shouldn't show structure - ideally)")

cat("Time used: ", proc.time(),"\n") # for ``statistical reasons''
