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

diffser.0 <- fracdiff:::diffseries.0
stopifnot(diffser.0   (1:20, d = 1) == c(-9.5, rep(1, 20-1)),
## FIXME: not at all true
  all.equal(diffseries(1:20, d = 1),   c(-9.5, rep(1, 20-1))),
          diffser.0   (-10:10, d = 0) == -10:10,
  all.equal(diffseries(-10:10, d = 0), -10:10))

set.seed(123)
## example(diffseries)
mem.l3 <- fracdiff.sim(80, d = 0.3)
mGPH <- fdGPH(mem.l3$series)
r <- diffser.0(mem.l3$series, d = mGPH$d)
r. <- diffseries(mem.l3$series, d = mGPH$d)
print(r, digits = 4)
print(acf(r)) # shouldn't show structure - ideally

plot(r. ~ r); abline(0,1, col=2) ## systematic difference ... depends on 'd' only ???

cat("Time used: ", proc.time(),"\n") # for ``statistical reasons''
