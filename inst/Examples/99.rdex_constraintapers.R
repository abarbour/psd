# \dontrun{
##
## Constrain tapers
##
data(bsm)
##load("data/bsm/bsm.rda")
X <- bsm$pinyon_dat$CH0[1:5e3]
##
##
## spectrum
ntap <- 10
psd0 <- psdcore(X, ntaper=ntap, plotpsd=TRUE)
Xpsd <- psd0$spec
Xtap <- psd0$taper
#
kopt.slope <- riedsid(Xpsd, ntaper=Xtap, restrict.deriv="slope")
kopt.loess <- riedsid(Xpsd, ntaper=Xtap, restrict.deriv="loess")
kopt.super <- riedsid(Xpsd, ntaper=Xtap, restrict.deriv="friedman.super")
kopt       <- riedsid(Xpsd, ntaper=Xtap, restrict.deriv="none")
#
plot(kopt, log="y")
lines(kopt.super, col="blue", lwd=3)
lines(kopt.loess, col="dark green", lwd=3)
lines(kopt.slope, col="black")
##
##
# }
