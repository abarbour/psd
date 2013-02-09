#RDEX#\dontrun{
##
## Using prewhiten to improve spectral estimates
##
data(magsat)
dt <- 1
mts <- ts(magsat$clean, deltat=dt)
mts.slope <- mts + seq_along(mts)
# mean + trend
mts.p <- prewhiten(mts.slope)
l0 <- length(mts.p)
# AR model
mts.par <- prewhiten(mts,  AR.max=10)
l2 <- length(mts.par)
# AR fit truncates the series by a few terms, so zero pad
mts.par <- c(mts.par, zeros(abs(l2-l0)))
#
ntap <- 20
ylog <- "dB"
plot(psd <- psdcore(mts.p, ntaper=ntap), log=ylog, lwd=2, ylim=c(-5,35))
# remove the effect of AR model
psd.ar <- psdcore(mts.par, ntaper=ntap)
psd.ar$spec <- psd.ar$spec / mean(psd.ar$spec)
psd$spec <- psd$spec / psd.ar$spec
plot(psd, log=ylog, add=TRUE, lwd=2, col="red")
plot(psd.ar, log=ylog, add=TRUE, col="blue", lwd=2)
##
#RDEX#}
