#RDEX#\dontrun{
require(psd)
##
## Using prewhiten to improve spectral estimates
##
data(magnet)
dx <- 1
mts <- ts(magnet$clean, frequency=dx)
mts.slope <- mts + seq_along(mts)
# mean + trend
# Prewhiten by removing mean+trend, and
# AR model; fit truncates the series by 
# a few terms, so zero pad
mts <- prewhiten(mts.slope,  AR.max=10, zero.pad="rear")
mts.p <- mts$prew_lm
mts.par <- mts$prew_ar
#
ntap <- 20
ylog <- "dB"
plot(PSD <- psdcore(mts.p, ntaper=ntap), log=ylog, lwd=2, ylim=c(-5,35))
# remove the effect of AR model
PSD.ar <- psdcore(mts.par, ntaper=ntap)
PSD.ar$spec <- PSD.ar$spec / mean(PSD.ar$spec)
PSD$spec <- PSD$spec / PSD.ar$spec
plot(PSD, log=ylog, add=TRUE, lwd=2, col="red")
plot(PSD.ar, log=ylog, add=TRUE, col="blue", lwd=2)
##
#RDEX#}
