#RDEX#\dontrun{
##
## Using prewhiten to improve spectral estimates
##
data(magsat)
dx <- 1
mts <- ts(magsat$clean, frequency=dx)
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
plot(psd <- psdcore(mts.p, ntaper=ntap), log=ylog, lwd=2, ylim=c(-5,35))
# remove the effect of AR model
psd.ar <- psdcore(mts.par, ntaper=ntap)
psd.ar$spec <- psd.ar$spec / mean(psd.ar$spec)
psd$spec <- psd$spec / psd.ar$spec
plot(psd, log=ylog, add=TRUE, lwd=2, col="red")
plot(psd.ar, log=ylog, add=TRUE, col="blue", lwd=2)
##
#RDEX#}
