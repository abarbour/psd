\dontrun{#REX
library(psd)

##
## Taper constraint procedures
##

data(magnet)
X <- magnet$clean

##
## spectrum
PSD <- psdcore(X, ntaper=10, refresh=TRUE)
## optimize tapers
kopt <- riedsid(PSD)
kopt.loess  <- riedsid(PSD, c.method="loess.smooth")
# the preferred function:
kopt2 <- riedsid2(PSD)
#
plot(as.tapers(kopt2), ylim =c(0, 60))
lines(as.tapers(kopt.loess), col='black')
lines(as.tapers(kopt), col='black', lwd=2)

##
## To compare all the methods at once:
demo("ctap")

}#REX
