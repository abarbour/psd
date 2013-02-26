#RDEX#\dontrun{
require(psd)
##
## Taper constraint procedures
##
data(magnet)
X <- magnet$clean
##
## spectrum, then riedsid
kopt <- riedsid(PSD <- psdcore(X, ntaper=10, refresh=TRUE))
kopt.loess  <- riedsid(PSD, c.method="loess.smooth")
kopt.super  <- riedsid(PSD, c.method="friedman.smooth")
kopt.markov <- riedsid(PSD, c.method="markov.chain")
#
plot(kopt, log="y", ylim =c(.1, 3e2))
lines(kopt.super, log="y", col="red")
lines(kopt.loess, log="y", col="green")
lines(kopt.markov, log="y", col="orange")
##
##
## To compare all the methods at once:
demo("ctap")
##
#RDEX#}
