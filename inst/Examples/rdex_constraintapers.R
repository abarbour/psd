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
#
plot(kopt, log="y", ylim =c(.1, 3e2))
lines(kopt.loess, log="y", col="green")
##
##
## To compare all the methods at once:
demo("ctap")
##
#RDEX#}
