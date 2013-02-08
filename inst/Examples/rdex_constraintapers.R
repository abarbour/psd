#RDEX#\dontrun{
##
## Taper constraint procedures
##
data(magsat)
X <- magsat$clean
##
## spectrum, then riedsid
kopt <- riedsid(psd <- psdcore(X, ntaper=10, refresh=TRUE))
kopt.loess  <- riedsid(psd, c.method="loess.smooth")
kopt.super  <- riedsid(psd, c.method="friedman.smooth")
kopt.markov <- riedsid(psd, c.method="markov.chain")
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
