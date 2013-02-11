#RDEX#\dontrun{
##
## Pilot spectrum
##
data(magsat)
## simply calculate the pilot spectrum with a few tapers
plot(pilot_spec(xc <-  magsat$clean), log="dB")
## remove the effect of an AR model works exceptionally
## well for the MAGSAT data:
plot(pilot_spec(xc, remove.AR=10), log="dB", add=TRUE, col="red")
##
#RDEX#}
