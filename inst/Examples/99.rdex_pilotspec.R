#RDEX#\dontrun{
##
## Pilot spectrum
##
data(magnet)
## simply calculate the pilot spectrum with a few tapers
plot(pilot_spec(xc <-  magnet$clean), log="dB")
## remove the effect of an AR model works exceptionally
## well for the Project MAGNET data:
plot(pilot_spec(xc, remove.AR=10), log="dB", add=TRUE, col="red")
##
#RDEX#}
