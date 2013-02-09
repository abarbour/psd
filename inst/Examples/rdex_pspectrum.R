#RDEX#\dontrun{
##
## Adaptive multitaper PSD estimation
##
# load MAGSAT data
data(magsat)
#
# magsat is 1/km, but lets change it to 10
pspectrum(X<-magsat$clean, 10)
#
# turn off normalization
pspectrum(X, 10, Nyquist.normalize=FALSE)
#RDEX#}
