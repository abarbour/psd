\dontrun{
##
## Adaptive multitaper PSD estimation
##
# load MAGSAT data
data(magsat)
#
#  use the defaults
pspectrum(X<-magsat$clean)
#
# turn off normalization
pspectrum(X, Nyquist.normalize=FALSE)
}
