#RDEX#\dontrun{
require(psd)
##
## Normalization
##
set.seed(1234)
# timeseries with sampling frequency **not** equal to 1:
X <- ts(rnorm(1e3), frequency=20)
# spec.pgram: double sided
pgram <- spectrum(X)
# psdcore: single sided
PSD <- psdcore(X)
# note the normalization differences:
plot(pgram, log="dB", ylim=c(-40,10))
plot(PSD, add=TRUE, col="red", log="dB")
# A crude representation of integrated spectrum: 
#   should equal variance of white noise series (~= 1)
mean(pgram$spec)*max(pgram$freq)
mean(PSD$spec)*max(PSD$freq)
#
# normalize objects with class 'spec'
pgram <- normalize(pgram, src="spectrum")
PSD <- normalize(pgram, src="psd")
# replot them
plot(pgram, log="dB", ylim=c(-40,10))
plot(PSD, add=TRUE, col="red", log="dB")
#
# Again, integrated spectrum should be ~= 1:
mean(pgram$spec)*max(pgram$freq)
mean(PSD$spec)*max(PSD$freq)
#
#RDEX#}
