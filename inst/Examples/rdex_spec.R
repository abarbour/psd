\dontrun{#REX
library(psd)
##
## Objects with class 'spec'
##
set.seed(1234)
#
x <- spectrum(xn<-rnorm(10), plot=FALSE)
xdf <-as.data.frame(x)
str(xdf)
is.tapers(xdf$taper)
#
# tapers class is retained
#
x <- psdcore(xn)
xdf <- as.data.frame(x)
str(xdf)
is.tapers(xdf$taper)
}#REX
