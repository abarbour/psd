\dontrun{#REX
library(psd)
##
## Multitaper PSD estimation
##
set.seed(1234)
X <- rnorm(1e3)
#
# use the defaults, and appeal to plot.spec
plot(psdcore(X))
#
# use more tapers, compare to stats::spectrum, and clear 
# env data from the previous calculation
psdcore(X, ntaper=10, plot=TRUE, refresh=TRUE)
#
# change the sampling frequency to 20
psdcore(X, 20, 10, plot=TRUE, refresh=TRUE) 
}#REX
