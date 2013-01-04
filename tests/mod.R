x<-1:10
mc1a <- mod(1,2)
mc2a <- mod(1+x,2)
mc1b <- 1 %% 2
mc2b <- 1 + x %% 2
mc2c <- (1 + x) %% 2
all.equal(mc1a, mc1b)
all.equal(mc2a, mc2b)
all.equal(mc2a, mc2c)

#
require(stats)
# quick power spectral density
psd <- spectrum(rnorm(1e2), plot=FALSE)
# return is class 'spec'
is.spec(psd) # TRUE
# but the underlying structure is a list
class(psd) <- "list"
is.spec(psd) # FALSE
#

##
## Spline gradient
x <- seq(0,5*pi,by=pi/64)
y <- cos(x) #**2
splineGrad(x, y, TRUE)
y <- y + rnorm(length(y), sd=.1)
# unfortunately, the presence of
# noise will affect numerical derivatives
splineGrad(x, y, TRUE)
# so change the smoothing ()
splineGrad(x, y, TRUE, spar=0.2)
splineGrad(x, y, TRUE, spar=0.6)
splineGrad(x, y, TRUE, spar=1.0)
##
