#RDEX#\dontrun{
require(psd)
##
## Riedel-Sidorenko--Parker taper optimization
##
set.seed(1234)
# some params
nd <- 512 # num data
ntap <- 10 # num tapers
nrm <- 40 # sharpness of the peaks rel 2*variance
#
# create a pseudo spectrum
# with broad peaks
riex <- rnorm(nd) + nrm*abs(cos(pi*(x<-0:(nd-1))/180) + 1.2)
riex <- riex + 8*nrm*dcauchy(x, nd/3)
riex <- riex + 5*nrm*dnorm(x, nd/2)
# flat regions
riex[riex<25] <- 25
ried <- dB(riex, invert=TRUE)
#
# optimize tapers
rtap <- riedsid(riex, ntaper=ntap)
#
# plot
op <- par(no.readonly = TRUE)
par(mfrow=c(2,1), mar=rep(1.3,4), mai=rep(0.6,4))
# ... the mock spectrum
plot(riex, type="h", xaxs="i", ylim=c(0,200)) 
# ... the optimal tapers
plot(rtap, log="y") 
# original tapers:
lines(as.tapers(rep.int(ntap,nd)), col="blue")
par(op)
#RDEX#}
