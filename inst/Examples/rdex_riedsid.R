# \dontrun{
##
## Riedel-Sidorenko--Parker taper optimization
##
set.seed(1234)
# some params
nd <- 1e3 # num data
ntap <- 10 # num tapers
nrm <- 40 # sharpness of the peaks rel 2*variance
#
# create a pseudo spectrum
riex <- rnorm(nd) + nrm*abs(cos(pi*(0:(nd-1))/180) + 1.2)
#
# optimize tapers
rtap <- riedsid(riex, ntaper=ntap)
# ... again, but do not constrain
rtap2 <- riedsid(riex, ntaper=ntap, constrained=FALSE)
#
# plot
op <- par(no.readonly = TRUE)
par(mfrow=c(2,1), mar=rep(1.3,4), mai=rep(0.6,4))
# ... the mock spectrum
plot(riex, type="h", xaxs="i") 
# ... the optimal tapers
# unconstrained: note logarithmic peaks from noise in spectrum
# even when variance is 20 times smaller than peak amplitudes
plot(rtap2, log="y") 
# constrained (red curve): much better behaved
lines(rtap) 
# original tapers:
lines(as.taper(rep.int(ntap,nd)), col="blue")
par(op)
# }
