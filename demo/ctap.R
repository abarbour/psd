##
## Show usage of taper constaint procedures
##
# fake taper series
set.seed(1234)
nd <- 2e2
x <- rnorm(nd, sd=10) + c(1:(nd/2),(nd/2 + 1):2)
ntap <- as.tapers(rep(x,3))

# Raw tapers
op <- par(no.readonly = TRUE)
par(mfrow=c(3,2))
ylim <- c(0,200)
plot(ntap, ylim=ylim, main="Tapers: Raw")

## Constraints
## minspan
plot(minspan(ntap), ylim=ylim, main="minspan")

## friedman
plot(ntap, ylim=ylim, main="friedman")
lines(ctap_friedman(ntap), col="green", lwd=2) # too flat: span is too wide
lines(ctap_friedman(ntap, smoo.span=.1), col="red", lwd=2)

## loess
plot(ntap, ylim=ylim, main="loess") # OK, but a bit too smooth: reduce span
lines(ctap_loess(ntap), col="green", lwd=2)
lines(ctap_loess(ntap, loess.span=0.1), col="red", lwd=2)

## markov
plot(ntap, ylim=ylim, main="markov")
lines(ctap_markov(ntap), col="green", lwd=2) # unstable and poor: chain.width is too long
lines(ctap_markov(ntap, chain.width=40), col="red", lwd=2) # note peak enhancement

## first difference
plot(ntap, ylim=ylim, main="simple")
lines(ctap_simple(ntap, maxslope=5), col="green", lwd=1) #resembles noise: use min slope (1, the default)
lines(ctap_simple(ntap), col="red", lwd=2)
#lines(minspan(ntaps), col="red", lwd=2)
par(op)
