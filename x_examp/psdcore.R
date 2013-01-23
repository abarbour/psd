\dontrun{
X.d <- rnorm(1e3)
plot(psdcore(X.d, ntaper=10), log="dB", ylim=10*c(-1,1))
psd.n <- psdcore(X.d, ntaper=10, Nyquist.normalize=FALSE)
lines(psd.n$freq, 10*log10(psd.n$spec), col="red") # note normalization
abline(h=c(0, 3), col=c("black","red"), lwd=2)
#
# 10Hz sampling
plot(psdcore(X.d, X.frq=10, ntaper=10), log="dB", ylim=10*c(-0.3,1.7))
psd.n <- psdcore(X.d, X.frq=10, ntaper=10, Nyquist.normalize=FALSE)
lines(10*psd.n$freq, 10*log10(psd.n$spec), col="red") # note normalization
abline(h=c(10, 3), col=c("black","red"), lwd=2)
#
# if ntaper is a vector:
psdcore(X.d, ntaper=rep(10,length(X.d))
}
