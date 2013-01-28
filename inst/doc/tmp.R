#<<eval=TRUE, eval=TRUE>>=
#@
set.seed(1234)

ne <- no <- 128
x <- rnorm(no, mean = 0, sd = 1)
xv <- var(x)
X <- fft(x)
class(X)
length(X)

Sa <- Mod(X) # Amplitude spectrum
Sp <- Arg(X) # Phase spectrum

XC <- Conj(X)
all.equal(Se <- Sa**2, Se_2 <- Mod(XC * X), Se_2R <- Mod(X * XC))

nf <- ne/2 # number of frequencies
nyfreq <- seq.int(from=0,to=0.5,length.out=nf) # nyquist frequencies
S <- Se[1:nf]

fsamp <- 1 # sampling freq, Hz
freq <- fsamp*nyfreq # Hz
fNyq <- fsamp/2 # nyquist or shannon(?) freq
plot(nyfreq, S, type="h")
plot(nyfreq, Sn<-S/no, type="h")
#lines(freq, Sn_2<-S/no/fsamp, col="blue")
abline(h=mSn<-round(mean(Sn)),lwd=2)

S_sum <- function(S_ave, fr, fl=0) S_ave*(fr - fl)

S_sum(round(mean(S)), fNyq)/nf	#1
S_sum(mSn, fNyq)/xv		#1/2 improperly normed
S_sum(mSn, fNyq, -fNyq)/xv	#1
sum(Se/no)/nf/var(x)		#2 improperly normed
sum(Se/no)/nf/var(x)*fNyq	#1

fres <- fsamp/ne # or should it always be 1/ne?
