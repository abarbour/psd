#<<eval=TRUE, eval=TRUE>>=
#@
set.seed(1234)

N <- 128
x <- rnorm(N, mean = 0, sd = 1)
xv <- var(x)
X <- fft(x)
class(X)
length(X)

Sa <- Mod(X) # Amplitude spectrum
Sp <- Arg(X) # Phase spectrum

# spectral density function |Sa|**2
XC <- Conj(X)
all.equal(Se <- Sa*Sa, Se_c <- Mod(XC * X), Se_cR <- Mod(X * XC))

Nf <- N/2 # number of frequencies
nyfreq <- seq.int(from=0,to=0.5,length.out=Nf) # nyquist frequencies
S <- 2*Se[1:Nf]/N

fsamp <- 1 # sampling freq, Hz
freq <- fsamp*nyfreq # Nyq-freq to Hz
fNyq <- fsamp/2 # nyquist or shannon(?) freq

plot(nyfreq, S, type="h")
abline(h=print(mSn<-round(mean(S))),lwd=2) # should be ~2

# crude representation of the integrated spectrum, which
# should equal the variance of the original series (xv)
S_sum <- function(S_ave, fr, fl=0) S_ave*(fr - fl)

S_sum(mSn, fNyq)/xv		#1 
S_sum(mSn, fNyq, -fNyq)/xv	#2
sum(Se)/N/N			#1
mSe <-round(mean(Se))
S_sum(mSe, fNyq)/xv		#64
S_sum(mSe, fNyq)/Nf/xv		#1

fres <- fsamp/N # or should it always be 1/ne?

## now, what if the sampling changes to 10/second?
fsamp <- 10
fNyq <- fsamp/2
freq <- fsamp * nyfreq # Nyq-freq to Hz
S_sum(mSn, fNyq)/xv		#10
S_sum(S_ave=mSn*fsamp, fNyq)/xv

