library(fftw)
library(rbenchmark)
library(plyr)

dftbm <- function(nd, reps=3e2){
	x <- rnorm(1e4)
	bmd <- benchmark(replications=reps, fftw::FFT(x), stats::fft(x),
			 columns=c('test', 'elapsed', 'relative', 'replications'))
	bmd$log10_num_dat <- log10(nd)
	return(bmd)
}

benchdat <- plyr::ldply(lapply(X=10**c(1:6), FUN=dftbm))

## stats::fft is the clear winner for real, univariate series:
benchdat[order(benchdat$relative),]
