library(psd)

#sourceCpp(file="src/resample_fft.cpp")

# message("\t-->\tmultiple tapers, inconsistent")
# fftz1 <- complex(real=1:2, imaginary = 1:2)
# #try(resample_fft_rcpp(fftz1, 1:2))
# 
# message("\t-->\tsingle taper, consistent:")
# fftz2 <- complex(real=1:4, imaginary = 1:4)
# #try(resample_fft_rcpp(fftz2,1:2))
# 
# message("\t-->\tsingle taper, inconsistent:")
# fftz3 <- complex(real=1:8, imaginary = 1:8)
# #try(resample_fft_rcpp(fftz3,2))

data(magnet)
data(Tohoku)

x <- subset(Tohoku, epoch=="preseismic")$areal * 1e-6
n <- length(x)

# detrend and demean
t. <- 1L:n - (n + 1)/2
sumt2. <- n * (n**2 - 1)/12
x <- x - mean(x) - sum(x * t.) * t./sumt2.

message("\t-->\tProject Magnet example:")

#try(dev.off(dev.list()["RStudioGD"]))

pm1 <- psdcore(magnet$clean, verbose=TRUE)
pm2 <- psdcore(magnet$clean, verbose=TRUE, first.last=FALSE)
pmn <- spectrum(magnet$clean, plot=FALSE)
try(with(normalize(pmn, src='spectrum'), plot(freq, dB(spec), type='l')))
try(with(pm1, lines(freq, dB(spec), col='blue')))
try(with(pm2, lines(freq, dB(spec), col='red')))

print(summary(dB(pm2$spec) - dB(pm1$spec)))

message("\t-->\tRandom-noise example:")

plot(p1 <- psdcore(rnorm(n), verbose=TRUE), log='dB')

message("\t-->\tTohoku example:")

pmn <- spectrum(x, plot=FALSE)
try(plot(normalize(pmn,src = 'spectrum'), log='dB'))

p1 <- psdcore(x, verbose=TRUE) #, first.last=FALSE)
#str(pl <- psd_envGet("last_psdcore_psd"))
#print(table(is.finite(pl)))
try(with(p1, lines(freq, dB(spec), col='blue')))

p2 <- psdcore(x, ntaper = 10, verbose=TRUE) #, first.last=FALSE)
#str(pl <- psd_envGet("last_psdcore_psd"))
#print(table(is.finite(pl)))
try(with(p2, lines(freq, dB(spec), col='red')))

message("\t-->\tdone.")
