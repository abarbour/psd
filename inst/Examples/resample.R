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

data(Tohoku, package='psd')
x <- subset(Tohoku, epoch=="preseismic")$areal #* 1e-6
#x <- c(x[1],x)
n <- length(x)

# detrend and demean
t. <- 1L:n - (n + 1)/2
sumt2. <- n * (n**2 - 1)/12
x <- x - mean(x) - sum(x * t.) * t./sumt2.

message("\t-->\tTohoku example:")

try(dev.off(dev.list()["RStudioGD"]))

plot(p1 <- psdcore(rnorm(n), verbose=TRUE), log='dB')

try(rm(p))
p <- psdcore(x, verbose=TRUE) #, first.last=FALSE)
#str(pl <- psd_envGet("last_psdcore_psd"))
#print(table(is.finite(pl)))
try(plot(p, log='dB'))

try(rm(p))
p <- psdcore(x, ntaper = 5, verbose=TRUE) #, first.last=FALSE)
#str(pl <- psd_envGet("last_psdcore_psd"))
#print(table(is.finite(pl)))
try(with(p, lines(freq, dB(spec), col='red')))

message("\t", psd_envGet('len_orig')/2, " ",
        psd_envGet('len_even')/2, " ",
        length(psd_envGet('fft_even_demeaned_padded'))/4, " "
)

message("\t-->\tdone.")
