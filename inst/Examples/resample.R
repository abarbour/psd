library(Rcpp)

#sourceCpp(file="src/resample_fft.cpp")

message("\t-->\tmultiple tapers, inconsistent")
fftz1 <- complex(real=1:2, imaginary = 1:2)
#try(resample_fft_rcpp(fftz1, 1:2))

message("\t-->\tsingle taper, consistent:")
fftz2 <- complex(real=1:4, imaginary = 1:4)
#try(resample_fft_rcpp(fftz2,1:2))

message("\t-->\tsingle taper, inconsistent:")
fftz3 <- complex(real=1:8, imaginary = 1:8)
#try(resample_fft_rcpp(fftz3,2))

data(Tohoku, package='psd')
x <- subset(Tohoku, epoch=="preseismic")$areal * 1e6
#x <- c(x[1],x)
n <- length(x)

message("\t-->\tTohoku example:")

plot(p <- psdcore(x))

message(psd_envGet('len_orig')/2, " ",
        psd_envGet('len_even')/2, " ",
        length(psd_envGet('fft_even_demeaned_padded'))/4, " "
)

message("\t-->\tdone.")
