library(multitaper)
library(rbenchmark)
# refresh latest core functionality
source('~/nute.processing/development/rlpSpec/rsrc/func_psdcore.R')

initEnv(refresh=TRUE)
nd <- 1e3 # differences in processing time reduce with orders of magnitude increases
X.d <- arima.sim(list(order = c(1,1,0), ar = 0.9),n=nd)
X.p <- prewhiten(X.d, 10, plot=FALSE, verbose=FALSE)
nt <- 8
ntaps <- rep.int(nt,nd)

PSD <- .psdcore.default
# not really any improvement in speed for byte-compiled version:
#library(compiler)
#PSDc <- cmpfun(.psdcore.default)

print(benchmark(PSD(X.d, ntaper=nt, plot=FALSE, force_calc=TRUE),
                PSD(X.d, ntaper=ntaps, plot=FALSE, force_calc=TRUE),
                multitaper::spec.mtm(X.d,k=nt,plot=FALSE),
                replications=5))