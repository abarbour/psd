library(multitaper)
library(rbenchmark)
library(compiler)
source('~/nute.processing/development/rlpSpec/rsrc/dev_func_psdcore.R')

initEnv(refresh=TRUE)
nd <- 1e3 # differences in processing time reduce with orders of magnitude increases
X.d <- arima.sim(list(order = c(1,1,0), ar = 0.9),n=nd)
X.p <- prewhiten(X.d, 10, plot=FALSE, verbose=FALSE)
nt <- 8
ntaps <- rep.int(nt,nd)

PSD <- ..dev_psdcore.default
#PSDc <- cmpfun(..dev_psdcore.default)

print(benchmark(PSD(X.d, ntaper=nt, plot=FALSE, force_calc=TRUE),
                #PSDc(X.d, ntaper=nt, plot=FALSE, force_calc=TRUE),
                multitaper::spec.mtm(X.d,k=nt,plot=FALSE),
                replications=5))