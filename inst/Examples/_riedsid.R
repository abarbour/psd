library(psd)
library(psd)

data(magnet)
dat <- magnet$clean[1:200]

pmo <- spectrum(dat, plot=FALSE)
ro <- riedsid(pmo)
print(ro)

pspectrum_basic <- function(x, initap=20, niter=5, plot=TRUE){
  
  message("Pilot spectrum (", initap, " tapers)")
  psd <- psdcore(x, ntaper=initap)
  kopt <- psd[['taper']]
  nf <- length(kopt)
  
  if (plot) plot(kopt, type='l', ylim=c(initap,2.1*initap))
  
  message("Iterative refinement of spectrum (", niter, " iterations)")
  for (iter in seq_len(niter)){
    message("\tstage ", iter)
    # find optimal tapers
    kopt <- riedsid(psd, kopt)
    print(tail(kopt))
    # update spectrum
    psd <- psdcore(x, ntaper=kopt)
    # plot
    lines(kopt, col=iter+1)
  }
  return(psd)
}

pm <- pspectrum_basic(dat)

stop()

nit <- 8
pm <- pspectrum(dat, niter=nit, plot=TRUE)
ah <- get_adapt_history()
tapstg <- ah[['stg_kopt']]

print(sapply(tapstg, tail))

tapmat <- t(plyr::ldply(tapstg))

itseq <- seq_len(nit+1)
itcol <- ((1:9 - 1) %% 6) + 1
matplot(tapmat, type='l', col=itcol)
text(itseq * 100, 300 + 0*itseq, itseq, font=2, cex=2, col=itcol)