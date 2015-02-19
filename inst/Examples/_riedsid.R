library(psd)
library(psd)

data(magnet)
dat <- magnet$clean[1:200]

pmo <- spectrum(dat, plot=FALSE)
ro <- riedsid(pmo)
print(ro)

pspectrum_basic <- function(x, initap=20, n.iter=5, plot=TRUE){
  
  message("Pilot spectrum (", initap, " tapers)")
  psd <- psdcore(x, ntaper=initap)
  if (plot) plot(psd)
  kopt <- psd[['taper']]
  nf <- length(kopt)
  
  message("Iterative refinement of spectrum (", n.iter, " iterations)")
  for (iter in seq_len(n.iter)){
    message("\tstage ", iter)
    # find optimal tapers
    kopt <- riedsid(psd, kopt)
    print(tail(kopt))
    # update spectrum
    psd <- psdcore(x, ntaper=kopt)
    if (plot) plot(psd, col=iter+1, add=TRUE)
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