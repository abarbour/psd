library(plyr)
library(psd)

data(magnet)
dat <- magnet$clean[1:300]

pmo <- spectrum(dat, plot=FALSE)
ro <- riedsid(pmo)
#print(ro)

riedsid2 <- function(psd, ntaper, ...) UseMethod("riedsid2")
riedsid2.spec <- function(psd, ntaper, ...){
  stopifnot(is.spec(psd))
  Pspec <- psd[['spec']]
  Tapseq <- psd[['freq']]
  if (missing(ntaper)){
    if (inherits(psd, 'amt')){
      psd[['taper']]
    } else {
      rep.int(1L, length(Pspec))
    }
  }
  riedsid2(Pspec, ntaper, ...)
}
riedsid2.default <- function(psd, ntaper, constrained=TRUE, verbose=TRUE, ...){
  psd <- as.vector(psd)
  ntaper <- as.vector(ntaper)
  eps <- 1e-78
  nf <- length(psd)
  nt <- length(ntaper)
  if (nt == 1) ntaper <- rep(ntaper, nf)
  nspan <- ceiling( pmin( nf/2, 7*ntaper/5 ) )
  nadd <- 1 + max(nspan)
  # Create log psd, and pad to handle begnning and end values
  ist <- nadd:2
  iend <- (nf - 1):(nf - nadd)
  S <- as.numeric(c(psd[ist], psd, psd[iend])) + eps
  Y <- log(S)
  DFUN <- function(j){
    j1 <- j - nspan[j] + nadd - 1
    j2 <- j + nspan[j] + nadd - 1
    u <- j1:j2 - (j1 + j2)/2
    L <- j2 - j1 + 1
    CC <- 12
    #
    uzero <- (L^2 - 1)/CC
    #
    # first deriv
    dY <- u  %*%  Y[j1:j2] * CC / (L*(L*L - 1))
    # second deriv
    d2Y <- (u*u - uzero)  %*%  Y[j1:j2] * 360 / (L*(L^2 - 1)*(L^2 - 4))
    #
    return(c(eps=eps, d2Y=d2Y, dYsq=dY*dY))
  }
  yders <- vapply(X=seq_len(nf), FUN=DFUN, FUN.VALUE=c(1,1,1))
  kopt <- round(480**0.2 / abs(colSums(yders))**0.4)
  if (constrained) kopt <- constrain_tapers(tapvec = kopt, verbose = verbose, ...)
  return(as.tapers(kopt))
}

koo <- riedsid2(pmo$spec, 10, constrained = FALSE)
koo.c <- ctap_simple(koo)
koo.cr <- ctap_simple_rcpp(koo)
koo.l <- ctap_loess(koo, loess.span=0.1)
print(all.equal(koo.c, koo.cr))

plot(koo, ylim=c(-5,110), main='Taper constraint methods')
lines(koo.c)
lines(koo.cr+1)
lines(koo.c - koo.cr)
lines(koo.l, col='blue')

pspectrum_basic <- function(x, initap=20, niter=5, plot=TRUE, verbose=TRUE, ...){
  
  if (verbose) adapt_message(0)
  P <- psdcore(x, ntaper=initap, preproc = FALSE, first.last=FALSE, refresh=TRUE)
  ko <- P[['taper']]
  nf <- length(ko)
  
  if (plot) plot(ko, type='l', ylim=c(initap,5.1*initap), main=paste0("Kopt\ninitial tapers: ", initap, ", iterations:", niter))
  
  # Iterate on optimal tapers, and resample spectrum
  if (verbose & niter > 0) message("Iterative refinement of spectrum (", niter, " iterations)")
  for (iter in seq_len(niter)){
    if (verbose) adapt_message(iter)
    # find optimal tapers
    ko <- riedsid2(P, ko, verbose=FALSE)
    # update spectrum
    P  <- psdcore(x, ntaper=ko, preproc = FALSE, first.last=FALSE)
    # plot
    if (plot) lines(ko, col=iter+1, lwd=2)
  }
  return(P)
}

message("\n++++>\tuniform tapered result\n")
pm0 <- pspectrum_basic(dat, initap=2, niter = 0)
message("\n++++>\tadaptively tapered result\n")
pm <- pspectrum_basic(dat)

plot(pm0)
lines(pm, col='red')

stop()

nit <- 8
pm <- pspectrum(dat, niter=nit, plot=TRUE)
ah <- get_adapt_history()
tapstg <- ah[['stg_kopt']]

print(sapply(tapstg, tail))

tapmat <- t(ldply(tapstg))

itseq <- seq_len(nit+1)
itcol <- ((1:9 - 1) %% 6) + 1
matplot(tapmat, type='l', col=itcol)
text(itseq * 100, 300 + 0*itseq, itseq, font=2, cex=2, col=itcol)