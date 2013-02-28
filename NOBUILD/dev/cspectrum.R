rm(list=ls())
set.seed(1234)

# something like the methods in ifultools/src/fra_sdf.c
#
# * @param result      Pointer to a universal matrix of type
# *   MUTIL\_DCOMPLEX containing the cross-spectral density function estimates.
# *   For an M-column input matrix, the number columns in this result
# *   matrix will be M*(M+1)/2, which corresponds to all unique cross-spectral
# *   combinations that can be formed. Let Sij = conj(Xi(f)) * Xj(f)
# *   be the cross-spectral density estimate of the ith and jth time series
# *   where Xk(f) is the DFT of the kth time series. If M = 3,
# *   the result will contain the following columns: [S00 S01 S02 S11 S12 S22].
# *   The memory for this matrix is automatically allocated by the function.
# *
# * @see frauniv_spectral_density_function_direct
#
x <- y <- ts(rnorm(50))
X <- cbind(x/10,10*y)
nser <- ncol(X)
sser <- seq_len(nser)
# all possible combinations (max two columns)
allcomb <- expand.grid(sser, sser)
# columnwise products
cprods <- apply(allcomb, 1, prod)
ucprods <- unique(cprods)
prodseq <- seq_len(length(ucprods))
# get non-redundant rows (from the perspective of indices) and swap columns
ijmat <- allcomb[match(ucprods, cprods),2:1]
#
# DFT
Xf <- unclass(mvfft(X))
SFUN <- function(ccomb){
  print(ccomb)
  mn <- as.numeric(ijmat[ccomb,])
  print(mn)
  m <- mn[1]
  n <- mn[2]
  Xi <- Xf[,m]
  Xj <- Xf[,n]
  print(class(Xi))
  Sij <- cbind(Conj(Xi) * Xj)
  print(class(Sij))
  print(colnames(Sij))
  colnames(Sij) <- paste0("S",m,n,".")
  print(colnames(Sij))
  return(Sij)
}
#
Sall <- sapply(X=prodseq, FUN=SFUN, simplify=TRUE, USE.NAMES=FALSE)
head(Sall)
head(simplify2array(Sall,F))

par(mfcol=c(2,1))
matplot(20*log10(Mod(Sall)),type="l")
matplot((Arg(Sall)*180/pi),type="l", ylim=180*c(-1,1))

#
# spec.pgram
#
#  N <- nrow(x)
#  Nspec <- floor(N/2)
#  x <- rbind(x, matrix(0, nrow = (NewN - N), ncol = ncol(x)))
#  xfft <- mvfft(x)
#  pgram <- array(NA, dim = c(N, ncol(x), ncol(x)))
#  for (i in 1L:ncol(x)) {
#    for (j in 1L:ncol(x)) {
#      pgram[, i, j] <- xfft[, i] * Conj(xfft[, j])/(N0 * xfreq)
#      pgram[1, i, j] <- 0.5 * (pgram[2, i, j] + pgram[N, i, j])
#    }
#  }
#  spec <- matrix(NA, nrow = Nspec, ncol = nser)
#  for (i in 1L:nser) spec[, i] <- Re(pgram[1L:Nspec, i, i])
#  coh <- phase <- matrix(NA, nrow = Nspec, ncol = nser * (nser - 1)/2)
#  for (i in 1L:(nser - 1)) {
#    for (j in (i + 1):nser) {
#      coh[, i + (j - 1) * (j - 2)/2] <- Mod(pgram[, i, j])^2/(spec[, i] * spec[, j])
#      phase[, i + (j - 1) * (j - 2)/2] <- Arg(pgram[, i, j])
#    }
#  }
#  for (i in 1L:nser) spec[, i] <- spec[, i]/u2
#  spec <- drop(spec)
