#' coherence
#'
#' Calculate coherence from the spectra and cross-spectra. This method is the same as
#' used in \code{spec.pgram}.
#'
#' @param pgram \code{numeric array} must be multivariate
#'
#'
#' @return list of coherence. For multivariate time series, a 
#' matrix containing the squared coherency between different series. 
#' Column i + (j - 1) * (j - 2)/2 contains the squared coherency 
#' between columns i and j of x, where i < j.
#'
#' @export
#'
# [ ] add methods for spec, ...
coherence <- function(pgram) {
  
  dims = dim(pgram)
  n_ser <- dims[2]
  n_r <- dims[1]
  
  if (n_ser == 1) {
    warning('Cannot calculate coherency for a single periodogram')
    return(NA)
  }
  
  coh <- matrix(NA_real_, nrow = n_r, ncol = n_ser * (n_ser - 1L) / 2L)
  
  for (i in 1L:(n_ser-1)) {
    
    for (j in (i + 1):n_ser) {
      
      ind <- i + (j - 1) * (j - 2) / 2
      
      coh[, ind]   <- Re(Mod(pgram[, i, j])^2 / 
                           (pgram[, i, i] * pgram[, j, j]))
      
    }
  }
  
  return(coh)
}


#' phase
#'
#' Calculate phase from the spectra and cross spectrum. This method is the same as
#' used in \code{spec.pgram}.
#'
#' @param pgram \code{numeric array} must be multivariate
#'
#'
#' @return list of phase. For multivariate time series a matrix containing the cross 
#' spectrum phase between different series. The format is the same as \code{\link{coherence}}.
#'
#' @export
#'
# [ ] add methods for spec, etc
phase <- function(pgram) {
  
  dims = dim(pgram)
  n_ser <- dims[2]
  n_r <- dims[1]
  
  if (n_ser == 1) {
    warning('Cannot calculate phase for a single periodogram')
    return(NA)
  }
  
  
  phase <- matrix(NA_real_, nrow = n_r, ncol = n_ser * (n_ser - 1L) / 2L)
  
  for (i in 1L:(n_ser-1)) {
    
    for (j in (i + 1):n_ser) {
      
      ind <- i + (j - 1) * (j - 2) / 2
      
      phase[, ind] <- Arg(pgram[, i, j])
      
    }
  }
  
  return(phase)
}