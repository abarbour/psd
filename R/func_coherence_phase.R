#' coherence
#'
#' calculate coherency from the spectra and cross-spectra
#'
#' @param pgram \code{numeric array} must be multivariate
#'
#'
#' @return list of coherency
#'
#' @export
#'
#'
coherence <- function(pgram) {

    dims = dim(pgram)
    n_ser <- dims[2]
    n_r <- dims[1]
    
    if (n_ser == 1) {
      warning('Cannot calculate coherency for a single periodogram')
      return(NA)
    }
    
    coh <- matrix(NA_real_, nrow = n_r, ncol = (n_ser - 1))
    
    for (i in 2L:(n_ser)) {
      coh[, i-1]   <- Re(Mod(pgram[, 1, i])^2 / (pgram[, 1, 1] * pgram[, i, i]))
    }
    
    return(coh)
  }

#' phase
#'
#' calculate phase from the spectra and cross-spectra
#'
#' @param pgram \code{numeric array} must be multivariate
#'
#'
#' @return list of phase
#'
#' @export
#'
phase <- function(pgram) {
  
  dims = dim(pgram)
  n_ser <- dims[2]
  n_r <- dims[1]
  
  if (n_ser == 1) {
    warning('Cannot calculate phase for a single periodogram')
    return(NA)
  }
  
  phase <- matrix(NA_real_, nrow = n_r, ncol = (n_ser - 1))
  
  for (i in 2L:(n_ser)) {
  
    phase[, i-1] <- Arg(pgram[, 1, i])
 
  }
  
  return(phase)
}