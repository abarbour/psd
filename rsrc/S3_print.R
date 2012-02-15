###
###  S3 (or S4) methods for PRINTING
###
print.summary.psd <- function(x, ...) {
  ##
  ## Printing of a summary.psd class
  ##
  ## Args:  
  ##
  ## Returns:  
  ##
  ## TODO(abarbour):
  ##
  cat("\n>>>> ADAPTIVE SINE-MULTITAPER PSD SUMMARY <<<<\n")
  # [ ] variance reduction? 
  # [ ] number of iterations
  cat("\tCall:\t", x$call[1]) 
  cat("\n\tNumber of frequencies:\t", x$nfreq[1])
  cat("\n\tQuantiles:\n")
  freq <- x$freq
  freq_spacing <- x$df
  psd <- x$psd
  ntap <- x$ntap
  print(rbind(freq, freq_spacing, psd, ntap))
  cat("\n")
}
print.psd <- function(psd, ...){
  ##
  ## Print a 'psd' class object
  ##
  ## Args:  
  ##
  ## Returns:	
  ##
  ## TODO(abarbour):
  ##
  # if psd is a structure the
  # a structure has attributes, whereas a list has names of data, which may
  # have attributes: so the final psd should be a list of data with attributes
  res <- cbind(psd$freq, psd$psd, psd$ntap)
  print(head(res))
}
###