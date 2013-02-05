###
###  S3 (or S4) methods for SUMMARIES
###
summary.psd <- function(object, ...) {
  ##
  ##  Form a summary of a 'psd' class object
  ##
  ## Args:  
  ##
  ## Returns:  
  ##
  ## TODO(abarbour):
  ##
  # from doc example
  #   se <- sqrt(diag(object$vcov)) 
  #   tval <- coef(object) / se
  #   TAB <- cbind(Estimate = coef(object), StdErr = se,
  #                t.value = tval, p.value = 2*pt(-abs(tval), df=object$df))
  #   res <- list(call=object$call, coefficients=TAB)
  #
  res <- list(call=object$call,
              df = quantile(diff(sort(d$freq))),
              nfreq = length(object$freq),
              rfreq = range(object$freq),
              freq = quantile(object$freq),
              psd = quantile(object$psd),
              ntap = quantile(object$ntap))
  class(res) <- "summary.psd" 
  res
}
###