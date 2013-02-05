#' Reports whether x is an object with class 'rlpspec'
#' @param x An object to test
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @seealso \code{\link{as.rlpspec}}
#' @examples
#' is.rlpspec(1:10)
#' is.rlpspec(as.rlpspec(1:10))
is.rlpspec <- function(x) inherits(x, "rlpspec")

#' The 'rlpspec' S4 class.
#'
#' This class very nearly resembles the structure of a 'spec'
#' object; however, it has slots to store information about tapers,
#' the adaptive procedure, and any  ancillary processing (e.g.
#' prewhitening, etc) which may have been done on the original
#' series.
#'
# \section{Slots}{
#   \describe{
#     \item{\code{tapers}:}{Object of class \code{"integer"}, containing 
#     a vector of tapers.}
#   }
# }
#'
#' @name rlpspec
#' @rdname rlpspec
#' @aliases rlpspec-class
#' @exportClass rlpspec
#' @author Andrew Barbour <andy.barbour@@gmail.com>
rlpspec <- setClass("rlpspec",
                  # if slots, add __tion('taper="integer",...)
                    # replicate a 'spec' object, and then add
                    # slots for rlpspec specific information
                  representation=representation(freq="numeric",
                                                spec="numeric",
                                                coh="numeric",
                                                phase="numeric",
                                                kernel="numeric",
                                                df="numeric",
                                                bandwidth="numeric",
                                                n.used="integer",
                                                orig.n="integer",
                                                series="character",
                                                snames="character",
                                                method="character",
                                                taper="taper",
                                                pad="numeric",
                                                detrend="logical",
                                                demean="logical",
                                                # more, specific to rlpspec
                                                orig.sps="numeric", # integer?
                                                prewhiten="logical",
                                                adapt.stage="integer",
                                                taper.optim="logical",
                                                nyquist.normalized="logical"),
                    prototype = c(0,
                                  0,
                                  NULL,
                                  NULL,
                                  NULL,
                                  Inf,
                                  0,
                                  0L,
                                  0L,
                                  NULL,
                                  NULL,
                                  "Adaptive multitaper power spectral density",
                                  taper(),
                                  0,
                                  TRUE,
                                  TRUE,
                                  # <rlpspec params>
                                  1,
                                  FALSE, #prew
                                  0L,
                                  FALSE,
                                  FALSE)
                    )
# function (x, spans = NULL, kernel = NULL, taper = 0.1, pad = 0, 
#           fast = TRUE, demean = FALSE, detrend = TRUE, plot = TRUE, 
#           na.action = na.fail, ...)
# spg.out <- list(freq = freq, 
#                 spec = spec, 
#                 coh = NULL, 
#                 phase = NULL, 
#                 kernel = NULL, 
#                 df = Inf, # degrees freedom
#                 bandwidth = bandwidth, 
#                 n.used = N, 
#                 orig.n = N0, 
#                 series = series, 
#                 snames = colnames(x), 
#                 method = ifelse(TRUE, "Adaptive multitaper spectral density", "orig"), 
#                 taper = taper, 
#                 pad = 0, 
#                 detrend = TRUE, 
#                 demean = FALSE)
# class(spg.out) <- "spec"###

#' Coerce an object into an 'rlpspec' object
#'
#' @param X An object to coerce
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @seealso \code{\link{is.rlpspec}}, \code{\link{rlpspec}}
as.rlpspec <- function(X) UseMethod("as.rlpspec")
as.rlpspec.default <- function(X){
  class(X) <- "rlpspec"
  return(X)
}
as.rlpspec.spec <- function(X){
  class(X) <- "rlpspec"
  return(X)
}
as.rlpspec.rlpspec <- function(X, verbose=TRUE){
  if (verbose) message(sprintf("Object %s is already class 'rlpspec'\n",deparse(substitute(X))))
  return(X)
}

# as.spec <- function(...) UseMethod("as.spec")
# as.spec.rlpspec <- function(...){
#   NULL
# }

# as.spec
# function (x, spans = NULL, kernel = NULL, taper = 0.1, pad = 0, 
#           fast = TRUE, demean = FALSE, detrend = TRUE, plot = TRUE, 
#           na.action = na.fail, ...)
# spg.out <- list(freq = freq, 
#                 spec = spec, 
#                 coh = NULL, 
#                 phase = NULL, 
#                 kernel = NULL, 
#                 df = Inf, # degrees freedom
#                 bandwidth = bandwidth, 
#                 n.used = N, 
#                 orig.n = N0, 
#                 series = series, 
#                 snames = colnames(x), 
#                 method = ifelse(TRUE, "Adaptive multitaper spectral density", "orig"), 
#                 taper = taper, 
#                 pad = 0, 
#                 detrend = TRUE, 
#                 demean = FALSE)
# class(spg.out) <- "spec"###
###  S3 (or S4) methods for PRINTING
###

#' @title Generic methods for 'rlpspec' objects
#' @keywords methods generic
#' @name rlpspec-methods
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @aliases rlpspec
#' @rdname rlpspec-methods
#' @seealso \code{\link{is.rlpspec}}, \code{\link{as.rlpspec}}
NULL

#
# PRINT
#

#' @rdname rlpspec-methods
#' @name print
#' @S3method print summary.rlpspec
print.summary.rlpspec <- function(x, ...) {
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
#' @rdname rlpspec-methods
#' @name print
#' @S3method print rlpspec
print.rlpspec <- function(psd, ...){
  # if psd is a structure the
  # a structure has attributes, whereas a list has names of data, which may
  # have attributes: so the final psd should be a list of data with attributes
  #
  #d <- list(psd=c(1,3,4,5,2), freq=c(0,2,4,5,1), ntap=c(3,4,5,5,6), call="pspectrum"); class(d)<-"psd"
  #
  res <- cbind(psd$freq, psd$psd, psd$ntap)
  print(head(res))
}

#
# SUMMARY
#

#' @rdname rlpspec-methods
#' @name summary
#' @S3method summary rlpspec
summary.rlpspec <- function(object, ...) {
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
  class(res) <- "summary.rlpspec" 
  res
}

#
# PLOTTING
#

#' @rdname rlpspec-methods
#' @name plot
#' @S3method plot rlpspec
plot.rlpspec <- function(psd.df, logx=TRUE, xlabel=NULL, ylabel=NULL, niter=NULL, showplot=TRUE,...){
  ##
  ## Plot the results of the PSD estimation
  ##
  ## Args:  
  ##
  ## Returns:  
  ##
  ## TODO(abarbour): 
  ## [ ] convert this to a method for class 'psd'
  ## [ ] MOD: http://learnr.wordpress.com/2009/05/26/ggplot2-two-or-more-plots-sharing-the-same-legend/
  ##
  require(sfsmisc, quietly=TRUE, warn.conflicts=FALSE)  # for nice labels
  require(ggplot2, quietly=TRUE, warn.conflicts=FALSE)  # for plot engine
  
  dims <- dim.data.frame(psd.df)
  nrow <- dims[1]
  nvar <- dims[2]
  pltdf <- psd.df[2:(nrow-1),]
#   if (is.null(niter) && exists("Niter")){niter <- Niter}
  if (exists("num_iter",envir=psdenv) & is.null(niter)){niter <- envGet("num_iter")}
  # percent spectral uncertainty from  number of tapers
  # sigma^2 ~ 6 S^2 / 5 K
  #   A*(1-(seA-1)
  pltdf$sigma2 <- (pltdf$psd**2)*1.2/pltdf$ntaper
  
  ## plot setup and label breaks
  pltdf$x <- pltdf$f
  #
  if (is.null(ylabel)){ylabel <- "PSD, log10 units**2 * N * dT"}
  if (is.null(xlabel)){xlabel <- "Frequency, 1/N/dT"}
  if (logx){
    xlabel <- sprintf("%s log10",xlabel)
    pltdf$x <- log10(pltdf$x)
  }
  g <- ggplot(pltdf, aes(x=x, y=psd))#, group=src))
  atY <- 10^seq(-4, 10, by=1)
  atYL <- sfsmisc::axTexpr(2, at=atY, drop.1=TRUE)
  # plot grobs
  p <- g +
    # std err
    geom_ribbon(size=0.25, colour="black", aes(ymax=(1+sqrt(sigma2)/psd), ymin=(1-sqrt(sigma2)/psd), fill="sig")) +
    # tapers
    geom_ribbon(size=0.25, colour="black", aes(ymax=ntaper, ymin=envGet("init_tap"), fill="tap")) +
    # psd
    geom_path(size=0.5)+
    scale_x_continuous(xlabel)+ #, expand=c(0,0))+
    scale_y_log10(ylabel, breaks=atY, labels=atYL, expand=c(0,0))+
    scale_fill_discrete("",
                        breaks=c("sig","tap"), 
                        labels=c("Uncertainty: +- Sigma", "Optim. tapers: Rel. pilot spec."))+
    theme_bw()+
    opts(title=sprintf("Adaptive Sine-multitaper PSD Estimation\n%i iterations",niter),
         legend.position=c(0.83, 0.90))
  if (showplot){print(p)}
  return(invisible(p))
}
# end plot.rlpspec
##
##
##
