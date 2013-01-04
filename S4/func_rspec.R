#' Reports whether x is an object with class 'rspec'
#' @param x An object to test
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @seealso \code{\link{as.rspec}}
#' @examples
#' is.rspec(1:10)
#' is.rspec(rlp <- newRspec())
#' is(1:10, "rspec")
#' is(rlp, "rspec")
is.rspec <- function(x) inherits(x, "rspec")
#is.rspec <- function(x) is(x, 'rspec')

#' The 'rspec' S4 class.
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
#' @name newRspec
#' @rdname rspec
#' @aliases rspec-class, rspec
#' @exportClass rspec
#' @export
#' @author Andrew Barbour <andy.barbour@@gmail.com>
newRspec <- setClass("rspec",
     representation=representation(freq="numeric",
                                   spec="numeric",
                                   coh="numeric",
                                   phase="numeric",
                                   kernel="numeric",
                                   df="numeric",
                                   bandwidth="numeric",
                                   n.used="numeric",
                                   orig.n="numeric",
                                   series="character",
                                   snames="character",
                                   method="character",
                                   taper="taper",
                                   pad="numeric",
                                   detrend="logical",
                                   demean="logical",
                                   # more, specific to rspec
                                   orig.sps="numeric", # integer?
                                   prewhiten="logical",
                                   adapt.stage="integer",
                                   taper.optim="logical",
                                   nyquist.normalized="logical"),
     prototype = prototype(coh=numeric(0),
                           phase=numeric(0),
                           kernel=numeric(0),
                           df=Inf,
                           snames="",
                           method = "Adaptive multitaper spectral density",
                           taper=as.taper(0),
                           detrend=FALSE,
                           demean=FALSE,
                           # <rspec params>
                           orig.sps=1,
                           prewhiten=FALSE,
                           adapt.stage=integer(0),
                           taper.optim=FALSE,
                           nyquist.normalized=FALSE),
     validity = function(object) {
       length(object@freq) == length(object@spec)
     }
)


nonull <- function(psd) UseMethod("nonull")
nonull.spec <- function(psd){
  stopifnot(inherits(psd, 'spec', FALSE))
  # spec.pgram may return NULL for these:
  psd$coh <- as.numeric(psd$coh)
  psd$phase <- as.numeric(psd$phase)
  psd$kernel <- as.numeric(psd$kernel)
  psd$snames <- as.character(psd$snames)
  psd$taper <- as.taper(psd$taper)
  return(psd)
}

#' upconvert an S3 class 'spec' to the S4 class 'rspec'
# @S3method as.rspec spec
as.rspec.spec <- function(psd){
  stopifnot(inherits(psd, 'spec', FALSE))
  psd <- nonull.spec(psd)
  S4spec <- newRspec() # or: S4spec <- new("rspec")
  spec_slots <- slotNames(S4spec)
  spec_slots <- spec_slots[match(names(psd), spec_slots, nomatch=0)]
  for (what in spec_slots){
    #message(what)
    if ((what=="taper")==TRUE){
      FUN <- as.taper
    } else {
      FUN <- as.vector
    }
    #message(what)
    slot(S4spec, what) <- FUN(unlist(psd[what]))
  }
  return(S4spec)
}

###
###  S3 (or S4) methods for PRINTING
###

#' @title Generic methods for 'rspec' objects
#' @keywords methods generic
#' @name rspec-methods
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' @aliases rspec
#' @rdname rspec-methods
#' @seealso \code{\link{is.rspec}}, \code{\link{as.rspec}}
NULL

#
# PRINT
#

#' @rdname rspec-methods
#' @name print
#' @S3method print summary.rspec
print.summary.rspec <- function(x, ...) {
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
#' @rdname rspec-methods
#' @name print
#' @S3method print rspec
print.rspec <- function(psd, ...){
  # if psd is a structure the
  # a structure has attributes, whereas a list has names of data, which may
  # have attributes: so the final psd should be a list of data with attributes
  #
  #d <- list(psd=c(1,3,4,5,2), freq=c(0,2,4,5,1), ntap=c(3,4,5,5,6), call="pspectrum"); class(d)<-"psd"
  #
  res <- cbind(psd$freq, psd$psd, psd$ntap)
  print(head(res))
}

#' @rdname rspec-methods
#' @name print
#' @exportMethod print rspec
#' @docType methods
setMethod("print", signature("rspec"), function(x){str(x)})

#' @rdname rspec-methods
#' @name summary
#' @exportMethod summary rspec
#' @docType methods
setMethod("summary", signature("rspec"), function(object){
  frqs <- object@freq
  psds <- object@spec
  ntap <- object@taper
  res <- rbind(
              freq.spacing = quantile(diff(frqs)),
              freq = quantile(frqs),
              psd = quantile(psds),
              ntap = quantile(ntap)
              )
  return(res)
  }
)

#
# PLOTTING
#

# # @rdname rspec-methods
# # @name plot
# # @S3method plot rspec
# plot.rspec <- function(psd.df, logx=TRUE, xlabel=NULL, ylabel=NULL, niter=NULL, showplot=TRUE,...){
#   ##
#   ## Plot the results of the PSD estimation
#   ##
#   ## Args:  
#   ##
#   ## Returns:  
#   ##
#   ## TODO(abarbour): 
#   ## [ ] convert this to a method for class 'psd'
#   ## [ ] MOD: http://learnr.wordpress.com/2009/05/26/ggplot2-two-or-more-plots-sharing-the-same-legend/
#   ##
#   require(sfsmisc, quietly=TRUE, warn.conflicts=FALSE)  # for nice labels
#   require(ggplot2, quietly=TRUE, warn.conflicts=FALSE)  # for plot engine
#   
#   dims <- dim.data.frame(psd.df)
#   nrow <- dims[1]
#   nvar <- dims[2]
#   pltdf <- psd.df[2:(nrow-1),]
# #   if (is.null(niter) && exists("Niter")){niter <- Niter}
#   if (exists("num_iter",envir=psdenv) & is.null(niter)){niter <- envGet("num_iter")}
#   # percent spectral uncertainty from  number of tapers
#   # sigma^2 ~ 6 S^2 / 5 K
#   #   A*(1-(seA-1)
#   pltdf$sigma2 <- (pltdf$psd**2)*1.2/pltdf$ntaper
#   
#   ## plot setup and label breaks
#   pltdf$x <- pltdf$f
#   #
#   if (is.null(ylabel)){ylabel <- "PSD, log10 units**2 * N * dT"}
#   if (is.null(xlabel)){xlabel <- "Frequency, 1/N/dT"}
#   if (logx){
#     xlabel <- sprintf("%s log10",xlabel)
#     pltdf$x <- log10(pltdf$x)
#   }
#   g <- ggplot(pltdf, aes(x=x, y=psd))#, group=src))
#   atY <- 10^seq(-4, 10, by=1)
#   atYL <- sfsmisc::axTexpr(2, at=atY, drop.1=TRUE)
#   # plot grobs
#   p <- g +
#     # std err
#     geom_ribbon(size=0.25, colour="black", aes(ymax=(1+sqrt(sigma2)/psd), ymin=(1-sqrt(sigma2)/psd), fill="sig")) +
#     # tapers
#     geom_ribbon(size=0.25, colour="black", aes(ymax=ntaper, ymin=envGet("init_tap"), fill="tap")) +
#     # psd
#     geom_path(size=0.5)+
#     scale_x_continuous(xlabel)+ #, expand=c(0,0))+
#     scale_y_log10(ylabel, breaks=atY, labels=atYL, expand=c(0,0))+
#     scale_fill_discrete("",
#                         breaks=c("sig","tap"), 
#                         labels=c("Uncertainty: +- Sigma", "Optim. tapers: Rel. pilot spec."))+
#     theme_bw()+
#     opts(title=sprintf("Adaptive Sine-multitaper PSD Estimation\n%i iterations",niter),
#          legend.position=c(0.83, 0.90))
#   if (showplot){print(p)}
#   return(invisible(p))
# }
# # end plot.rspec

# setGeneric("plot", function(x, y, ...) standardGeneric("plot")) 
# 
#' @rdname rspec-methods
#' @name plot
#' @exportMethod plot rspec
#' @docType methods
#thank god:
#  http://www.soph.uab.edu/Statgenetics/Events/Rshort/060227-8-s4slides.pdf
setGeneric("plot", function (x, y, ...) standardGeneric("plot"))
setMethod("plot",c("rspec","missing"),function(x, y, 
                                                 type="b",
                                                 null.line=TRUE,
                                                 xlab="frequency",
                                                 ylab="power spectral density",
                                                 y.dB=TRUE,
                                                 main=NULL,...){
  par(pty="s")
  frq <- x@freq
  spc <- x@spec
  if (y.dB) spc <- 10*log10(spc)
  if (length(frq)>1e3){message("too many values to see points effectively: type=line"); type <- "l"}
  if (x@nyquist.normalized){
    xlabp <- "Hz"
  } else {
    xlabp <- "Nyq"
  }
  xlab <- sprintf("%s [%s]", xlab, xlabp)
  ylab <- sprintf("%s [counts^2 / %s]", ylab, xlabp)
  plot(frq, spc, type=type, xlab=xlab, ylab=ylab, ...)
  if (null.line){
    #v<-c(0,0.5)
    v<-range(frq)
    abline(v=v, lty=3, col="blue", lwd=1.5)
    h <- 1
    if (y.dB) h <- 10*log10(c(1,2,10))
    abline(h=h, lty=3, col="red", lwd=1.5)
  }
  if (is.null(main)){
    main <- x@method
  }
  title(main=main)
})

#' @rdname rspec-methods
#' @name lines
#' @exportMethod lines rspec
#' @docType methods
setGeneric("lines", function (x, ...) standardGeneric("lines"))
setMethod("lines", c("rspec"), function(x,...){
  lines(x@freq, x@spec,...)
})
