#' Perform adaptive estimation 
#' of the power spectral density (PSD) using the sine multitapers in which
#' the number of tapers (and hence the resolution and uncertainty) 
#' vary according to spectral shape. 
#' The main function to be used 
#' is \code{\link{pspectrum}}.
#'
#' @details
#' In frequency ranges where the spectrum  (\eqn{S})
#' is relatively flat, more tapers are taken and so a higher accuracy is 
#' attained at the expense of lower frequency resolution. 
#' The program makes a pilot estimate of the spectrum, then uses
#' Riedel and Sidorenko's estimate of the MSE (minimum square error) value, 
#' which is based on an estimate of the second derivative of the PSD (\eqn{S''}). 
#' The process is repeated \code{niter} times with a default of \code{niter=5}. 
#' Further iteration may be necessary to reach convergence, or an acceptably low
#' spectral variance. Although the term "acceptable" is rather subjective, one can 
#' usually detect an unconverged state by a rather jagged appearence of the spectrum;
#' this is rather uncommon in our experience.
#'
#' \subsection{Adaptive estimation}{
#' The adaptive process used is as follows. A quadratic fit to the logarithm of the
#' PSD within an 
#' adaptively determined frequency band is used to find an estimate of the local second 
#' derivative of the spectrum. This is used in an equation like R-S eq (13) for 
#' the MSE taper number, with the difference that a parabolic weighting is applied with 
#' increasing taper order. Because the FFTs of the tapered series can be found by 
#' resampling the FFT of the original time series (doubled in length and padded with zeros) 
#' only one FFT is required per series, no matter how many tapers are used. 
#' The spectra associated with the sine tapers are weighted before averaging with a 
#' parabolically varying weight. The expression for the optimal number of tapers 
#' given by R-S must be modified since it gives an unbounded result near points 
#' where \eqn{S''} vanishes, which happens at many points in most spectra. 
#' This program restricts the rate of growth of the number of tapers so that a 
#' neighboring covering interval estimate is never completely contained in the next 
#' such interval.
#' }
#'
#' \subsection{Resolution and uncertainty}{
#' The sine multitaper adaptive process 
#' introduces a variable resolution and error in the frequency domain. 
#' See documentation for \code{\link{spectral_properties}} details on
#' how these are computed.
#' }
#'
#' @docType package
#' @name rlpSpec-package
#' @aliases rlpSpec spec.rlp
#' @title Adaptively estimate power spectral densities of an optimally tapered series.
#' 
#' @author Robert L. Parker and Andrew J. Barbour <andy.barbour@@gmail.com> 
#' 
#' @import Peaks RColorBrewer signal zoo
# not needed since specified in DESCRIPTION:DEPENDS stats utils graphics grDevices
#' @useDynLib rlpSpec
#'
#' @references Parker, R. L., \emph{PSD}, Program documentation. \emph{Maintained Software}, N.p. 11 Nov. 2011,
#' Web. 17 Jan. 2013, <\url{http://igppweb.ucsd.edu/\%7Eparker/Software/\#PSD}>.
#'
#' @references Percival, D. B., and A.T. Walden (1993),
#' Spectral analysis for physical applications,
#' \emph{Cambridge University Press}
#'
#' @references Prieto, G. A., R. L. Parker, D. J. Thomson, F. L. Vernon, and R. L. Graham  (2007), 
#' Reducing the bias of multitaper spectrum estimates,
#' \emph{Geophysical Journal International}, \strong{171}, 1269--1281,
#' doi: 10.1111/j.1365-246X.2007.03592.x
#' 
#' @references Riedel, K. S., & Sidorenko, A. (1995), 
#' Minimum bias multiple taper spectral estimation,
#' \emph{Signal Processing, IEEE Transactions on}, \strong{43}(1), 188--195.
#
#' @references Riedel, K. S. (1996),
#' Adaptive smoothing of the log-spectrum with multiple tapering,
#' \emph{Signal Processing, IEEE Transactions on}, \strong{44}(7), 1794--1800.
#'
#' @references Walden, A. T., and  E. J. McCoy, and D. B. Percival (1995),
#' The effective bandwidth of a multitaper spectral estimator,
#' \emph{Biometrika}, \strong{82}(1), 201--214.
# \url{http://biomet.oxfordjournals.org/content/82/1/201}
#'
# @seealso \code{\link{pspectrum}}, \code{\link{psdcore}}
#'  
NULL

##
## Datasets
##

#' A single line of MAGSAT horizontal field intensity.
#' 
#' The Magnetic Satellite (MAGSAT) mission 
#' provided a wealth of airborne-magnetometer data
#' spanning the globe (Langel el al., 1982).  
#' This dataset represents a single track of horizontal field
#' intensities (a very small subset of the full collection!).
#'
#' \subsection{Raw and Clean Sets}{
#' There are non-real data points in raw MAGSAT series; these are instrumental
#' artefacts, and can severely affect
#' power spectral density (PSD) estimates.  
#' A clean series has been included
#' so that a comparison of PSDs may be made.
#'
#' Some command like \code{subset(magsat, abs(mdiff) > 0)}
#' can be used to identify the rows where edits have been made.
#' }
#' 
#' @name magsat
#' @docType data
#' @format A dataframe with 2048 observations on the following 4 variables.
#'
#' \describe{
#' \item{\code{km}}{Relative along-track distance, in kilometers. The first observation is at zero kilometers.}
#' \item{\code{raw}}{Raw intensities, in nanoTesla.}
#' \item{\code{clean}}{Edited raw intensites, in nanoTesla}
#' \item{\code{mdiff}}{The difference between \code{clean} and \code{raw} intensities, in nanoTesla.}
#' }
#'
#' @references Langel, R., G. Ousley, J. Berbert, J. Murphy, and M. Settle, (1982),
#' The MAGSAT mission, \emph{Geophysical Research Letters}, \strong{9} (4), 243-245.
#' @source Data subsetter: \url{http://ftpbrowser.gsfc.nasa.gov/magsat.html} 
# need login from outside browser?
# @source Entire MAGSAT archive: \url{ftp://spdf.gsfc.nasa.gov/pub/data/magsat/}
#' @keywords datasets magsat
#' @examples
#' data(magsat)
#' summary(magsat)
NULL
