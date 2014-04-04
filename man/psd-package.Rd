\docType{package}
\name{psd-package}
\alias{psd}
\alias{psd-package}
\alias{spec.psd}
\title{Adaptively estimate power spectral densities of an optimally tapered series.}
\description{
Perform adaptive estimation of the power spectral density
(PSD) using the sine multitapers in which the number of
tapers (and hence the resolution and uncertainty) vary
according to spectral shape. The main function to be used
is \code{\link{pspectrum}}.
}
\details{
In frequency ranges where the spectrum (\eqn{S}) is
relatively flat, more tapers are taken and so a higher
accuracy is attained at the expense of lower frequency
resolution. The program makes a pilot estimate of the
spectrum, then uses Riedel and Sidorenko's estimate of the
MSE (minimum square error) value, which is based on an
estimate of the second derivative of the PSD (\eqn{S''}).
The process is repeated \code{niter} times with a default
of \code{niter=5}. Further iteration may be necessary to
reach convergence, or an acceptably low spectral variance.
Although the term "acceptable" is rather subjective, one
can usually detect an unconverged state by a rather jagged
appearence of the spectrum; this is rather uncommon in our
experience.

\subsection{Adaptive estimation}{ The adaptive process used
is as follows. A quadratic fit to the logarithm of the PSD
within an adaptively determined frequency band is used to
find an estimate of the local second derivative of the
spectrum. This is used in an equation like R-S eq (13) for
the MSE taper number, with the difference that a parabolic
weighting is applied with increasing taper order. Because
the FFTs of the tapered series can be found by resampling
the FFT of the original time series (doubled in length and
padded with zeros) only one FFT is required per series, no
matter how many tapers are used. The spectra associated
with the sine tapers are weighted before averaging with a
parabolically varying weight. The expression for the
optimal number of tapers given by R-S must be modified
since it gives an unbounded result near points where
\eqn{S''} vanishes, which happens at many points in most
spectra. This program restricts the rate of growth of the
number of tapers so that a neighboring covering interval
estimate is never completely contained in the next such
interval. }

\subsection{Resolution and uncertainty}{ The sine
multitaper adaptive process introduces a variable
resolution and error in the frequency domain. See
documentation for \code{\link{spectral_properties}} details
on how these are computed. }
}
\author{
Robert L. Parker and Andrew J. Barbour
<andy.barbour@gmail.com>
}
\references{
Parker, R. L., \emph{PSD}, Program documentation.
\emph{Maintained Software}, N.p. 11 Nov. 2011, Web. 17 Jan.
2013,
<\url{http://igppweb.ucsd.edu/\%7Eparker/Software/\#PSD}>.

Percival, D. B., and A.T. Walden (1993), Spectral analysis
for physical applications, \emph{Cambridge University
Press}

Prieto, G. A., R. L. Parker, D. J. Thomson, F. L. Vernon,
and R. L. Graham (2007), Reducing the bias of multitaper
spectrum estimates, \emph{Geophysical Journal
International}, \strong{171}, 1269--1281, doi:
10.1111/j.1365-246X.2007.03592.x

Riedel, K. S., & Sidorenko, A. (1995), Minimum bias
multiple taper spectral estimation, \emph{Signal
Processing, IEEE Transactions on}, \strong{43}(1),
188--195.

Riedel, K. S. (1996), Adaptive smoothing of the
log-spectrum with multiple tapering, \emph{Signal
Processing, IEEE Transactions on}, \strong{44}(7),
1794--1800.

Walden, A. T., and E. J. McCoy, and D. B. Percival (1995),
The effective bandwidth of a multitaper spectral estimator,
\emph{Biometrika}, \strong{82}(1), 201--214.
}
\seealso{
\code{\link{pspectrum}}
}

