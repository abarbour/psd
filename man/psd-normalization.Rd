\name{psd-normalization}
\alias{normalization}
\alias{normalize}
\alias{normalize.default}
\alias{normalize.list}
\alias{normalize.spec}
\alias{psd-normalization}
\title{Normalization of power spectral density estimates.}
\usage{
normalize(Spec, Fsamp = 1, src = NULL, verbose = TRUE, ...)

\method{normalize}{default}(Spec, Fsamp = 1, src = NULL, verbose = TRUE,
  ...)

\method{normalize}{list}(Spec, Fsamp = 1, src = NULL, verbose = TRUE, ...)

\method{normalize}{spec}(Spec, Fsamp = 1, src = NULL, verbose = TRUE, ...)
}
\arguments{
  \item{Spec}{spectrum to normalize}

  \item{Fsamp}{sampling frequency}

  \item{src}{character string; the source of the spectrum
  estimator}

  \item{verbose}{logical; should messages be given?}

  \item{...}{(unused) additional parameters}
}
\value{
An object with its spectral values normalized accordingly.
}
\description{
Normalize power spectral densities from various estimators
into single-sided spectra.
}
\details{
Normalizations commonly encountered for power spectra
depend on it's assumed sidedness: whether the spectrum is
either single- or double-sided. The normalizations
performed here enforce single-sidedness, and correct as
necessary.

Frequencies are assumed to be based on the Nyquist
frequency (half the sampling rate).  For example: If a
series \eqn{X} has sampling frequency \eqn{F_S}, then the
PSD frequencies will span \eqn{[0,F_S/2]}.

For amplitudes, improper normalization can can introduce
errant factors of either 1/2 or \eqn{F_S} into the
estimates, depending on the assumed sidedness. These
factors can be accounted for with the \code{src} argument,
which defaults to normalizing a double-sided spectrum.
}
\section{Spectrum sidedness and the \code{src} argument}{
  \subsection{\code{"double.sided"} or \code{"spectrum"}}{

  These spectra assume frequency range of
  \eqn{[-F_S/2,F_S/2]}, and so are normalized by scaling by
  a factor of two upwards. Some estimators producing
  double-sided spectra: \itemize{
  \item{\code{stats::spectrum}}{}
  \item{\code{RSEIS::mtapspec}}{} } }

  \subsection{\code{"single.sided"} or \code{"psd"}}{ As
  mentioned before, these spectra assume frequency range of
  \eqn{[0,F_S/2]} and are scaled only by the inverse of the
  sampling rate. Some estimators producing single-sided
  spectra: \itemize{ \item{\code{\link{psdcore}}}{} } }
}
\examples{
#RDEX#\dontrun{
require(psd)
##
## Normalization
##
set.seed(1234)
# timeseries with sampling frequency **not** equal to 1:
X <- ts(rnorm(1e3), frequency=20)
# spec.pgram: double sided
pgram <- spectrum(X)
# psdcore: single sided
PSD <- psdcore(X)
# note the normalization differences:
plot(pgram, log="dB", ylim=c(-40,10))
plot(PSD, add=TRUE, col="red", log="dB")
# A crude representation of integrated spectrum: 
#   should equal variance of white noise series (~= 1)
mean(pgram$spec)*max(pgram$freq)
mean(PSD$spec)*max(PSD$freq)
#
# normalize objects with class 'spec'
pgram <- normalize(pgram, src="spectrum")
PSD <- normalize(pgram, src="psd")
# replot them
plot(pgram, log="dB", ylim=c(-40,10))
plot(PSD, add=TRUE, col="red", log="dB")
#
# Again, integrated spectrum should be ~= 1:
mean(pgram$spec)*max(pgram$freq)
mean(PSD$spec)*max(PSD$freq)
#
#RDEX#}
}
\author{
A.J. Barbour <andy.barbour@gmail.com>
}
\seealso{
\code{\link{psdcore}}, \code{\link{spectral_properties}}
}
\keyword{normalization}
\keyword{prewhiten}
\keyword{spectrum-estimation}

