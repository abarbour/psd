% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/func_norm.R
\name{psd-normalization}
\alias{psd-normalization}
\alias{normalization}
\alias{normalize}
\alias{normalize.list}
\alias{normalize.spec}
\alias{normalize.amt}
\title{Normalization of power spectral density estimates.}
\usage{
normalize(Spec, ...)

\method{normalize}{list}(Spec, ...)

\method{normalize}{spec}(
  Spec,
  Fsamp = 1,
  src = c("spectrum", "double.sided", "psd", "single.sided"),
  verbose = TRUE,
  ...
)

\method{normalize}{amt}(Spec, ...)
}
\arguments{
\item{Spec}{spectrum to normalize}

\item{...}{(unused) additional parameters}

\item{Fsamp}{sampling frequency}

\item{src}{character string; the source of the spectrum estimator}

\item{verbose}{logical; should messages be given?}
}
\value{
An object with its spectral values normalized accordingly.
}
\description{
Normalize power spectral densities from
various estimators into single-sided spectra.
}
\details{
Normalizations commonly encountered for power spectra 
depend on it's assumed sidedness: whether the
spectrum is either single- or double-sided.
The normalizations performed here enforce single-sidedness, and correct
as necessary.

Frequencies are assumed to be based on the Nyquist frequency (half the 
sampling rate).  For example: If a series \eqn{X} has sampling frequency \eqn{F_S},
then the PSD frequencies will span \eqn{[0,F_S/2]}.

For amplitudes, improper normalization can can introduce errant factors
of either 1/2 or \eqn{F_S} into the estimates, depending on the assumed sidedness.  
These factors can be accounted for with the \code{src}
argument, which defaults to normalizing a double-sided spectrum.
}
\section{Spectrum sidedness and the \code{src} argument}{

\subsection{\code{"double.sided"} or \code{"spectrum"}}{

These spectra assume frequency range of \eqn{[-F_S/2,F_S/2]}, and so are normalized
by scaling by a factor of two upwards.
Some estimators producing double-sided spectra: 
\itemize{
\item{\code{stats::spectrum}}{}
\item{\code{RSEIS::mtapspec}}{}
}
}

\subsection{\code{"single.sided"} or \code{"psd"}}{
As mentioned before, 
these spectra assume frequency range of \eqn{[0,F_S/2]} and
are scaled only by the inverse of the sampling rate.
Some estimators producing single-sided spectra: 
\itemize{
\item{\code{\link{psdcore}}}{}
}
}
}

\examples{
\dontrun{#REX
library(psd)

##
## Normalization
##

# timeseries with sampling frequency **not** equal to 1:
set.seed(1234)
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
mean(pgram[['spec']]) * max(pgram[['freq']])
mean(PSD[['spec']]) * max(PSD[['freq']])

# normalize 
pgram <- normalize(pgram, src="spectrum")
PSD <- normalize(pgram, src="psd")
# replot them
plot(pgram, log="dB", ylim=c(-40,10))
plot(PSD, add=TRUE, col="red", log="dB")

# Again, integrated spectrum should be ~= 1:
mean(pgram[['spec']]) * max(pgram[['freq']])
mean(PSD[['spec']]) * max(PSD[['freq']])


}#REX
}
\seealso{
\code{\link{psdcore}}, \code{\link{spectral_properties}}
}
\author{
A.J. Barbour
}
