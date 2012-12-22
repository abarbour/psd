\name{ctap_loess}
\alias{constrain_taper_loess_smooth}
\alias{ctap_loess}
\alias{ctap_loess.taper}
\title{ctap_loess}
\usage{
  ctap_loess(tapvec, tapseq = NULL, loess.span = 0.3,
    loess.degree = 1, verbose = TRUE, ...)

  ctap_loess.taper(tapvec, tapseq = NULL, loess.span = 0.3,
    loess.degree = 1, verbose = TRUE, ...)
}
\arguments{
  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{tapseq}{vector; positions or frequencies --
  necessary for smoother methods}

  \item{loess.span}{(\code{'loess.smooth'}) scalar; the
  span used in \code{loess}}

  \item{loess.degree}{(\code{'loess.smooth'}) scalar; the
  polynomial degree}

  \item{verbose}{boolean; should warnings and messages be
  given?}

  \item{...}{optional argments; unused}
}
\value{
  an object with class 'taper', with a constrained taper
  numbers
}
\description{
  Constrain tapers with the Loess smoother.
}
\details{
  This method uses \code{stats::loess} to smooth the taper
  vector; is can be very slow thanks to quadratic scaling.
}
\note{
  The results obtained by \strong{\code{'loess.smooth'}}
  are strongly dependent on the loess tuning parameters;
  hence, some effort should be given to understand their
  effect, and/or re-tuning the default parameters.
}
\examples{
\dontrun{
ntap <- as.taper(1:10)
ctap_loess(ntap)
}
}
\seealso{
  \code{\link{constrain_tapers}}, \code{\link{loess}}
}

