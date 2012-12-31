\name{ctap_friedman}
\alias{constrain_taper_friedman_smooth}
\alias{ctap_friedman}
\alias{ctap_friedman.taper}
\title{ctap_friedman}
\usage{
  ctap_friedman(tapvec, tapseq = NULL, smoo.span = 0.3,
    smoo.bass = 2, verbose = TRUE, ...)

  ctap_friedman.taper(tapvec, tapseq = NULL,
    smoo.span = 0.3, smoo.bass = 2, verbose = TRUE, ...)
}
\arguments{
  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{tapseq}{vector; positions or frequencies --
  necessary for smoother methods}

  \item{smoo.span}{scalar; fraction of the observations in
  the span of the running lines smoother}

  \item{smoo.bass}{scalar; controls the smoothness of the
  fitted curve}

  \item{verbose}{logical; should warnings and messages be
  given?}

  \item{...}{optional argments; unused}
}
\value{
  an object with class 'taper', with a constrained taper
  numbers
}
\description{
  Constrain tapers with the Friedman 'super-smoother'.
}
\details{
  This method uses the base Fiedman super-smoother.
}
\note{
  The results obtained by \strong{\code{'friedman.smooth'}}
  are generally poor; hence, the method may be removed in
  future releases.
}
\examples{
\dontrun{
ntap <- as.taper(1:10)
ctap_friedman(ntap)
}
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{constrain_tapers}}, \code{\link{supsmu}}
}

