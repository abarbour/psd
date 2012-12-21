\docType{methods}
\name{constrain_tapers}
\alias{.ctap_friedman.default}
\alias{.ctap_loess.default}
\alias{.ctap_markov.default}
\alias{.ctap_simple.default}
\alias{constrain_tapers}
\title{constrain_tapers}
\usage{
  constrain_tapers(tapvec, tapseq = NULL,
    constraint.method = c("simple.slope", "markov.chain", "loess.smooth", "friedman.smooth", "none"),
    min.tapers = 1, verbose = TRUE, ...)

  \method{'simple.slope'}{constrain_tapers} (tapvec,
    tapseq = NULL, maxslope = 1L)

  \method{'markov.chain'}{constrain_tapers} (tapvec,
    tapseq = NULL, chain.width = 3 * length(tapvec))

  \method{'loess.smooth'}{constrain_tapers} (tapvec,
    tapseq, loess.span = 0.3, loess.degree = 1)

  \method{'friedman.smooth'}{constrain_tapers} (tapvec,
    tapseq, smoo.span = 0.3, smoo.bass = 2)
}
\arguments{
  \item{tapvec}{vector; the number of tapers at each
  frequency}

  \item{tapseq}{vector; positions or frequencies --
  necessary for smoother methods}

  \item{constraint.method}{method to use for constraints on
  tapers numbers}

  \item{min.tapers}{integer; the minimum number of tapers}

  \item{verbose}{boolean; should warnings and messages be
  given?}

  \item{...}{optional argments passed to the
  constraint.method}

  \item{maxslope}{(\code{'simple.slope'}) integer;
  constrain based on this maximum first difference}

  \item{chain.width}{(\code{'markov.chain'}) scalar; the
  width the MC should consider for the change probability}

  \item{loess.span}{(\code{'loess.smooth'}) scalar; the
  span used in \code{loess}}

  \item{loess.degree}{(\code{'loess.smooth'}) scalar; the
  polynomial degree}

  \item{smoo.span}{(\code{'friedman.smooth'}) scalar;
  fraction of the observations in the span of the running
  lines smoother}

  \item{smoo.bass}{(\code{'friedman.smooth'}) scalar;
  controls the smoothness of the fitted curve}
}
\value{
  \code{NULL}
}
\description{
  Constrain the number of tapers applied

  \strong{\code{'simple.slope'}} constrains the
  first-difference between neighboring spectra.

  \strong{\code{'markov.chain'}} uses a Morkov Chain based
  on the theory of quantum-well probability in gamma-ray
  spectroscopy.

  \strong{\code{'loess.smooth'}} uses \code{stats::loess}
  for the smoothing.

  \strong{\code{'friedman.smooth'}} uses builtin Fiedman
  super-smoother.
}
\details{
  Refines the number of tapers; the method which it does
  this may be chosen by the user.

  \strong{\code{'simple.slope'}}: Constrain tapers with
  slopes from first differencing

  \strong{\code{'markov.chain'}}: Constrain tapers with a
  Markov Chain, based on quantum-well probability

  \strong{\code{'loess.smooth'}}: Constrain tapers with the
  Loess smoother (slow)

  \strong{\code{'friedman.smooth'}}: Constrain tapers with
  the Friedman 'super-smoother'
}
\note{
  The results obtained by \strong{\code{'friedman.smooth'}}
  are generally poor; hence, the method may be removed in
  future releases.
}
\seealso{
  \code{\link{riedsid}}

  \code{\link{loess}}

  \code{\link{supsmu}}
}

