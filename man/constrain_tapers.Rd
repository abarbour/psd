\name{constrain_tapers}
\alias{constrain_tapers}
\alias{constrain_tapers.taper}
\title{constrain_tapers}
\usage{
  constrain_tapers(tapvec, tapseq = NULL,
    constraint.method = c("simple.slope", "markov.chain", "loess.smooth", "friedman.smooth", "none"),
    min.tapers = 1, verbose = TRUE, ...)

  constrain_tapers.taper(tapvec, tapseq = NULL,
    constraint.method = c("simple.slope", "markov.chain", "loess.smooth", "friedman.smooth", "none"),
    min.tapers = 1, verbose = TRUE, ...)
}
\arguments{
  \item{tapvec}{vector; the number of tapers at each
  frequency}

  \item{tapseq}{vector; positions or frequencies --
  necessary for smoother methods}

  \item{constraint.method}{method to use for constraints on
  tapers numbers}

  \item{min.tapers}{integer; the minimum number of tapers}

  \item{verbose}{logical; should warnings and messages be
  given?}

  \item{...}{optional argments passed to the
  constraint.method}
}
\description{
  Apply constraints on the number of tapers
}
\details{
  Refines the number of tapers; the method by which it does
  this may be chosen by the user.
}
\author{
  Andrew Barbour <andy.barbour@gmail.com> ported original
  by R.L.Parker to R and C.
}
\seealso{
  \code{\link{ctap_simple}}, \code{\link{ctap_markov}},
  \code{\link{riedsid}}
}

