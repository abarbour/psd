\name{constrain_tapers}
\alias{constrain_tapers}
\title{constrain_tapers}
\usage{
  constrain_tapers(tapvec, tapseq = NULL,
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

  \item{verbose}{boolean; should warnings and messages be
  given?}

  \item{...}{optional argments passed to the
  constraint.method}
}
\description{
  Constrain the number of tapers applied
}
\details{
  Refines the number of tapers; the method which it does
  this may be chosen by the user.
}
\seealso{
  \code{\link{riedsid}}, \code{\link{ctap_simple}}
}

