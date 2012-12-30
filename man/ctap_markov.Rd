\name{ctap_markov}
\alias{constrain_taper_markov_chain}
\alias{ctap_markov}
\alias{ctap_markov.taper}
\title{ctap_markov}
\usage{
  ctap_markov(tapvec, tapseq = NA,
    chain.width = round(5 * length(tapvec)), ...)

  ctap_markov.taper(tapvec, tapseq = NA,
    chain.width = round(5 * length(tapvec)), ...)
}
\arguments{
  \item{tapvec}{'taper' object; the number of tapers at
  each frequency}

  \item{tapseq}{(unused) vector; positions or frequencies
  -- necessary for smoother methods}

  \item{chain.width}{(\code{'markov.chain'}) scalar; the
  width the MC should consider for the change probability}

  \item{...}{(unused) optional argments}
}
\value{
  an object with class 'taper', with a constrained taper
  numbers
}
\description{
  Constrain tapers with a Markov Chain, based on
  quantum-well probability.
}
\details{
  This method uses a Markov Chain, based on the theory of
  quantum-well probability, commonly used in gamma-ray
  spectroscopy.

  The main function behind this method is
  \link[Peaks]{SpectrumSmoothMarkov} which calculates the
  probability that the number of tapers would have changed
  (from it's previous value); it is very fast.  Details of
  the theory behind this algorithm may be found in Morhac
  (2008) and Silagadze (1996).
}
\note{
  The results obtained by \strong{\code{'markov.chain'}}
  can produce unstable results if successively executed on
  taper vectors with even modest sized serially-correlated
  peaks; the method should be used with caution, unless the
  intention is to enhance peaks.  The \code{'chain.width'}
  parameter constrols the broadness of the a priori
  distribution.  As a rule of thumb: the smaller the
  parameter is, the shorter the tails become.
}
\examples{
\dontrun{
ntap <- as.taper(1:10)
ctap_markov(ntap)
}
}
\author{
  Andrew Barbour <andy.barbour@gmail.com> adapted algorithm
  from Morhac (2008).
}
\references{
  Morhac, M. (2008), Peaks: Peaks, \emph{R package},
  \strong{version 0.2}

  Silagadze, Z.K. (1996), A new algorithm for automatic
  photopeak searches, \emph{Nucl. Instrum. Meth. A},
  \strong{376} 451.

  \url{http://arxiv.org/abs/hep-ex/9506013}
}
\seealso{
  \code{\link{constrain_tapers}},
  \code{\link[Peaks]{SpectrumSmoothMarkov}}
}

