\name{riedsid}
\alias{riedsid}
\alias{riedsid.default}
\title{riedsid}
\usage{
  riedsid(spec, ntaper, tapseq = NULL, c.method = NULL,
    ...)

  riedsid.default(spec, ntaper, tapseq = NULL,
    c.method = NULL, ...)
}
\arguments{
  \item{spec}{vector; the spectral values used to optimize
  taper numbers}

  \item{ntaper}{scalar or vector; number of tapers to apply
  optimization}

  \item{tapseq}{vector; representing positions or
  frequencies (same length as spec)}

  \item{c.method}{string; constraint method to use}

  \item{...}{optional argments passed to
  \code{constrain_tapers}}
}
\value{
  \code{NULL}
}
\description{
  Riedel & Sidorenko taper optimization
}
\details{
  Estimates optimal number of tapers at each frequency of
  given spec, based on Riedel-Sidorenko MSE recipe, and
  other tweaks by RLP.
}
\seealso{
  \code{\link{psdcore}}, \code{\link{constrain_tapers}}
}

