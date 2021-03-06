% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{riedsid_rcpp}
\alias{riedsid_rcpp}
\title{replaces time consuming portion of riedsid2}
\usage{
riedsid_rcpp(PSD, ntaper, riedsid_column = 0L)
}
\arguments{
\item{PSD}{vector or class \code{'amt'} or \code{'spec'}; the spectral values used to optimize taper numbers}

\item{ntaper}{scalar or vector; number of tapers to apply optimization}

\item{riedsid_column}{scalar integer; which column to use in multivariate optimization. If the value is 0 the maximum number of tapers for all columns is chosen. If the value is < 0 the minimum number of tapers for all columns is chosen. If the value is 1, 2, 3, etc. the number of tapers is based on the column selected.}
}
\value{
kopt vector
}
\description{
replaces time consuming portion of riedsid2
}
