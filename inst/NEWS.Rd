% borrowed some formatting from 
% http://svn.r-project.org/R/trunk/doc/NEWS.Rd
% -*- coding: utf-8 -*-
%%\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}
%%\newcommand{\Rlogo}{\if{html}{\figure{../../html/logo.jpg}}\if{latex}{\figure{logo.jpg}{options: width=0.5in}}}
\newcommand{\GITpkg}{\href{https://www.github.com/abarbour/#1}{github}}
\newcommand{\CRANpkg}{\href{https://cran.r-project.org/package=#1}{CRAN}}
\newcommand{\CODEcov}{\href{https://codecov.io/github/abarbour/#1?branch=master}{codecov.io}}
\newcommand{\Appvey}{\href{https://ci.appveyor.com/project/abarbour/#1}{AppVeyor}}
\newcommand{\Trav}{\href{https://travis-ci.org/abarbour/#1}{Travis-CI}}
\newcommand{\psd}{\pkg{psd}}
\newcommand{\gitpsd}{\GITpkg{psd} page}
\newcommand{\citation}{\href{https://doi.org/10.1016/j.cageo.2013.09.015}{10.1016/j.cageo.2013.09.015}}

\name{psd}
\title{News for package \pkg{psd}: \CRANpkg{psd} -- \GITpkg{psd} -- \Trav{psd} -- \Appvey{psd} -- \CODEcov{psd}}

\encoding{UTF-8}

% template vim:6Y from below here
%
%\section{CHANGES IN \psd{} VERSION X.Y-Z (XXXX-YY-ZZ)}{
%  \itemize{ \item \code{} }
%}
%

\section{CHANGES IN \psd{} VERSION 2.1.0 (2020-06-28)}{
  \itemize{
    \item Fixed an issue related to how psd initiates its internal options.
    \item New flag in \code{psdcore} indivating if calculations are multivariate or not; this
    will be used in the future to make plotting functions more robust.
    \item With the addition of a new function \code{na_locf},
    the dependency on \pkg{zoo} has been removed.
  }
}

\section{CHANGES IN \psd{} VERSION 2.0.0 (2020-06-21)}{
  \itemize{
    \item Major bump which adds multivariate calculations (i.e., the cross spectrum)
    using the same optimization feature.
  }
}

\section{CHANGES IN \psd{} VERSION 1.2.1 (2019-03-21)}{
  \itemize{
    \item Minor bump to fix incorrect orcid no. in description file
  }
}

\section{CHANGES IN \psd{} VERSION 1.2.0 (2019-03-20)}{
  \itemize{
    \item
    Significant speed improvements in spectral
    derivative computation implemented in \code{riedsid2} 
    thanks to \href{https://github.com/jkennel}{@jkennel}; 
    this will be most apparent for very long
    timeseries.
    \item \code{riedsid} is deprecated
    \item \code{niter} in \code{pspectrum} default set down from 5 to 3
    \item Modernized revdep checks
  }
}

\section{CHANGES IN \psd{} VERSION 1.1.0 (2017-07-13)}{
  \itemize{
    \item \code{ctap_simple_rcpp} replaces \code{ctap_simple}; associated .c code deprecated 
  	\item New \code{.travis.yml} and \CODEcov{psd} unit-test coverage 
  	\item Windows CI with \Appvey{psd}
  	\item Fixing new NAMESPACE requirements with \code{importFrom} statements
  	\item Reverse-dependency checking in \code{revdep/}
  }
}

%
% 1.0 -- the era of clean code and Rcpp goodness begins
%
\section{VERSION 1.0-1 (2015-03-19)}{
  \itemize{
	  \item Addresses three minor c++ issues associated with type overloads (not seen until CRAN-build-checks):
	  \enumerate{
		  \item \code{resample_fft.cpp:92:22}
		  \item \code{resample_fft.cpp:104:31}
	  	\item \code{resample_fft.cpp:264:34}
	  }
	  \item Prevents test-failure when fftw is not available.  Thanks again, fftw!
  }
}

\section{VERSION 1.0-0 (2015-03-17)}{
  \itemize{
    
    \item This is a major version bump that addresses performance issues by
    converting known bottlenecks to c++ though the 
    \href{https://cran.r-project.org/package=Rcpp}{\pkg{Rcpp}} package.
    There has also been quite a lot of cleanup and rearranging of
    documentation, and making sure methods are consistent.
    Some attempts to be backwards compatible have been made, but
    don't be surprised if results from previous versions are different --
    sorry!

    \item \code{psdcore} now uses \code{resample_fft_rcpp} which
    is a lightning-fast (by comparison with pure-R) implementation
    of the fft-resampling/reweighting method we use; this new function
    leverages the power of 
    \href{https://cran.r-project.org/package=RcppArmadillo}{\pkg{RcppArmadillo}}.
    
    \item \code{ctap_simple} has been superceded by \code{ctap_simple_rcpp}; in
    the process of rewriting in c++, a minor bug in the c implementation was
    found (and fixed).

    \item \Trav{psd} is now used to track builds; see also the \gitpsd{}.
    
    \item Issues/pull-requests can be entered on the \href{https://github.com/abarbour/psd/issues}{Issues} page.
    
    \item Taper constraints: \code{ctap_simple_rcpp}/\code{ctap_simple} and \code{ctap_loess} gain
    their own man-page; \code{ctap_markov} and \code{ctap_friedman} were
    previously made defunct, but are now permanently removed from the code-base.

  }
}


\section{VERSION 0.4-3 (2014-10-16)}{
  \itemize{ 
    \item Thanks to David Myer for catching a number of mistakes
    in the code and documentation,
    including a bug in \code{riedsid} from an improper choice
    of logarithm base.
    
    \item Removed the \code{Local.loss} argument from \code{riedsid}
  }
}

\section{VERSION 0.4-2 (2014-08-05)}{
  \itemize{ 
    \item Removed unused variables in \code{ctap_simple.c} 
    \item Removed deprecated \code{@S3method} statements.
  }
}

\section{VERSION 0.4-1 (2014-04-15)}{
  \itemize{
  	\item We have published a paper in Computers & Geosciences regarding
	\psd{}, and the citation is up to date: \code{citation('psd')}
	If you cannot access the publication (\citation{}), 
	contact Andy for a reprint.
  	\item Fixed an example which was causing CRAN-buildcheck errors.
	\item Fixed a small error in the normalization vignette -- thanks 
	      to Richard Gaal for spotting this.
  }
}
%
%
\section{VERSION 0.4-0 (2013-09-19)}{
  \itemize{
    \item The taper-constraint methods \code{ctap_markov} and \code{ctap_friedman} are now defunct.
    \item Removed 'units.sty' dependency in vignettes (was causing NOTES and unhappiness).
    \item \code{VignetteEngine} added, and \code{:::}'s removed, in order
    to ensure compatibility with R-devel.
    \item Added citation information for upcoming paper in Computers and Geoscience.
  }
}
%
%
\section{VERSION 0.3-2 (2013-03-27)}{
  \itemize{
    \item Revised discussions regarding AR response spectrum.
    \item Fixed \code{psdcore}: it was applying parabolic weights incorrectly.
    \item Changed argument name \code{ntap_pilot} in \code{pspectrum} to \code{ntap.init}.
    \item Changed \code{ceiling} in \code{minspan} to \code{round}.
    \item Fixed \code{km} field in \code{magnet} dataset.
    \item Fixed referencing s.t. \code{?plot} isn't ambiguous.
    \item The default number of iterations and pilot tapers
          was changed to 3 and 7 respectively.
  }
}
%
\section{VERSION 0.3-1 (2013-03-12)}{
  \itemize{
    \item This marks the first release to CRAN: \psd{}.
          We have included three vignettes, one of which is an overview of
          the general functionality of the program. These may be listed with
	  \code{vignette(pack="psd")}, but we recommend consulting
	  \code{vignette("psd_overview", package="psd")} to begin.
    \item Altered the environment manipulation features to prevent .GlobalEnv
    modification.
    \item Fixed CRAN-check failure with \code{usepackage{color}} in vignettes.
    \item Removed empty \code{.Last.lib} function to prevent CRAN-check note.
  }
}
\section{VERSION 0.3-0 (2013-03-11)}{
  \itemize{
    \item Added \code{Tohoku} data.
  }
}
\section{VERSION 0.1-0}{
	\itemize{
		\item First working version
	}
}
