% borrowed some formatting from 
% http://svn.r-project.org/R/trunk/doc/NEWS.Rd
% -*- coding: utf-8 -*-
%%\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}
%%\newcommand{\Rlogo}{\if{html}{\figure{../../html/logo.jpg}}\if{latex}{\figure{logo.jpg}{options: width=0.5in}}}
\newcommand{\GITpkg}{\href{http://www.github.com/abarbour/#1}{github}}
\newcommand{\CRANpkg}{\href{http://cran.r-project.org/package=#1}{CRAN}}
\newcommand{\psd}{\pkg{psd}}
\newcommand{\gitpsd}{\GITpkg{psd} page}
\newcommand{\citation}{\href{http://dx.doi.org/10.1016/j.cageo.2013.09.015}{10.1016/j.cageo.2013.09.015}}

\name{psd}
\title{News for package \pkg{psd}: \CRANpkg{psd} -- \GITpkg{psd} -- \href{https://travis-ci.org/abarbour/psd}{Travis}}

\encoding{UTF-8}

% template vim:6Y from below here
%
%\section{CHANGES IN \psd{} VERSION X.Y-Z (XXXX-YY-ZZ)}{
%  \itemize{ \item \code{} }
%}
%

%
% 1.0 -- the era of clean code and Rcpp goodness begins
%
\section{CHANGES IN \psd{} VERSION 1.0-0 (2015-03-17)}{
  \itemize{
    
    \item This is a major version bump that addresses performance issues by
    converting known bottlenecks to c++ though the 
    \href{http://cran.r-project.org/package=Rcpp}{\pkg{Rcpp}} package.
    There has also been quite a lot of cleanup and rearranging of
    documentation, and making sure methods are consistent.
    Some attempts to be backwards compatible have been made, but
    don't be surprised if results from previous versions are different --
    sorry!

    \item \code{psdcore} now uses \code{resample_fft_rcpp} which
    is a lightning-fast (by comparison with pure-R) implementation
    of the fft-resampling/reweighting method we use; this new function
    leverages the power of 
    \href{http://cran.r-project.org/package=RcppArmadillo}{\pkg{RcppArmadillo}}.
    
    \item \code{ctap_simple} has been superceded by \code{ctap_cimple_rcpp}; in
    the process of rewriting in c++, a minor bug in the c implementation was
    found (and fixed).

    \item \href{https://travis-ci.org/abarbour/psd}{Travis-CI} is now used to track builds; see 
    also the \gitpsd{}.
    
    \item Issues/pull-requests can be entered on the \href{https://github.com/abarbour/psd/issues}{Issues} page.
    
    \item Taper constraints: \code{ctap_simple_rcpp}/\code{ctap_simple} and \code{ctap_loess} gain
    their own man-page; \code{ctap_markov} and \code{ctap_friedman} were
    previously made defunct, but are now permanently removed from the code-base.

  }
}


\section{CHANGES IN \psd{} VERSION 0.4-3 (2014-10-16)}{
  \itemize{ 
    \item Thanks to David Myer for catching a number of mistakes
    in the code and documentation,
    including a bug in \code{riedsid} from an improper choice
    of logarithm base.
    
    \item Removed the \code{Local.loss} argument from \code{riedsid}
  }
}

\section{CHANGES IN \psd{} VERSION 0.4-2 (2014-08-05)}{
  \itemize{ 
    \item Removed unused variables in \code{ctap_simple.c} 
    \item Removed deprecated \code{@S3method} statements.
  }
}

\section{CHANGES IN \psd{} VERSION 0.4-1 (2014-04-15)}{
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
\section{CHANGES IN \psd{} VERSION 0.4-0 (2013-09-19)}{
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
\section{CHANGES IN \psd{} VERSION 0.3-2 (2013-03-27)}{
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
\section{CHANGES IN \psd{} VERSION 0.3-1 (2013-03-12)}{
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
\section{CHANGES IN \psd{} VERSION 0.3-0 (2013-03-11)}{
  \itemize{
    \item Added \code{Tohoku} data.
  }
}
