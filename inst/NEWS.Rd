% borrowed some formatting from 
% http://svn.r-project.org/R/trunk/doc/NEWS.Rd
% -*- coding: utf-8 -*-
%%\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}
%%\newcommand{\Rlogo}{\if{html}{\figure{../../html/logo.jpg}}\if{latex}{\figure{logo.jpg}{options: width=0.5in}}}
\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
\newcommand{\psd}{\CRANpkg{psd}}
\newcommand{\citation}{\href{http://dx.doi.org/10.1016/j.cageo.2013.09.015}{10.1016/j.cageo.2013.09.015}}

\name{psd}
\title{News for package 'psd'}
\encoding{UTF-8}

% template vim:5Y
%
%\section{CHANGES IN \psd{} VERSION X.Y-Z (XXXX-YY-ZZ)}{
%  \itemize{ \item \code{} }
%}
%

\section{CHANGES IN \psd{} VERSION 0.4-1 (2014-04-15)}{
  \itemize{
  	\item We have published a paper in Computers \& Geosciences regarding
	\psd{}, and the citation is up to date: \code{citation('psd')}
	If you cannot access the publication, contact A.~J.~Barbour
	for a reprint.
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
