% borrowed some formatting from 
% http://svn.r-project.org/R/trunk/doc/NEWS.Rd
% -*- coding: utf-8 -*-
%%\newcommand{\PR}{\Sexpr[results=rd]{tools:::Rd_expr_PR(#1)}}
%%\newcommand{\Rlogo}{\if{html}{\figure{../../html/logo.jpg}}\if{latex}{\figure{logo.jpg}{options: width=0.5in}}}
\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

\name{psd}
\title{News for package 'psd'}
\encoding{UTF-8}
\section{CHANGES IN psd VERSION 0.3-1 (2013-03-12)}{
  \itemize{
    \item This marks the first release to CRAN: \CRANpkg{psd}.
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
\section{CHANGES IN psd VERSION 0.3-0 (2013-03-11)}{
  \itemize{
    \item Added \code{Tohoku} data.
  }
}
