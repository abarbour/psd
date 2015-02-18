# Makefile to use knitr for package vignettes
# (modified from knitr package source)
#
# DO NOT EXCLUDE FROM BUILD!

TEXS = fftw.tex normalization.tex psd_overview.tex
PDFS = fftw.pdf normalization.pdf psd_overview.pdf

all: $(PDFS)

clean:
	rm -rf *-concordance.tex *.bbl *.blg *.aux *.out *.log *.spl *.md figure/ *.idx *.ilg *.ind *.toc
	rm -f $(TEXS)

%.pdf: %.Rnw
	$(R_HOME)/bin/Rscript -e "if (getRversion() < '3.0.0') knitr::knit2pdf('$*.Rnw') else tools::texi2pdf('$*.tex')"

%.html: %.Rmd
	$(R_HOME)/bin/Rscript -e "if (getRversion() < '3.0.0') knitr::knit2html('$*.Rmd')"

