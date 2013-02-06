# Makefile to use knitr for package vignettes

R_HOME = /usr

PDFS = fftw.pdf normalization.pdf rlpSpec_overview.pdf

all: $(PDFS)

clean:
	rm -rf *-concordance.tex *.bbl *.blg *.aux *.out *.log *.spl *.md figure/ *.idx *.ilg *.ind *.toc

%.pdf: %.Rnw
	$(R_HOME)/bin/R -q -e "library('knitr');knit('$*.Rnw')"
	$(R_HOME)/bin/R CMD pdflatex $*.tex
	$(R_HOME)/bin/R CMD bibtex $*
	$(R_HOME)/bin/R CMD makeindex $*.idx
	$(R_HOME)/bin/R CMD pdflatex $*.tex
	$(R_HOME)/bin/R CMD pdflatex $*.tex
