# Makefile to use knitr for package vignettes
# (modified from knitr package source)
#
# DO NOT EXCLUDE FROM BUILD!

TEXS = fftw.tex normalization.tex psd_overview.tex multivariate.tex
PDFS = fftw.pdf normalization.pdf psd_overview.pdf multivariate.pdf

all: $(PDFS)

clean:
	rm -rf *-concordance.tex *.bbl *.blg *.aux *.out *.log *.spl *.md figure/ *.idx *.ilg *.ind *.toc
	rm -f *.synctex.gz
	rm -f $(TEXS)

allclean: clean
	rm -f $(PDFS)
	rm -f figure/*.tex
	rm -f figure/*.pdf

refresh: allclean all

%.pdf: %.Rnw
	$(R_HOME)/bin/Rscript -e "if (getRversion() >= '3.0.0') knitr::knit2pdf('$*.Rnw') else tools::texi2pdf('$*.tex')"

%.html: %.Rmd
	$(R_HOME)/bin/Rscript -e "if (getRversion() >= '3.0.0') knitr::knit2html('$*.Rmd')"

