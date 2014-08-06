all: NEWS.pdf

NEWS.pdf: NEWS.Rd
	R CMD Rd2pdf --force NEWS.Rd
clean:
	rm -rf NEWS.pdf .Rd2pdf*/
