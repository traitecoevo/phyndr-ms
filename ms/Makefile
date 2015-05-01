all: phyndr-ms.pdf

phyndr-ms.pdf: *.tex *.bib
	pdflatex phyndr-ms
	bibtex phyndr-ms
	pdflatex phyndr-ms
	pdflatex phyndr-ms

tidy:
	rm -f *.aux *.log *.bbl *.blg

clean: tidy
	rm -f phyndr-ms.pdf
