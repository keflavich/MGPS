# latex Makefile
ifndef texpath
texpath=/Library/TeX/texbin/
endif
PDFLATEX=${texpath}pdflatex -halt-on-error -synctex=1 --interaction=nonstopmode
SKIPERR=${texpath}pdflatex --interaction=nonstopmode
LATEX=${PDFLATEX}
BIBTEX=bibtex
DVIPS=dvips
PS2PDF=ps2pdf
SHELL=/bin/bash

all: gitstuff.tex pilotpaper

.PHONY: pilotpaper
pilotpaper: 
	echo "texpath: ${texpath}"
	python make.py --texpath=${texpath}
	#python make_apjform.py
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=pilotpaper_compressed.pdf pilotpaper.pdf
	python paper2arxiv.py --apj
	python paper2arxiv.py
	#python make_apjform.py

.PHONY: diff
diff:
	python parse_macros.py pilotpaper.tex pilotpaper_fordiff.tex
	#python parse_macros.py pilotpaper_feb13.tex submitted_fordiff.tex
	latexdiff submitted_fordiff.tex pilotpaper_fordiff.tex  > diff.tex
	${SKIPERR} diff.tex
	${BIBTEX} diff
	${SKIPERR} diff.tex
	${BIBTEX} diff
	${SKIPERR} diff.tex
	gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=diff_compressed.pdf diff.pdf

.PHONY: referee
referee:
	echo "texpath: ${texpath}"
	python make.py --referee --texpath=${texpath}


gitstuff.tex:
	echo "%%% This file is generated by the Makefile." > gitstuff.tex
	git log -1 --date=short --format="format:\\newcommand{\\githash}{%h}\\newcommand{\\gitdate}{%ad\\xspace}\\newcommand{\\gitauthor}{%an\\xspace}" >> gitstuff.tex
