# latex Makefile
texpath=/Library/TeX/texbin/
PDFLATEX=TEXMFOUTPUT=../.. ${texpath}pdflatex -halt-on-error -synctex=1 --interaction=nonstopmode
BIBTEX=${texpath}bibtex


all: clean MGPS_2020B_proposal

MGPS_2020B_proposal: 
	${PDFLATEX} MGPS_2020B_proposal.tex
	${BIBTEX} MGPS_2020B_proposal
	${PDFLATEX} MGPS_2020B_proposal.tex
	${BIBTEX} MGPS_2020B_proposal
	${PDFLATEX} MGPS_2020B_proposal.tex
	bibexport -o extracted MGPS_2020B_proposal.aux
	gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER \
       -dFirstPage=1 -dLastPage=4 \
       -dMaxBitmap=500000000 \
       -dPDFSETTINGS=/prepress \
       -sOutputFile=MGPS_2020B_proposal_4pg.pdf MGPS_2020B_proposal.pdf

clean: 
	@rm -f *.aux *.bbl *.blg *.dvi *.log

