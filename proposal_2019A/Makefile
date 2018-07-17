# latex Makefile
texpath=/usr/texbin/
PDFLATEX=TEXMFOUTPUT=../.. ${texpath}pdflatex -halt-on-error -synctex=1 --interaction=nonstopmode
BIBTEX=${texpath}bibtex


all: clean MGPS_2019_proposal

MGPS_2019_proposal: 
	${PDFLATEX} MGPS_2019_proposal.tex
	${BIBTEX} MGPS_2019_proposal
	${PDFLATEX} MGPS_2019_proposal.tex
	${BIBTEX} MGPS_2019_proposal
	${PDFLATEX} MGPS_2019_proposal.tex
	gs -sDEVICE=pdfwrite -dNOPAUSE -dBATCH -dSAFER \
       -dFirstPage=1 -dLastPage=4 \
       -dMaxBitmap=500000000 \
       -dPDFSETTINGS=/prepress \
       -sOutputFile=MGPS_2019_proposal_4pg.pdf MGPS_2019_proposal.pdf

clean: 
	@rm -f *.aux *.bbl *.blg *.dvi *.log

