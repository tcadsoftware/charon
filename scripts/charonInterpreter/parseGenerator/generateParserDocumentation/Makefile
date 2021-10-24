
PDFLATEX=pdflatex
MKINDEX=makeindex
BIBTEX=bibtex

BASE 		:= Charon_generateParsers

CLASS_FILES	:= ./SANDreport.cls

BIBFILE		= references  Charon_bibliography

GRAPH_FILES	= 
#./examples/pn-diode/pndiode-geom.pdf \
#./examples/pn-diode/pndiode-mesh.pdf ./misc/distrib_tree.pdf

#--------------------------------------------------------------------
TEX_ALL_FILES	:= $(wildcard *.tex)

AUX		:= $(addsuffix .aux, $(BASE))
BIB		:= $(addsuffix .bib, $(BIBFILE))
IST		:= $(addsuffix .ist, $(BASE))
BBL		:= $(addsuffix .bbl, $(BASE))
PDF		:= $(addsuffix .pdf, $(BASE))
ALL		:= $(AUX) $(BBL) $(PDF)
LOG		:= $(addsuffix .log, $(BASE))

#TEXMF		:= .:~/.tex:/usr/share/texmf
TEXINPUTS	:= .:./misc:$(TEXINPUTS)
TEX_SCR_OUT	:= make.out

pdf:	$(PDF)
all:	pdf

$(AUX):	$(TEX_ALL_FILES) $(BIB)
$(BBL):	$(TEX_ALL_FILES) $(BIB)
$(PDF):	$(TEX_ALL_FILES) $(BIB) $(CLASS_FILES)

$(PDF):	%.pdf: %.tex %.bbl $(GRAPH_FILES)
	echo "Building PDF file $(PDF_TARGET)...";\
	TEXINPUTS=$(TEXINPUTS) TEXMF=$(TEXMF) $(PDFLATEX) -shell-escape $*.tex; \
	$(BIBTEX) $(BASE) ; \
	echo "Rerun to get cross-references"; \
	TEXINPUTS=$(TEXINPUTS) TEXMF=$(TEXMF) $(PDFLATEX)  -shell-escape $*.tex; \


$(BBL): %.bbl: %.tex
	TEXINPUTS=$(TEXINPUTS) $(PDFLATEX) -shell-escape  $*.tex ; \
	$(BIBTEX) $* ; \
        if (sed -n '/citation/d' $*.aux;) then \
                $(BIBTEX) $* ; \
        fi

distclean:
	rm -f *.aux *.log *.bbl *.blg *.lof *.idx *.lot *.out *.toc *.ind *.ilg *.ent

clean:
	rm -f *.aux *.log *.bbl *.blg *.lof *.idx *.lot *.out *.toc *.ind *.ilg *.ent

clean4ht:
	rm -f *.png *.css *.html *.4ct *.4tc *.idv *.lg *.xref *.dvi *.tmp
