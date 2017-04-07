
# TeX sources of the documentation
DOC_SOURCE =	\
				tutorial-rules.tex \
				tutorial-viz.tex \
				tutorial-groups.tex

#				faq.tex


# TeX sources source files
DOC_SOURCE_SRC =	\
					mhchem.sty \
					def-include.tex \
					gml-include.tex \
					Figures/aldose-ketose-rx.cht \
					Figures/aldose-ketose-rx.eps \
					Figures/aldose-ketose-rx.fig \
					Figures/aldose-ketose-rx.pdf \
					Figures/bromination-rx.cht \
					Figures/bromination-rx.eps \
					Figures/bromination-rx.fig \
					Figures/bromination-rx.pdf \
					Figures/cannizzaro-ITS.cht \
					Figures/cannizzaro-ITS.eps \
					Figures/cannizzaro-ITS.fig \
					Figures/cannizzaro-ITS.pdf \
					Figures/cannizzaro-rx.cht \
					Figures/cannizzaro-rx.eps \
					Figures/cannizzaro-rx.fig \
					Figures/cannizzaro-rx.pdf \
					Figures/diels-alder-rx.cht \
					Figures/diels-alder-rx.eps \
					Figures/diels-alder-rx.fig \
					Figures/diels-alder-rx.pdf \
					Figures/keto-enol-rx.cht \
					Figures/keto-enol-rx.eps \
					Figures/keto-enol-rx.fig \
					Figures/keto-enol-rx.pdf \
					Figures/tbilogo.eps \
					Figures/tbilogo.pdf \
					Figures/univielogo.eps \
					Figures/univielogo.pdf \
					Figures/unifrlogo.eps \
					Figures/unifrlogo.png \
					Figures/NAD.pdf \
					Figures/NAD-groups.pdf \
					Figures/NAD-groups-new.pdf \
					Figures/NADH.pdf \
					Figures/NADH-groups.pdf \
					Figures/NADH-groups-new.pdf \
					Figures/CONH2.pdf \
					Figures/CONH2-new.pdf \
					Figures/lactat-dehydrogenase.pdf \
					Figures/depictMol.pdf \
					Figures/depictMolIndex.pdf \
					Figures/molGraph.pdf \
					Figures/molGraphIndex.pdf \
					$(DOC_SOURCE:.tex=.tex.in)


################################################################################
# IF DOCUMENTATION HAS TO BE COMPILED
if build_docs
################################################################################

# documentation target files 
DOC_TARGETS = $(DOC_SOURCE:.tex=.pdf)

# temporary files created during documentation compilation
DOC_TARGET_TMP_FILES =	$(DOC_SOURCE:.tex=.aux) \
						$(DOC_SOURCE:.tex=.bbl) \
						$(DOC_SOURCE:.tex=.blg) \
						$(DOC_SOURCE:.tex=.log) 

# files to delete with clean
MOSTLYCLEANFILES =	$(DOC_TARGETS) $(DOC_TARGET_TMP_FILES)

# compilation of TeX file via pdflatex and removal of temporary files
.tex.pdf: $(DOC_SOURCE) $(DOC_SOURCE_SRC)
	pdflatex $<
	pdflatex $<
	$(RM) $(DOC_TARGET_TMP_FILES)
	
all : $(DOC_SOURCE) $(DOC_SOURCE_SRC) $(DOC_TARGETS)

EXTRA_DIST = $(DOC_SOURCE_SRC) $(DOC_TARGETS) 

################################################################################
else
################################################################################

# DO NOTHING

EXTRA_DIST = $(DOC_SOURCE_SRC) 

################################################################################
endif
################################################################################

