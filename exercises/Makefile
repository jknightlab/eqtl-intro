WD:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
RSCRIPT=docker run -v  $(WD):/exercises -v $(GENOTYPE_DIR):/data/genotypes \
        -v $(ANALYSIS_DIR):/data/analysis -w /exercises humburg/eqtl-intro Rscript
PANDOC=docker run -v  $(WD):/exercises -v $(TOP_DIR)include/:/exercises/include \
       -w /exercises humburg/eqtl-intro pandoc

PANDOC_COMMON = --standalone --smart --normalize --toc --highlight-style=tango \
             --bibliography=include/references.bib
PANDOC_PDF = -t latex --latex-engine=pdflatex --template=include/report.latex \
           -V geometry:margin=2cm -V documentclass:article -V classoption:a4paper \
           -H include/captions.tex
PANDOC_HTML = -t html5 --self-contained --webtex --template=include/report.html5 \
            -c include/buttondown.css


EXERCISES=questions/simple_eqtl.md questions/covar_eqtl.md questions/pca_eqtl.md
SOLUTIONS=solutions/simple_eqtl.Rmd solutions/covar_eqtl.Rmd solutions//pca_eqtl.Rmd

all: exercises.md exercises_and_solutions.md pdf html
exercises.md: exercises.Rmd $(EXERCISES)
	$(RSCRIPT) exercises.r
exercises_and_solutions.md: exercises.Rmd $(EXERCISES) $(SOLUTIONS)
	$(RSCRIPT) exercises.r --solution
exercises.html: exercises.md ../include/buttondown.css ../include/report.html5 ../include/references.bib
	$(PANDOC) $(PANDOC_COMMON) $(PANDOC_HTML) -s exercises.md -o exercises.html
	chown $(UNAME):$(GNAME) exercises.html
exercises_and_solutions.html: exercises_and_solutions.md ../include/buttondown.css ../include/report.html5 ../include/references.bib
	$(PANDOC) $(PANDOC_COMMON) $(PANDOC_HTML) -s exercises_and_solutions.md -o exercises_and_solutions.html
	chown $(UNAME):$(GNAME) exercises_and_solutions.html
exercises.pdf: exercises.md ../include/report.latex ../include/captions.tex ../include/references.bib
	$(PANDOC) $(PANDOC_COMMON) $(PANDOC_PDF) -s exercises.md -o exercises.pdf
	chown $(UNAME):$(GNAME) exercises.pdf
exercises_and_solutions.pdf: exercises_and_solutions.md ../include/report.latex ../include/captions.tex ../include/references.bib
	$(PANDOC) $(PANDOC_COMMON) $(PANDOC_PDF) -s exercises_and_solutions.md -o exercises_and_solutions.pdf
	chown $(UNAME):$(GNAME) exercises_and_solutions.pdf


.PHONY: pdf html all
pdf: exercises.pdf exercises_and_solutions.pdf
html: exercises.html exercises_and_solutions.html