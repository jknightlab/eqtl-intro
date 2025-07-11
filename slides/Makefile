SHELL := /bin/bash
PATH = $(realpath $(lastword $(MAKEFILE_LIST)))
WD:=$(shell dirname '$(PATH)')
DOCKER = /usr/bin/docker
CHOWN = /usr/bin/chown
PANDOC=$(DOCKER) run -v  "$(WD)":/slides -v '$(TOP_DIR)/include/':/slides/include \
       -w /slides humburg/eqtl-intro pandoc

PANDOC_COMMON = --standalone --toc --highlight-style=tango \
             --bibliography=include/references.bib
PANDOC_SLIDES = -t revealjs --css include/slides.css --toc-depth=1 -V toc-header="Overview" \
                --template=include/default.revealjs --slide-level=2 --mathjax

all: slides              
eqtl-analysis.html: eqtl-analysis.md ../include/slides.css ../include/default.revealjs figure/*
	$(PANDOC) $(PANDOC_COMMON) $(PANDOC_SLIDES) -s eqtl-analysis.md -o eqtl-analysis.html
	$(CHOWN) $(UNAME):$(GNAME) eqtl-analysis.html

.PHONY: slides all
slides: eqtl-analysis.html