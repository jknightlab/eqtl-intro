export TOP_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
export GENOTYPE_DIR = $(HOME)/eqtl_course/genotypes
export ANALYSIS_DIR = $(HOME)/eqtl_course/analysis
export UNAME=humburg
export GNAME=$(UNAME)
BUILD_DIR=docker
SIM_DIR=simulation
SIM_DATA=$(SIM_DIR)/output
DATA_OUT=$(BUILD_DIR)/data
SIM_OUT=$(DATA_OUT)/simulated
EX_DIR=exercises
SLIDE_DIR=slides
WEB_DIR=/usr/local/www/data/eqtl-intro

MKDIR=mkdir -p

SIM_FILES=$(SIM_OUT)/sim_genotypes.tab \
          $(SIM_OUT)/sim_covariates.tab \
          $(SIM_OUT)/sim_expression1.tab \
          $(SIM_OUT)/sim_expression2.tab
EX_FILES=$(EX_DIR)/exercises.md $(EX_DIR)/exercises_and_solutions.md \
         $(EX_DIR)/exercises.html $(EX_DIR)/exercises_and_solutions.html \
         $(EX_DIR)/exercises.pdf $(EX_DIR)/exercises_and_solutions.pdf
SLIDE_FILES=$(SLIDE_DIR)/eqtl-analysis.html

all: data exercises slides deploy

$(SIM_OUT)/sim_genotypes.tab: $(SIM_DATA)/sim_genotypes.tab
$(SIM_OUT)/sim_covariates.tab: $(SIM_DATA)/sim_covariates.tab
$(SIM_OUT)/sim_expression1.tab: $(SIM_DATA)/sim_expression1.tab
$(SIM_OUT)/sim_expression2.tab: $(SIM_DATA)/sim_expression2.tab

$(SIM_DATA)/%:
	$(MAKE) -C $(SIM_DIR)
$(SIM_OUT)/%:
	$(MKDIR) $(SIM_OUT)
	cp -f $< $@
	chown -R $(UNAME):$(GNAME) $(SIM_OUT)

$(EX_DIR)/%:
	$(MAKE) -C $(EX_DIR)
	chown $(UNAME):$(GNAME) $(EX_FILES)
	
<<<<<<< HEAD
.PHONY: exercises simulation
exercises: $(EX_FILES)
simulation: $(SIM_FILES)
=======
$(SLIDE_DIR)/%:
	$(MAKE) -C $(SLIDE_DIR)
	
.PHONY: data exercises slides all deploy
data: $(SIM_FILES)
exercises: $(EX_FILES)
slides: $(SLIDE_FILES)
deploy: 
	cp $(SLIDE_FILES) $(WEB_DIR)/
	cp -r $(SLIDE_DIR)/figure $(WEB_DIR)/
	cp $(TOP_DIR)/include/slides.css $(WEB_DIR)/include/
	mv $(WEB_DIR)/eqtl-analysis.html $(WEB_DIR)/index.html
>>>>>>> master
