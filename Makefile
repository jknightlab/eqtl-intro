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
GITHUB_DIR=/home/humburg/git/eqtl-intro-page

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
	
$(SLIDE_DIR)/%:
	$(MAKE) -C $(SLIDE_DIR)
	
.PHONY: data exercises slides all deploy deploy_github
data: $(SIM_FILES)
exercises: $(EX_FILES)
slides: $(SLIDE_FILES)
deploy: 
	cp $(SLIDE_FILES) $(WEB_DIR)/
	cp -r $(SLIDE_DIR)/figure $(WEB_DIR)/
	cp $(TOP_DIR)/include/slides.css $(WEB_DIR)/include/
	mv $(WEB_DIR)/eqtl-analysis.html $(WEB_DIR)/index.html
deploy_github:
	cp $(SLIDE_FILES) $(GITHUB_DIR)/slides/
	cp -r $(SLIDE_DIR)/figure $(GITHUB_DIR)/slides
	cp $(TOP_DIR)/include/slides.css $(GITHUB_DIR)/slides/include/
	cp $(EX_DIR)/exercises.html $(GITHUB_DIR)/exercises/ 
	cp $(EX_DIR)/exercises_and_solutions.html $(GITHUB_DIR)/exercises/
	cd $(GITHUB_DIR)
	git add .
	git commit -m "Updated page content"
	git push
	