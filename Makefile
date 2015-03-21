export TOP_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
export GENOTYPE_DIR = $(HOME)/eqtl_course/
export UNAME=humburg
export GNAME=$(UNAME)
BUILD_DIR=docker
SIM_DIR=simulation
SIM_DATA=$(SIM_DIR)/output
DATA_OUT=$(BUILD_DIR)/data
SIM_OUT=$(DATA_OUT)/simulated
EX_DIR=exercises

MKDIR=mkdir -p

SIM_FILES=$(SIM_OUT)/sim_genotypes.tab \
          $(SIM_OUT)/sim_covariates.tab \
          $(SIM_OUT)/sim_expression1.tab \
          $(SIM_OUT)/sim_expression2.tab
EX_FILES=$(EX_DIR)/exercises.md $(EX_DIR)/exercises_and_solutions.md \
         $(EX_DIR)/exercises.html $(EX_DIR)/exercises_and_solutions.html \
         $(EX_DIR)/exercises.pdf $(EX_DIR)/exercises_and_solutions.pdf

all: simulation exercises

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
	
.PHONY: exercises simulation
exercises: $(EX_FILES)
simulation: $(SIM_FILES)