BUILD_DIR=docker
SIM_DIR=simulation
SIM_DATA=$(SIM_DIR)/output
DATA_OUT=$(BUILD_DIR)/data
SIM_OUT=$(DATA_OUT)/simulated

MKDIR=mkdir -p

SIM_FILES=$(SIM_OUT)/sim_genotypes.tab \
          $(SIM_OUT)/sim_covariates.tab \
          $(SIM_OUT)/sim_expression1.tab \
          $(SIM_OUT)/sim_expression2.tab
          
all: $(SIM_FILES)

$(SIM_OUT)/sim_genotypes.tab: $(SIM_DATA)/sim_genotypes.tab
$(SIM_OUT)/sim_covariates.tab: $(SIM_DATA)/sim_covariates.tab
$(SIM_OUT)/sim_expression1.tab: $(SIM_DATA)/sim_expression1.tab
$(SIM_OUT)/sim_expression2.tab: $(SIM_DATA)/sim_expression2.tab

$(SIM_DATA)/%:
	$(MAKE) -C $(SIM_DIR)
$(SIM_OUT)/%:
	$(MKDIR) $(SIM_OUT)
	cp -f $< $@
