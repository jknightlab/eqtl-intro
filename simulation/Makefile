WD:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
vpath %.tab output
UNAME=humburg
GNAME=$(UNAME)
RSCRIPT=docker run -v  $(WD):/simulation -w /simulation humburg/eqtl-intro Rscript

all: sim_genotypes.tab sim_expression1.tab sim_covariates.tab sim_expression2.tab
output/sim_genotypes.tab: simGenotype.R
	$(RSCRIPT) simGenotype.R
	chown $(UNAME):$(GNAME) output/sim_genotypes.tab
output/sim_expression1.tab: simExpression1.R sim_genotypes.tab
	$(RSCRIPT) simExpression1.R
	chown $(UNAME):$(GNAME) output/sim_expression1.tab
output/sim_covariates.tab: simCovariate.R
	$(RSCRIPT) simCovariate.R
	chown $(UNAME):$(GNAME) output/sim_covariates.tab
output/sim_expression2.tab: simExpression2.R sim_expression1.tab sim_covariates.tab
	$(RSCRIPT) simExpression2.R
	chown $(UNAME):$(GNAME) output/sim_expression2.tab