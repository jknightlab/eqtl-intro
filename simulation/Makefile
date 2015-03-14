WD:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
vpath %.tab output
UNAME=humburg
GNAME=$(UNAME)
RSCRIPT=docker run -v  $(WD):/simulation -w /simulation humburg/eqtl-intro Rscript

all: sim_genotypes.tab sim_expression1.tab 
output/sim_genotypes.tab: simGenotype.R
	$(RSCRIPT) simGenotype.R
	chown $(UNAME):$(GNAME) output/sim_genotypes.tab
output/sim_expression1.tab: simExpression1.R
	$(RSCRIPT) simExpression1.R
	chown $(UNAME):$(GNAME) output/sim_expression1.tab