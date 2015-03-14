WD:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
UNAME=humburg
GNAME=$(UNAME)
RSCRIPT=docker run -v  $(WD):/simulation -w /simulation humburg/eqtl-intro Rscript

output/sim_genotypes.tab: simGenotype.R
	$(RSCRIPT) simGenotype.R
	chown $(UNAME):$(GNAME) output/sim_genotypes.tab