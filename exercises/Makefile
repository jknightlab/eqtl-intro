WD:=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
RSCRIPT=docker run -v  $(WD):/exercises -w /exercises humburg/eqtl-intro Rscript

EXERCISES=questions/simple_eqtl.md questions/covar_eqtl.md
SOLUTIONS=solutions/simple_eqtl.Rmd solutions/covar_eqtl.Rmd

all: exercises.md exercises_and_solutions.md
exercises.md: exercises.Rmd $(EXERCISES)
	$(RSCRIPT) exercises.r
exercises_and_solutions.md: exercises.Rmd $(EXERCISES) $(SOLUTIONS)
	$(RSCRIPT) exercises.r --solution
