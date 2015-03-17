#!/usr/bin/env Rscript

library(knitr)
library(optparse)

optList <- list(make_option("--input", default="exercises.Rmd", help="Input file name."),
		make_option("--solution", default=FALSE, action="store_true",
				help="Include solutions in output."))

parser <- OptionParser("%prog [options]", option_list=optList)
opt <- parse_args(parser)

solution <- opt$solution
output <- "exercises"
if(solution){
	output <- paste0(output, "_and_solutions")
}
output <- paste0(output, ".md")
knit(opt$input, output)