## Generate gene expression data that depends on genotype as well
## as a number of (hidden) covariates.
expr <- readr::read_tsv("output/sim_expression1.tab")
cvrt <- readr::read_tsv("output/sim_covariates.tab")
coef <- seq(2.5, 0.5, length.out=ncol(cvrt)-1)
expr2 <- cbind(expr[1], expr[-1] + colSums(coef*t(cvrt[-1])))
names(expr2) <- c("sample", paste("gene", 1:(ncol(expr2)-1), sep="_"))
write.table(expr2, file="output/sim_expression2.tab", sep="\t", quote=FALSE, row.names=FALSE)