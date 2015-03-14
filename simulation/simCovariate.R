## Generate a set of covariates to use when generating
## gene expression data.
set.seed(102)
n <- 300; nc <- 20
cvrt <- matrix(rnorm(n*nc, mean=1), ncol = nc)
cvrt <- data.frame(paste("sample", 1:n, sep="_"), cvrt)
names(cvrt) <- c("sample", paste("var", 1:nc, sep="_"))
write.table(cvrt, file="output/sim_covariates.tab", sep="\t", quote=FALSE, row.names=FALSE)