## Simulate gene expression for one gene per SNP.
## Expression is base + snp + noise
set.seed(101)
snp <- readr::read_tsv("output/sim_genotypes.tab")
errorMeans <- rnorm(nrow(snp))
expr <- sapply(snp[-1], function(g) 7 + 1.5*g) + 
		matrix(rnorm(nrow(snp)*(ncol(snp)-1), mean=errorMeans, sd=0.5), nrow=nrow(snp))
expr <- data.frame(snp$sample, expr)
names(expr) <- c("sample", paste("gene", 1:(ncol(expr)-1), sep="_"))
write.table(expr, file="output/sim_expression1.tab", sep="\t", quote=FALSE, row.names=FALSE)