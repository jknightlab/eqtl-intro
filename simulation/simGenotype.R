## Simulate 10 SNPs with different minor allele frequencies for 300 individuals
set.seed(100)
n <- 300
maf <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
prob <- cbind((1-maf)^2, 2*maf*(1-maf), maf^2)
geno <- apply(prob, 1, function(p) sample(0:2, size=n, prob=p, replace=TRUE))
colnames(geno) <- paste("snp", 1:ncol(geno), sep="_")
rownames(geno) <- paste("sample", 1:n, sep="_")
write.table(geno, file="output/sim_genotypes.tab", sep="\t", quote=FALSE)