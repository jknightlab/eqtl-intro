# Associations between SNPs and gene expression - A simple example
We will investigate the properties of a small simulated data set 
consisting of genotypes for 10 SNPs and expression values for 10 genes
from 300 individuals. Genotypes are encoded as 0, 1, and 2, indicating 
the number of copies of the second allele.

Genotypes are located in the file `/data/simulated/sim_genotypes.tab` and
gene expression values can be found in `/data/simulated/expression1.tab`.

## Questions

#. What are the minor allele frequencies of the different SNPs in the data set?
#. Consider pairs of SNPs and genes such that *snp_1* is paired with *gene_1*,
*snp_2* with *gene_2* and so on.
	i. Create a plot showing gene expression by genotype for one of the SNP/gene pairs.
    ii. For each SNP/gene pair fit a linear regression model to obtain an 
      estimate of the genotype effect on gene expression and compute the
      95% confidence intervals for the ten SNP effects.
    iii. Create a plot to compare the estimated coefficients and their 95% confidence
      intervals to 1.5, the true value of $\beta$. What do you observe?