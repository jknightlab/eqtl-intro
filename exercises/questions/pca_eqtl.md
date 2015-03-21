# Using principle components as covariates {.exercise}
We will explore the use of principle components as 
covariates in linear models of gene expression to account
for unknown sources of variation.

Gene expression data are located in */data/monocytes/expression/ifn_expression.tab.gz*
Genotypes are located in */data/genotypes/genotypes.tab.gz* (provided during course)

These data are part of the dataset published in
Fairfax, Humburg, Makino, et al.
Innate Immune Activity Conditions the Effect of Regulatory Variants upon 
Monocyte Gene Expression. Science (2014).
doi:[10.1126/science.1246949](http://doi.org/10.1126/science.1246949). 

In addition to the primary datasets a few files with annotations for SNPs and
genes is available in the */data/monocytes/annotation* directory:

snp_loc_hg19.tab
  : Genomic location of SNPs.

probe_loc_hg19.tab
  : Genomic location of gene expression probes.
  
probeAnnotations.tab
  : Further annotations for gene expression probes, including associated gene symbols.
  
All coordinates refer to the hg19 reference build.

## Exercises
#. Determine the dimensions of this dataset. How many genes, SNPs and samples are included?
#. Principle components of the expression data.
    i. Compute the principle components.
    ii. Create a plot of the variances for the first 20 PCs.
    iii. How much of the total variance is explained by the first 20 PCs?
#. 