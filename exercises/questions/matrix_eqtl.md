# Genome-wide eQTL analysis {.exercise}

In this set of exercises we'll use Matrix-eQTL to conduct a larger scale 
scan for SNP/gene interactions. To reduce the computing time required
the data has been restricted to chromosome 9.

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

#. Use Matrix-eQTL to carry out a *cis*/*trans* eQTL analysis. 
   i. Use a 1MB window around probes as local association region and a 
      p-value threshold of $10^{-3}$ and $10^{-5}$ for $cis$ and $trans$ associations 
      respectively.
   ii. Repeat the analysis with 10 PCs included as covariates. How do the numbers of
      reported *cis* and *trans* associations change?
   iii. Replace the Probe IDs in the Matrix-eQTL results with the corresponding
      gene names.
   iv. Find the results for SNP rs4077151 and compare them to the result from the previous 
      exercise.

   