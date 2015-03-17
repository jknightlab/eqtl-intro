# Associations between SNPs and gene expression - Confounding variation {.exercise}
In this example we investigate the effect that the presence of
other sources of variation has on our ability to detect the
genotypic effects of interest.

This exercise uses the same simulated genotypes as the previous one
(`/data/simulated/sim_genotypes.tab`). The gene expression data is
located in `/data/simulated/sim_expression2.tab`. The later parts of 
the exercise also requires a number of covariates located in
`/data/simulated/sim_covariates.tab`

## Questions

#. Create a plot of gene expression by genotype for one of the SNP/gene pairs.
   How does this compare to the plot from the previous exercise?
#. Carry out a simple  eQTL analysis for the matched SNP/gene pairs.
    i. For each SNP/gene pair fit a linear regression model to obtain an 
       estimate of the genotype effect on gene expression and compute the
       95% confidence intervals for the ten SNP effects.
    ii. Create a plot that compares the estimates of effect size obtained
       above to the true value of 1.5. How does this compare to the results
       from the previous example?
#. Using the additional variables contained in the covariates file,
   fit another set of models.
    i. For each gene fit a model that incorporates the corresponding
       SNP as well as the first five variables from the covariates file.
    ii. Create the same plot of effect size estimates as before for this 
       extended model. How do they compare?
    iii. Repeat the above analysis with all covariates included in the model.
    iv. Create a plot of gene expression by genotype illustrating the effect.