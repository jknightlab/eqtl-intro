---
title: Exercises for *Introduction to eQTL Analysis*
author: Peter Humburg
---




# Prerequisites {-}
## Using the Docker image
All exercises assume the use of the docker container `humburg/eqtl-intro`
to provide the required data as well as the necessary software environment.
This requires a working Docker installation^[installation instructions are available
from the [Docker website](https://docs.docker.com/installation/).].
The docker image can be obtained from DockerHub via

```sh
docker pull humburg/eqtl-intro
```
To run the RStudio server run

```sh
docker run -p 8787:8787 humburg/eqtl-intro
```
RStudio is then accessible at `localhost:8787` or, when using
[boot2docker](http://boot2docker.io/) via the IP address indicated 
by `boot2docker ip`.

## Included data
The image includes a number of simulated and real data sets used for these 
exercises. All data are provided as tab-separated files (typically with a column header).
Files are located in directories below `/data`. All simulated data are located in
`/data/simulated`. Real data can be found in `/data/genotyping`, `/data/expression`
and `/data/annotation` for genotyping, gene expression and annotation data 
respectively.


# Associations between SNPs and gene expression - A simple example {.exercise}
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


