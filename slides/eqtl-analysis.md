---
title: Introduction to eQTL analysis
author: Peter Humburg
---

# *Hands-on* : Setup
## Docker image
Make sure you have the latest version of the docker image.

```sh
docker pull humburg/eqtl-intro
```

## Directory structure
* Create a working directory for this course. (`eqtl_course`)
* Within that directory, create two sub-directories
    * `genotypes`
    * `analysis`
* *Mac*: create directories under `/Users/`
* *Windows*: create directories under `C:\Users\`

## Additional data
The docker image contains most of the data. Further genotyping data is available from

```sh
ftp://galahad.well.ox.ac.uk
```
Username: eqtl_course

Download the data and save it into `genotypes`

<div class="notes">
Need to provide password separately.
</div>

## Start the RStudio server    

Windows
  :
  
    docker run -p 8787:8787 
      -v /c/Users/user/eqtl_course/genotpes:/data/genotypes
      -v /c/Users/user/eqtl_course/analysis:/data/analysis 
      humburg/eqtl-intro
    

IP address of the server usually is 192.168.59.103.

Use `boot2docker ip` to check if necessary. 

## Start the RStudio server

Mac
  :
  
    
    docker run -p 8787:8787 
       -v /Users/user/eqtl_course/genotpes:/data/genotypes 
       -v /Users/user/eqtl_course/analysis:/data/analysis 
       humburg/eqtl-intro

IP address of the server usually is 192.168.59.103.

Use `boot2docker ip` to check if necessary. 

## Start the RStudio server    

Linux
  :
  
    docker run -p 8787:8787 
      -v /home/user/eqtl_course/genotpes:/data/genotypes
      -v /home/user/eqtl_course/analysis:/data/analysis
      -e USER=$USER  -e USERID = $UID 
      humburg/eqtl-intro
    

IP address of the server usually is 127.0.0.1 (or *localhost*).
    
## Using RStudio

Access the RStudio interface at *http://yourip:8787*. 

    * Username: rstudio
    * Password: rstudio

# What are eQTL?
## Quantitative trait loci
QTL are regions of the genome associated with quantitative traits

* height
* BMI
* lung capacity
* ...

## Expression quantitative trait loci {data-transition="none"}
<div class="multicolumn">
<div>
If the trait of interest is the expression of a gene, we talk about eQTL.

* **Associations can be local (*cis*)**
* or distant (*trans*)

</div>
<div>
![](figure/cis_eqtl.png)
</div>
</div> 

## Expression quantitative trait loci {data-transition="none"}
<div class="multicolumn">
<div>
If the trait of interest is the expression of a gene, we talk about eQTL.

* Associations can be local (*cis*)
* **or distant (*trans*)**

</div>
<div>
![](figure/trans_eqtl.png)
</div>
</div> 

## Expression quantitative trait loci {data-transition="none"}
<div class="multicolumn">
<div>
If the trait of interest is the expression of a gene, we talk about eQTL.

* Associations can be local (*cis*)
* **or distant (*trans*)**

</div>
<div>
![](figure/trans_network.png)
</div>
</div> 


# Detecting eQTL -- A simple model
## How do we know whether a locus is associated with the expression of a gene? {data-transition="none"}

<div class="multicolumn">
<div>
![](figure/genomic_locus.png)

</div>
<div>
* **Any given locus contains multiple SNPs.**
* Can determine genotypes for these SNPs in a (large) number of individuals.
* Measure gene expression for genes of interest.
* Assess the evidence that expression varies with genotype.
</div> 
</div>

## How do we know whether a locus is associated with the expression of a gene? {data-transition="none"}

<div class="multicolumn">
<div>

 sample	   snp_1   snp_2   snp_3   ...
--------   -----   -----   -----  -----
sample 1    AA      AA      AB     ...
sample 2    AB      AB      AA     ...
sample 3    AB      BB      AB     ...
...        ...     ...      ...    ...

</div>
<div>
* Any given locus contains multiple SNPs.
* **Can determine genotypes for these SNPs in a (large) number of individuals.**
* Measure gene expression for genes of interest.
* Assess the evidence that expression varies with genotype.
</div> 
</div>

## How do we know whether a locus is associated with the expression of a gene? {data-transition="none"}

<div class="multicolumn">
<div>

 sample	   gene_1   gene_2   gene_3   ...
--------   ------  ------    ------  -----
sample 1    7.3     12.8      6.5     ...
sample 2    10.9    9.6       8.8     ...
sample 3    9.5     10.7      15.1    ...
...         ...     ...       ...     ...

</div>
<div>
* Any given locus contains multiple SNPs.
* Can determine genotypes for these SNPs in a (large) number of individuals.
* **Measure gene expression for genes of interest.**
* Assess the evidence that expression varies with genotype.
</div> 
</div>

## How do we know whether a locus is associated with the expression of a gene? {data-transition="none"}

<div class="multicolumn">
<div id="snp-expr">
![](figure/snp_expression.png)
</div>
<div>
* Any given locus contains multiple SNPs.
* Can determine genotypes for these SNPs in a (large) number of individuals.
* Measure gene expression for genes of interest.
* **Assess the evidence that expression varies with genotype.**
</div> 
</div>


## Linear additive model
<div class="multicolumn">
<div id="snp-expr">
![](figure/snp_expression.png)
</div>
<div>
Different alleles of a SNP may exhibit a dosage effect.

* Using the *AA* genotype as baseline, each copy of the *B* allele changes 
  expression by a fixed amount.
* Implies a linear relationship between the mean
  gene expression and the number of *B* alleles.
* Estimating the change in expression due to the *B* allele to quantify 
  the SNP's contribution to gene expression. 
</div> 
</div>

## Linear regression {data-transition="none"}
$$Y = \beta_0 + \beta X + \varepsilon$$

Y
  : Response variable (here: vector of expression values for gene of interest).

X
  : Explanatory variable (here: vector of genotypes (coded as 0, 1, 2) 
    for the SNP under consideration).

## Linear regression {data-transition="none"}
$$Y = \beta_0 + \beta X + \varepsilon$$

$\beta_0$
  : Intercept (here: mean expression for *AA* genotype).

$\beta$
  : Regression coefficient; the effect of X on the mean of Y (here: change in 
    mean gene expression for each copy of the *B* allele).

$\varepsilon$
  : Residuals; the difference between observed values of $Y$ and the estimated
    mean of $Y|X$
    
## Linear regression -- Assumptions

Residuals are

* independent
* normally distributed with mean 0 and constant variance

. . .

It is implied that

* There is a linear relationship between $X$ and $Y$.
* Values of $Y$ for each value of $X$ are normally distributed.
* There is only one source of variation not explained by $X$.


## Linear regression -- Robustness {data-transition="none"}
### Independence of residuals
Lack of independence can produce misleading results.

What could cause this?

## Linear regression -- Robustness {data-transition="none"}
### Constant variance of residuals (homoskedacity)
Violation of this assumption will lead to incorrect p-values and confidence intervals.

## Linear regression -- Robustness {data-transition="none"}
### Normality
* Estimates and their confidence intervals and p-values are fairly robust.
* But beware of long tailed distributions.
* Prediction can become problematic (but we are not interested in that here).

## Linear regression -- Robustness {data-transition="none"}
### Linearity
<div class="left">
If the true relationship between $Y$ and $X$ is non-linear conclusions may be
misleading.

When might this occur with eQTL data?
</div>

## Linear regression -- Robustness {data-transition="none"}
### Linearity
<div class="left">
If the true relationship between $Y$ and $X$ is non-linear conclusions may be
misleading.

When might this occur with eQTL data?
</div>
<div class="right">
![](figure/dominant.png)
</div>

# *Hands-on* : Simple linear regression    
## Data

Genotypes
  : */data/simulated/sim_genotypes.tab*

Gene expression
  : */data/simulated/sim_expression1.tab*
    
## Warm-up

> #. Load the data into R.
> #. How is the data formatted?
> #. Determine minor allele frequencies.

<div class="notes">

```r
geno <- readr::read_tsv("/data/simulated/sim_genotypes.tab")
expr <- readr::read_tsv("/data/simulated/sim_expression1.tab")

maf <- colMeans(geno[-1])/2
```
</div>

##  Plotting the data
* Choose a SNP/gene pair (snp_1 / gene_1, snp_2 / gene_2, ...)
* Create a plot showing gene expression by genotype for this pair.

<div class="notes">
Assign SNP/gene pairs to participants to ensure each is handled at least once.

```r
library(ggplot2)
genoLong <- tidyr::gather(geno, snp, genotype, -sample)
exprLong <- tidyr::gather(expr, gene, expression, -sample)
dataLong <- cbind(genoLong, exprLong["expression"])
dataLong$genotype <- as.factor(dataLong$genotype) 
ggplot(dataLong, aes(genotype, expression)) +
		geom_jitter(colour="darkgrey", 
				position=position_jitter(width=0.25)) +
		geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") + 
		facet_wrap(~snp) + theme_bw()
```
</div>

## Fitting a simple linear regression

#. Fit a simple linear regression for the SNP/gene pair of your choice.
#. Compute the 95% confidence interval for the genotype effect.

<div class="notes">
Collect results from participants into a file and load into R.

Code fit all models is below.

```r
fit <- mapply(function(e, g) lm(e ~ g), 
		expr[-1], geno[-1], SIMPLIFY=FALSE)
betaHat <- sapply(fit, coef)[2,]

ci <- sapply(fit, confint, "g")
rownames(ci) <- c("lower", "upper")
```
</div>

## Examining the results

#. Create diagnostic plots.
#. Plot estimated genotype effects by minor allele frequency.

<div class="notes">
Display R's build-in diagnostics by plotting the `fit` object.

The true SNP effect is 1.5

```r
estimates <- data.frame(estimate=betaHat, t(ci), maf=maf)
fig <- ggplot(estimates, aes(x=maf)) +  
		geom_hline(yintercept=0, linetype="longdash") + 
		geom_errorbar(aes(ymin=lower, ymax=upper)) +
		geom_point(aes(y=estimate))  + theme_bw()
fig <- fig + geom_hline(yintercept=1.5)
```
</div>

# Detecting eQTL -- Not *that* simple
## Additional sources of variation

* Gene expression is subject to many sources of variation.
* A model with a single SNP as explanatory variable is unlikely to be sufficient.

. . .

Use multiple regression to obtain better estimates of SNP effects.

## Multiple regression
$$Y = \beta_0 + \sum_{i=1}^n \beta_i X_i + \varepsilon$$ 

* Similar to simple linear regression but incorporates multiple 
  explanatory variables.
* $\beta_i$ are interpreted as "change of $Y$ due to a unit change in $X_i$ when
  *all other explanatory variables are held constant*".
  
> * Explanatory variables are assumed to be uncorrelated with each other.

# *Hands-on* : Covariates
## Data

Genotypes
  : */data/simulated/sim_genotypes.tab*

Gene expression
  : */data/simulated/sim_expression2.tab*
  
Covariates
  : */data/simulated/sim_covariates.tab*
  
<div class="notes">
geno <- readr::read_tsv("/data/simulated/sim_genotypes.tab")
expr <- readr::read_tsv("/data/simulated/sim_expression2.tab")
covar <- readr::read_tsv("/data/simulated/sim_covariates.tab")
</div>

## Plotting the data
Create a plot of gene expression by genotype for you SNP/gene pair of choice.

How does this compare to the plot from the previous exercise.

<div class="notes">

```r
library(ggplot2)
genoLong <- tidyr::gather(geno, snp, genotype, -sample)
exprLong <- tidyr::gather(expr, gene, expression, -sample)
dataLong <- cbind(genoLong, exprLong["expression"])
dataLong$genotype <- as.factor(dataLong$genotype) 
ggplot(dataLong, aes(genotype, expression)) +
		geom_jitter(colour="darkgrey", 
				position=position_jitter(width=0.25)) +
		geom_boxplot(outlier.size=0, alpha=0.6, fill="grey") + 
		facet_wrap(~snp) + theme_bw()
```

There is very little evidence of a SNP effect in these plots.
</div>
    
## Simple linear regression
Repeat the simple linear regession analysis with these data.

#. Fit a model for one of the SNP/gene pairs
#. Create diagnostic plots
#. Compute the confidence interval for the coefficient

How does this compare to the result from the previous analysis?

<div class="notes">

```r
simpleFit <- mapply(function(e, g) lm(e ~ g), 
		expr[-1], geno[-1], SIMPLIFY=FALSE)
simpleBetaHat <- sapply(simpleFit, coef)[2,]
simpleCI <- sapply(simpleFit, confint, "g")
rownames(simpleCI) <- c("lower", "upper")
```

Plot all results:

```r
maf <- colMeans(geno[-1])/2
estimates <- data.frame(estimate=simpleBetaHat, t(simpleCI), maf=maf)
ggplot(estimates, aes(x=maf)) + geom_hline(yintercept=1.5) + 
		geom_hline(yintercept=0, linetype="longdash") + 
		geom_errorbar(aes(ymin=lower, ymax=upper)) +
		geom_point(aes(y=estimate))  + theme_bw()
```
Note that the true SNP coefficient is still 1.5 but estimates are
considerably worse.
</div>

## Multiple linear regression
#. Use the first five variables contained in the covariates file as 
  covariates in your model.
#. Create diagnostic plots
#. Compute the confidence interval for the coefficient

<div class="notes">

```r
covarFit <- mapply(function(e, g, var) lm(e ~ g + var), 
		expr[-1], geno[-1], 
		MoreArgs=list(as.matrix(covar[2:6])), SIMPLIFY=FALSE)
covarBetaHat <- sapply(covarFit, coef)[2,]
covarCI <- sapply(covarFit, confint, "g")
rownames(covarCI) <- c("lower", "upper")
```

May also want to try including all 20 variables.
</div>

# Covariates -- Choose wisely
## (Multi)-colinearity
If $X_i$ and $X_j$ are correlated the estimates of $\beta_i$ and
$\beta_j$ will be biased. Several issues may occur:

* All the variance in $Y$ due to $X_i$ and $X_j$ is wholly attributed to
  one of the variables (say, $X_i$). 
    
    * Results in over estimation of $\beta_i$ and under estimation of $\beta_j$
* Estimates for both $\beta_i$ and $\beta_j$ are too low because all of the variance
  is explained by the other variable that is *held constant*.
* No sensible interpretation for $\beta_i$ and $\beta_j$.
* Model fitting may fail.

## Dealing with colinearity

* Check correlation between explanatory variables.
* Check variance inflation factor (VIF).
* Choose only one from each group of correlated variables.

Only really need to worry about variables of interest for
downstream analysis.

# If only we knew -- Covariates for real data
## Limited information

* Real data comes with varying amounts of additional information.
    * Sex and age are common
    * May have detailed phenotyping data
* Some variables will have big impact on gene expression
* Won't have information on all relevant variables.

## Use the data

* We don't have to know what the (non-genetic) sources of variation are
  as long as we can account for them.
* Dimensionality reduction techniques can identify major directions of variation
  from the data.
  
## Principle component analysis

* Transforms the data into set of (linearly) uncorrelated variables (principle components).
* Principle components (PCs) are ordered by the proportion of variance they explain.
* Can reduce the dimensionality of a dataset by considering only the first $k$ PCs.

How does that help us?

## Accounting for unknown sources of variation

* The major sources of variation in gene expression data are (usually) not genetic.
* Can remove non-genetic sources of variation by including the3 first $k$ PCs into
  the model.

## Practical considerations

* Not obvious how many PCs to include in model.
* Need to look out for PCs that *do* correlate with genotype.

. . .

* Beware of data formatting issues.
* Data should be centred and scaled prior to PCA.

# *Hands-on* : Dealing with real data

## Data

Gene expression
  : */data/monocytes/expression/ifn_expression.tab.gz*
  
Genotypes
  : */data/genotypes/genotypes.tab.gz*  
  (downloaded earlier) 

Subset of data published in  
Fairfax, Humburg, Makino, *et al.*  
**Innate Immune Activity Conditions the Effect of Regulatory Variants upon Monocyte 
Gene Expression**. Science (2014). doi:[ 10.1126/science.1246949](http://doi.org/10.1126/science.1246949).

# *Hands-on* : Scaling it up
 
# Interpreting results
