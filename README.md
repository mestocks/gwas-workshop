*Population Genomics Workshop 2016, University of Sheffield*
#Genome-wide association studies in structured populations
####Michael Stocks

By the end of this tutorial you should be able to:

* Perform a genome-wide association study
* Add covariates to the model
* Account for multiple testing
* Correct for population stratification

##Getting started

If you type `ls` and enter in the terminal then you should see two files. `gen_RUFF_qc.raw` was generated from `VCF` format using `plink` and contains the genotypic information for each marker and individual. `phe_RUFF.txt` is a tab-delimited file containing the phenotypes for each individual.

Type `R` and enter to begin an `R` session. You can exit at any time by typing `q()` to get back to the `unix` environment.

##Importing and exploring the data
Load the library:
```{r }
library(GenABEL)
```
Now import the data and assign the object to the variable `ruff.data`:
```{r }
ruff.data <- load.gwaa.data(phe = "gen_RUFF_qc.raw", gen = "phe_RUFF.raw", force = T)
```
A brief overview of the data are given for the traits or the markers using the following two commands:
```{r }
descriptives.trait(ruff.data)
descriptives.marker(ruff.data)
```

<p align="center">
<img src="https://github.com/mspopgen/genomics-workshop2016/blob/master/ruff-sys.png" width="640" align="center">
</p>

See [*Küpper et al. (2016)*](http://www.nature.com/ng/journal/v48/n1/full/ng.3443.html), and it's sister publication [*Lamichhaney et al. (2016)*](http://www.nature.com/ng/journal/v48/n1/full/ng.3430.html), for more details.




##Perform an association study

First let's perform a GWAS on the *Faeder* morph. That is, are there markers associated with an individual being either a *Faeder* or a non-*Faeder* (i.e. a *Satellite* or *Independent*) morph? This trait is indicated by the phenotype `fo` and is a categorical trait, so we can start by using a logistic regression:
```{r }
fo.MLR <- mlreg(fo ~ 1, data = ruff.data, trait = "binomial")
```
This then gives us a *p*-value for association between each marker and the `fo` trait. A breakdown of the 10 most significant associations can be produced using the `summary()` command:
```{r }
summary(fo.MLR)
```
You can also increase the number of the markers displayed in the summary:
```{r }
summary(fo.MLR, top = 30)
```
This gives information about the location (`Chromosome`, `Position`) of each marker and the significance of the association (`P1df`). We will deal with the `Pc1df` column in the population stratification section (more details on the columns are given at [http://www.genabel.org/GenABEL/scan.gwaa-class.html](http://www.genabel.org/GenABEL/scan.gwaa-class.html)). You can also check how many markers are below a certain significance threshold:
```{r }
sum(fo.MLR[, "P1df"] <= 0.0001)
```
##Correct for multiple tests

Due to the number of tests being performed (equal to the number of markers), we would expect some significant results by chance alone. There are numerous ways to do this (e.g. Bonferroni correction, FDR etc...), one which is to use permutations of the data to find a signifance cut-off threshold. This randomly shuffles the phenotypic values with respect to the individual genotypes at each marker. This creates independence between the trait and the markers and can be used to generate a suitable significance threshold. This can be done using the `qtscore` function, using the `times` option to specify 1000 permutations:
```{r }
fo.QT1k <- qtscore(fo ~ 1, data = ruff.data, trait = "binomial", times = 1000)
```
Using the `summary()` command, you can see that the *p*-values in the `P1df` column have been adjusted to account for multiple testing and indicate the proportion of permutations yielding a more significant *p*-value than that observed in the real data. 

##Population stratification

Structure can create artificial associations between markers and phenotypes. In natural populations, there are a number of ways to use non-causal markers to account for any population stratification. We will look at two methods:

1. The first method uses all markers in the dataset to measure the level of statification within the population. This value (known as &lambda)

##References

