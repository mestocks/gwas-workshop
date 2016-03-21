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
<p align="center">
<img src="https://github.com/mspopgen/genomics-workshop2016/blob/master/ruff-sys.png" width="640" align="center">
</p>

See [*Küpper et al. (2016)*](http://www.nature.com/ng/journal/v48/n1/full/ng.3443.html), and it's sister publication [*Lamichhaney et al. (2016)*](http://www.nature.com/ng/journal/v48/n1/full/ng.3430.html), for more details.

```{r }
library(GenABEL)
```

```{r }
ruff.data <- load.gwaa.data(phe = "gen_RUFF_qc.raw", gen = "phe_RUFF.raw", force = T)
```

```{r }
descriptives.trait(ruff.data)
descriptives.trait(ruff.data)
```


First let's perform a GWAS on the *Faeder* morph. That is, are there markers associated with an individual being either a *Faeder* or a non-*Faeder* (i.e. a *Satellite* or *Independent*) morph? This categorical trait is indicated by the phenotype `fo`:
```{r }
fo.MLR <- mlreg(fo ~ 1, data = ruff.data, trait = "binomial")
```

```{r }
summary(fo.MLR)

summary(fo.MLR, top = 30)
```

```{r }
sum(fo.MLR[, "P1df"] <= 0.0001)
```
##References

