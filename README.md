*Population Genomics Workshop 2016, University of Sheffield*
#Genome-wide association studies in structured populations
####Michael Stocks

The aim of this practical is to identify markers associated with male mating behaviour in the ruff. We will use the dataset of [*Küpper et al. (2016)*](http://www.nature.com/ng/journal/v48/n1/full/ng.3443.html) and attempt to identify markers underlying the "female-mimic" *Faeder* morph. For this, we will use the `R` package [GenABEL](http://genabel.org/), and by the end of the tutorial you should be able to:

* Perform a genome-wide association study
* Add covariates to the model
* Account for multiple testing
* Correct for population stratification

##Getting started

**make better**

If you type `ls` and enter in the terminal then you should see two files. `gen_RUFF_qc.raw` was generated from `VCF` format using `plink` and contains the genotypic information for each marker and individual. `phe_RUFF.txt` is a tab-delimited file containing the phenotypes for each individual.

Type `R` and enter to begin an `R` session. You can exit at any time by typing `q()` to get back to the `unix` environment.

##Importing and exploring the data
Load the library:
```{r }
library(GenABEL)
```
Now import the data and assign the object to the variable `ruff.data`:
```{r }
ruff.data <- load.gwaa.data(phe = "phe_RUFF.txt", gen = "gen_RUFF_qc.raw", force = T)
```
We'll now do some basic filtering to remove uninformative markers and to decrease the size of the dataset:
```{r }
qc <- check.marker(ruff.data, callrate = 0.66, p.level = 1e-5, perid.call = 0, extr.perid.call = 0, ibs.mrk = -1)
ruff.clean <- ruff.data[qc$idok, qc$snpok]
```
We will not go into the details of these filters but more information is available [here](http://genabel.org/GenABEL/check.marker.html).

A brief overview of the data are given for the traits or the markers using the following two commands:
```{r }
descriptives.trait(ruff.clean)
descriptives.marker(ruff.clean)
```
The data consists of 41 male ruff individuals. The phenotypic data is contained in a data frame given by the `phdata` method. Using the command `table(ruff.clean@phdata[, "morph"])`, you can see that of the three different morphs: 21 are *Independent*, 10 are *Satellite* and 10 are *Faeder*. Each of these morphs has a number of distinct physiological and behavioural phenotypes (see table below) that have been shown through breeding experiments to be Mendelian inherited traits (Lank *et al.* [1995](http://www.nature.com/nature/journal/v378/n6552/abs/378059a0.html), [2013](http://rsbl.royalsocietypublishing.org/content/9/6/20130653)). 

<p align="center">
<img src="https://github.com/mspopgen/genomics-workshop2016/blob/master/ruff-sys.png" width="640" align="center">
</p>

The genotypic data consists of markers derived from RAD sequencing, and have been mapped to the ruff draft genome. See [*Küpper et al. (2016)*](http://www.nature.com/ng/journal/v48/n1/full/ng.3443.html), and it's sister publication [*Lamichhaney et al. (2016)*](http://www.nature.com/ng/journal/v48/n1/full/ng.3430.html), for more details. After filtering we now have a dataset of 383,514 snps with which to look for associations with the trait of interest. 

##Perform an association study

We will start with some simple tests of association on just one of the traits, `fo`, which indicates whether an individual is a *Faeder* or non-*Faeder* morph. That is, are there markers associated with an individual being either a *Faeder* or a non-*Faeder* (i.e. a *Satellite* or *Independent*) morph? This trait is indicated by the phenotype `fo` and is a categorical trait, so we can start by using a logistic regression:
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

Due to the number of tests being performed (equal to the number of markers), we would expect some significant results by chance alone. There are numerous ways to do this (e.g. Bonferroni correction, FDR etc...), one which is to use permutations of the data to find a signifance cut-off threshold. This randomly shuffles the phenotypic values with respect to the individual genotypes at each marker. This creates independence between the trait and the markers and can be used to generate a suitable significance threshold. This can be done using the `qtscore` function, using the `times` option to specify 100 permutations:
```{r }
fo.QT1k <- qtscore(fo ~ 1, data = ruff.data, trait = "binomial", times = 100)
```
Using the `summary()` command, you can see that the *p*-values in the `P1df` column have been adjusted to account for multiple testing and indicate the proportion of permutations yielding a more significant *p*-value than that observed in the real data. 

##Population stratification

Structure can create artificial associations between markers and phenotypes. In natural populations, there are a number of ways to use non-causal markers to account for any population stratification. We will look at two methods:

1. The first method uses all markers in the dataset to measure the level of statification within the population. This value (often known as lambda) is then used to correct the *p*-values from the association test. The value of lamdba can be obtained with:
```{r }
lambda(ruff.MLR)
```
A value of 1 indicates no stratification, and everything above that means that there may be some form of structure in the population. The corrected *p*-values are given by the `summary()` command under the `Pc1df` heading.
2. Another method is to first estimate identity-by-state (IBS) and kinship information from the marker dataset, perform multidimensional scaling on these IBS coefficients, and then use these as covariants in the model. First, get the kinship matrix:
```{r }
ruff.kin <- ibs(ruff.clean)
```
Then perform the multidimensional scaling:
```{r }
ruff.mds <- cmdscale(as.dist(0.5 - ruff.kin))
```
And add these as covariates in the model:
```{r }
fo.QTibs <- qtscore(fo ~ mds[, 1] + mds[, 2], data = ruff.clean, trait.type = "binomial")
```
***and plot again***

##References

Aulchenko, Yurii S., Stephan Ripke, Aaron Isaacs, and Cornelia M. Van Duijn. "GenABEL: an R library for genome-wide association analysis." *Bioinformatics* 23, no. 10 (2007): 1294-1296.

Küpper, Clemens, Michael Stocks, Judith E. Risse, Natalie dos Remedios, Lindsay L. Farrell, Susan B. McRae, Tawna C. Morgan et al. "A supergene determines highly divergent male reproductive morphs in the ruff." *Nature genetics* (2016).

Lamichhaney, Sangeet, Guangyi Fan, Fredrik Widemo, Ulrika Gunnarsson, Doreen Schwochow Thalmann, Marc P. Hoeppner, Susanne Kerje et al. "Structural genomic changes underlie alternative reproductive strategies in the ruff (Philomachus pugnax)." *Nature genetics* 48, no. 1 (2016): 84-88.

Lank, David B., Constance M. Smith, Olivier Hanotte, Terry Burke, and Fred Cooke. "Genetic polymorphism for alternative mating behaviour in lekking male ruff Philomachus pugnax." *Nature* (1995).

Lank, David B., Lindsay L. Farrell, Terry Burke, Theunis Piersma, and Susan B. McRae. "A dominant allele controls development into female mimic male and diminutive female ruffs." *Biology letters* 9, no. 6 (2013): 20130653.
