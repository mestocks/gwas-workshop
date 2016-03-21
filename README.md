*Population Genomics Workshop 2016, University of Sheffield*

#Genomewide association studies in structured populations
###Michael Stocks


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
