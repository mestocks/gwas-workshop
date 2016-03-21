# genomics-workshop2016



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
