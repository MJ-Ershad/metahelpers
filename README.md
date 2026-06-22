# metahelpers

Convenience wrappers over [`metafor`](https://www.metafor-project.org/) and
[`clubSandwich`](https://cran.r-project.org/package=clubSandwich) for the
recurring case of **effect sizes nested within clusters** — multiple estimates
per study or per cohort — fit as a **three-level random-effects** meta-analysis
with cluster-robust inference.

The point is to turn the boilerplate of setting up `rma.mv(random = ~ 1 | cluster/effect)`,
computing CR2 standard errors, decomposing `I^2`, and building a prediction
interval into a few readable calls.

## Install

```r
# install.packages("remotes")
remotes::install_github("MJ-Ershad/metahelpers")
```

`metafor` and `clubSandwich` are Suggests — install them to use the modeling
functions.

## Usage

```r
library(metahelpers)

fit <- fit_three_level(
  data    = dat,
  yi      = yi,        # effect sizes
  vi      = vi,        # sampling variances
  cluster = study,     # higher level (study / cohort)
  robust  = TRUE       # CR2 cluster-robust SEs
)

summary(fit)
#> == metahelpers three-level fit ==
#> Three-level random-effects meta-analysis (metahelpers)
#>   60 effect sizes in 30 clusters
#>   pooled estimate = 0.41  (CR2 SE = 0.05, p = 1e-09)
#>   ...
#>   Multilevel I^2 (total) = 71.3%
#>   95% prediction interval: [0.02, 0.80]
```

## What you get

| function               | returns                                                   |
|------------------------|-----------------------------------------------------------|
| `fit_three_level()`    | the `rma.mv` model + CR2 `coef_test`, as one object       |
| `multilevel_i2()`      | total `I^2` split into between- and within-cluster shares |
| `prediction_interval()`| interval for the true effect of a new comparable study    |
| `tidy_rema()`          | a one-row data frame for stacking many models             |

`tidy_rema()` makes it easy to fit a model per subgroup and bind the rows:

```r
rows <- lapply(split(dat, dat$subgroup), function(d)
  tidy_rema(fit_three_level(d, yi, vi, study)))
do.call(rbind, rows)
```

## Tests

```r
devtools::test()   # synthetic 30-study dataset; mean recovery + tidy/I^2 checks
```

## References

- Cheung, M. W.-L. (2014). Modeling dependent effect sizes with three-level
  meta-analyses. *Psychological Methods*, 19(2), 211-229.
- Pustejovsky & Tipton (2018). Small-sample cluster-robust variance estimators.

## License

MIT (c) 2026 Mohamadjavad Ershadmanesh
