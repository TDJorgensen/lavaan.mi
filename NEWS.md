# lavaan.mi 0.1-0 (in development)

The `semTools` package has a **Missing Data** suite, which for years included utilities for fitting structural equation models to multiply imputed data sets with the `lavaan` package and automatically pooling the results. These multiple-imputation features continue to grow, justifying their own separate package (and current deprecation in `semTools`).

## Differences from `semTools`:

- The `runMI()` function does not exist.  The `semTools::runMI()` function had a `fun=` argument to specify which `lavaan` function to call (`lavaan()` or the `cfa()` and `sem()` wrappers).  Instead, analogous to `lavaan`, the basic function will be `lavaan.mi()`, with `cfa.mi()` and `sem.mi()` as wrappers around `lavaan.mi()`.
- Unlike `semTools::runMI()`, users may **not** impute data indirectly via `miPackage=` and `miArgs=` arguments.  This feature was originally included to mimic M*plus*, which allows for model-based imputation.  Because `semTools::runMI()` did not use the specified SEM as an imputation model, the mimicry of the M*plus* functionality was only cosmetic.  Even in `semTools::runMI()`, it was always better to impute (with full control) using dedicated imputation software first (e.g., the `Amelia` or `mice` package, or external software like [Blimp](https://www.appliedmissingdata.com/blimp)), to be analyzed with `lavaan.mi()`.
- Because users can pass a vector of multiple (e.g., robust) test statistics to the `lavaan(test=)` argument, there is a new `test=` argument for `lavaan::lavTestLRT()`, enabling users to select which specific $\Delta \chi^2$ test should be calculated for a model comparison.  Accordingly, the `lavTestLRT.mi()` function's `test=` argument has been renamed `pool.method=`, so that users can still pass a `test=` argument to `lavTestLRT.mi(...)`.


## New Features:

- The [newly proposed method](https://doi.org/10.5705/ss.202019.0314) for pooling LRT statistics is the new default (`pool.method = "D4"`) for `lavTestLRT.mi()`, replacing the old default (but `pool.method = "D3"` is still available).  The D4 pooling method performs similarly to D3 (see Grund et al.'s [simulation study](https://doi.org/10.31234/osf.io/d459g)) but is less computationally intensive.
- Because users can pass a vector of multiple (e.g., robust) test statistics to the `lavaan(test=)` argument, the user can choose which statistic can be pooled using `lavTestLRT.mi(..., pool.method = "D2", pool.robust = TRUE)`.  This is controlled via the new `scaled.test=` argument.  As an alternative to the standard $\chi^2$ statistic, the user can also pool [Browne's (1984)](https://doi.org/10.1111/j.2044-8317.1984.tb00789.x) residual-based statistic using, e.g., `lavTestLRT.mi(..., pool.method = "D2", standard.test = "browne.residual.nt")` with the default `pool.robust=FALSE`.
- New functions have been written to provide the same `lavaan` functionality for pooled results: `paramterEstimates.mi()`, `standardizedSolution.mi()`, and `lavResiduals.mi()`. 
- The `poolSat()` function implements Cai and colleagues' ([2012](https://doi.org/10.3102/1076998612458320), [2019](https://doi.org/10.1080/00273171.2018.1523000)) method of pooling the saturated-model results, then fitting SEM(s) to those pooled summary statistics.  This can be much more computationally efficient than fitting each SEM of interest to each imputation and pooling the results of each estimated SEM.  Instead, pooling only occurs for the summary statistics, and each SEM is estimated as it would be for complete data.


## Bug Fixes:

- When the `VGAM` package was loaded, it created a [conflict](https://github.com/simsem/semTools/issues/89) with finding the `anova() ` method for `lavaan.mi` objects in `semTools`. This was resolved in February 2022 when `VGAM` version 1.1-6 was sent to CRAN.

## Known Bugs:

- `fitMeasures()` returns an error for multigroup MLSEMs, unless no (S/C)RMR index is requested in the `fit.measures=` argument.


Track our progress or report issues on GitHub:

https://github.com/TDJorgensen/lavaan.mi
