# `lavaan.mi`

This is an R package whose primary purpose is to extend the functionality of the R package `lavaan`. When incomplete data have been multiply imputed, the `lavaan.mi` package allows users to fit structural equation models as they would using the `lavaan` package. 

The development version can be installed from GitHub using the [`remotes` package](https://cran.r-project.org/package=remotes):

```
remotes::install_github("TDJorgensen/lavaan.mi")
```


## Connecting `lavaan.mi` to `lavaan`

The top 4 lines of the table below indicate the `lavaan.mi::lavaan.mi()` function corresponds to the `lavaan::lavaan()` function, as do the wrappers `cfa.mi()`, `sem.mi()`, and `growth.mi()`, all of which return an object of class `lavaan.mi`.


|   `lavaan` Function      |   `lavaan.mi`    Function   |
|:-------------------------|:----------------------------|
|      `lavaan()`          |       `lavaan.mi()`         |
|        `cfa()`           |         `cfa.mi()`          |
|        `sem()`           |         `sem.mi()`          |
|      `growth()`          |       `growth.mi()`         |
|  `parameterEstimates()`  |  `parameterEstimates.mi()`  |
| `standardizedSolution()` | `standardizedSolution.mi()` |
|     `lavTestLRT()`       |      `lavTestLRT.mi()`      |
|     `lavTestWald()`      |      `lavTestWald.mi()`     |
|     `lavTestScore()`     |      `lavTestScore.mi()`    |
|     `modIndices()`       |      `modIndices.mi()`      |

The `data=` argument must be a `list` of `data.frame`s rather than a single `data.frame`, and the specified `lavaan::model.syntax` will be applied to each `data.frame` in the `list`.


## Inheritance from `lavaanList`

The `lavaan.mi` class inherits from class `lavaanList` (from the `lavaan` package), extending it with a few additional and customized slots (see `class?lavaanList` and `class?lavaan.mi` for details).

It is not wise to use "lazy loading" (i.e., calling `lavaan.mi::lavaan.mi()` without loading `library(lavaan.mi)` explicitly) because R may have trouble finding the correct methods in the package for applying generic functions (e.g., `summary`) to `lavaan.mi-class` objects.  This occasionally happened in the `semTools` package as well. When the correct methods are not found, then the generic function will revert to the methods written for `lavaanList-class` objects (e.g., `summary()` will only have an `est.ave` column for the average estimate across data sets).  In most cases (e.g., `fitted()`, `anova()`, `residuals()`, etc.), this will cause an error because few methods have been written for `lavaanList-class` objects.

Additional `lavaan` output can be obtained per imputation using the `FUN=` argument, which is applied to a `lavaan`-class object formed from the result of each imputation (prior to its being stored more succinctly in a `lavaanList`-class object).  The last `?lavaan.mi` help-page example demonstrates how to obtain and extract additional output.


## Methods for `lavaan` and `lavaan.mi` objects

Many remaining methods written for `lavaan`-class objects are also available for `lavaan.mi`-class objects (also see `class?lavaan.mi` for a list and descriptions of arguments):

- `show()`
- `summary()`
- `fitMeasures()`
- `anova()`
- `nobs()`
- `coef()`
- `vcov()`
- `fitted.values()` and alias `fitted()`
- `residuals()` and alias `resid()`


## Connecting `lavaan.mi` to `semTools`

Many of the functions in the `semTools` package continue to be available for `lavaan.mi`-class objects (e.g., `semTools::reliability()`), and functions such as `semTools::parcelAllocation()` and `semTools::plausibleValues()` are enhanced by allowing their output to be analyzed using `lavaan.mi()`.


