### Terrence D. Jorgensen & Yves Rosseel
### Last updated: 30 May 2024
### Pooled likelihood ratio test for multiple imputations
### Borrowed source code from lavaan/R/lav_test_LRT.R


## -------------
## Main function
## -------------

##' Likelihood Ratio Test for Multiple Imputations
##'
##' Likelihood ratio test (LRT) for lavaan models fitted to multiple imputed
##' data sets.
##'
##' The `"D2"` method is available using any estimator and test statistic.
##' When using a likelihood-based estimator, 2 additional methods are available
##' to pool the LRT.
##' \itemize{
##'   \item The Meng & Rubin (1992) method, commonly referred to as `"D3"`.
##'         This method has many problems, discussed in Chan & Meng (2022).
##'   \item The Chan & Meng (2022) method, referred to as `"D4"` by
##'         Grund et al. (2023), resolves problems with `"D3"`.
##' }
##' When `"D2"` is not explicitly requested in situations it is the only
##' applicable method, (e.g., DWLS for categorical outcomes), users are notified
##' that `pool.method` was set to `"D2"`.
##'
##' `pool.method = "Mplus"` implies `"D3"` and `asymptotic = TRUE`
##' (see Asparouhov & Muthen, 2010).
##'
##' Note that the `anova()` method simply calls `lavTestLRT.mi()`.
##'
##' @aliases lavTestLRT.mi
##' @importFrom lavaan lavTestLRT
##'
##' @param object An object of class [lavaan.mi-class]
##' @param ... Additional objects of class [lavaan.mi-class], as
##'   well as arguments passed to [lavaan::lavTestLRT()] when
##'   `pool.method = "D2"` and `pool.robust = TRUE`.
##' @param modnames Optional `character` of model names to use as row names
##'   in the resulting matrix of results (when more than 2 models are compared)
##' @param asANOVA `logical` indicating whether to return an object of
##'   class `"anova"`.  If `FALSE`, a numeric vector is returned for
##'   one (pair of) model(s), or a `data.frame` is returned for multiple
##'   pairs of models.
##' @param pool.method `character` indicating which pooling method to use.
##'   \itemize{
##'     \item `"D4"`, `"new.LRT"`, `"cm"`, or `"chan.meng"`
##'           requests the method described by Chan & Meng (2022).
##'           This is currently the default.
##'     \item `"D3"`, `"old.LRT"`, `"mr"`, or `"meng.rubin"`
##'           requests the method described by Meng & Rubin (1992).
##'     \item `"D2"`, `"LMRR"`, or `"Li.et.al"` requests the
##'           complete-data LRT statistic should be calculated using each
##'           imputed data set, which will then be pooled across imputations, as
##'           described in Li, Meng, Raghunathan, & Rubin (1991).
##'   }
##'   Find additional details in Enders (2010, chapter 8).

#FIXME: Remove these 2 (passed via ...)
## @param standard.test
## @param scaled.test

##' @param omit.imps `character` vector specifying criteria for omitting
##'   imputations from pooled results.  Can include any of
##'   `c("no.conv", "no.se", "no.npd")`, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (`"no.npd"`) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases. Specific imputation numbers can also be included in this
##'   argument, in case users want to  apply their own custom omission criteria
##'   (or simulations can use different numbers of imputations without
##'   redundantly refitting the model).
##' @param asymptotic `logical`. If `FALSE` (default), the pooled test
##'   will be returned as an *F*-distributed statistic with numerator
##'   (`df1`) and denominator (`df2`) degrees of freedom.
##'   If `TRUE`, the pooled *F* statistic will be multiplied by its
##'   `df1` on the assumption that its `df2` is sufficiently large
##'   enough that the statistic will be asymptotically \eqn{\chi^2} distributed
##'   with `df1`.
##' @param pool.robust `logical`. Ignored unless `pool.method = "D2"` and a
##'   robust test was requested. If `pool.robust = TRUE`, the robust test
##'   statistic is pooled, whereas `pool.robust = FALSE` will pool
##'   the naive test statistic (or difference statistic) and apply the average
##'   scale/shift parameter to it. The harmonic mean is applied to the scaling
##'   factor, whereas the arithmetic mean is applied to the shift parameter.
##'
##' @return
##'
##'   - When `asANOVA=TRUE`, returns an object of class [stats::anova] with a
##'     a test of model fit for a single model (`object`) or test(s) of the
##'     difference(s) in fit between nested models passed via `...` (either an
##'     `F` or \eqn{\chi^2} statistic, depending on the `asymptotic` argument),
##'     its degrees of freedom, its *p* value, and 2 missing-data diagnostics:
##'     the relative increase in variance (RIV = FMI / (1 \eqn{-} FMI)) and the
##'     fraction of missing information (FMI = RIV / (1 + RIV)).
##'   - When `asANOVA=FALSE`, returns a vector containing the LRT statistic for
##'     a single model or comparison of a single pair of models, or a
##'     `data.frame` of multiple model comparisons. Robust statistics will also
##'     include the average (across imputations) scaling factor and
##'     (if relevant) shift parameter(s), unless `pool.robust = TRUE`.
##'     When using `pool.method = "D3"` or `"D4"`, the vector for a single
##'     model also includes its average log-likelihood and information criteria.
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##'   Based on source code for [lavaan::lavTestLRT()] by Yves Rosseel.
##'
##' @references
##'   Asparouhov, T., & Muthen, B. (2010). *Chi-square statistics
##'   with multiple imputation*. Technical Report. Retrieved from
##'   <http://www.statmodel.com/>
##'
##'   Chan, K. W., & Meng, X. L. (2022). Multiple improvements of multiple
##'   imputation likelihood ratio tests. *Statistica Sinica, 32*,
##'   1489--1514. \doi{10.5705/ss.202019.0314}
##'
##'   Enders, C. K. (2010). *Applied missing data analysis*.
##'   New York, NY: Guilford.
##'
##'   Grund, S., Lüdtke, O., & Robitzsch, A. (2023). Pooling methods for
##'   likelihood ratio tests in multiply imputed data sets.
##'   *Psychological Methods, 28*(5), 1207--1221.
##'   \doi{10.1037/met0000556}
##'
##'   Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
##'   Significance levels from repeated *p*-values with multiply-imputed
##'   data. *Statistica Sinica, 1*(1), 65--92. Retrieved from
##'   <https://www.jstor.org/stable/24303994>
##'
##'   Meng, X.-L., & Rubin, D. B. (1992). Performing likelihood ratio tests with
##'   multiply-imputed data sets. *Biometrika, 79*(1), 103--111.
##'   \doi{10.2307/2337151}
##'
##'   Rubin, D. B. (1987). *Multiple imputation for nonresponse in surveys*.
##'   New York, NY: Wiley. \doi{10.1002/9780470316696}
##'
##' @seealso [lavaan::lavTestLRT()] for arguments that can be passed via \dots,
##' and use [lavaan::fitMeasures()] to obtain fit indices calculated from pooled
##' test statistics.
##'
##' @examples
##' data(HS20imps) # import a list of 20 imputed data sets
##'
##' ## specify CFA model from ?lavaan::cfa help page
##' HS.model <- '
##'   visual  =~ x1 + x2 + x3
##'   textual =~ x4 + x5 + x6
##'   speed   =~ x7 + x8 + x9
##' '
##'
##' fit1 <- cfa.mi(HS.model, data = HS20imps, estimator = "mlm")
##'
##' ## more constrained model: parallel indicators
##' HS.parallel <- '
##'   visual  =~ x1 + 1*x2 + 1*x3
##'   textual =~ x4 + 1*x5 + 1*x6
##'   speed   =~ x7 + 1*x8 + 1*x9
##' '
##'
##' fitp <- cfa.mi(HS.parallel, data = HS20imps, estimator = "mlm")
##'
##' ## Even more constrained model: orthogonal factors
##' fit0 <- cfa.mi(HS.parallel, data = HS20imps, estimator = "mlm",
##'                orthogonal = TRUE)
##'
##' ## By default, pool.method = "D4".
##' ## Must request an asymptotic chi-squared statistic
##' ## in order to accommodate a robust correction.
##' lavTestLRT.mi(fit1, fit0, fitp, asymptotic = TRUE)
##' ## or   anova(fit1, fit0, fitp, asymptotic = TRUE)
##'
##' ## Pass any argument to lavTestLRT()
##' lavTestLRT.mi(fit1, fit0, fitp, asymptotic = TRUE,
##'               method = "satorra.bentler.2010")
##'
##' ## For a single model, you can request a vector instead of an anova-class
##' ## table in order to see naive information criteria (only using D3 or D4),
##' ## which are calculated using the average log-likelihood across imputations.
##' lavTestLRT.mi(fit1, asANOVA = FALSE)
##'
##'
##'
##' ## When using a least-squares (rather than maximum-likelihood) estimator,
##' ## only the D2 method is available.  For example, ordered-categorical data:
##' data(binHS5imps) # import a list of 5 imputed data sets
##'
##' ## fit model using default DWLS estimation
##' fit1c <- cfa.mi(HS.model   , data = binHS5imps, ordered = TRUE)
##' fit0c <- cfa.mi(HS.parallel, data = binHS5imps, ordered = TRUE,
##'                 orthogonal = TRUE)
##'
##' ## Using D2, you can either robustify the pooled naive statistic ...
##' lavTestLRT.mi(fit1c, fit0c, asymptotic = TRUE, pool.method = "D2")
##' ## ... or pool the robust chi-squared statistic (NOT recommended)
##' lavTestLRT.mi(fit1c, fit0c, asymptotic = TRUE, pool.method = "D2",
##'               pool.robust = TRUE)
##'
##' ## When calculating fit indices, you can pass lavTestLRT.mi() arguments:
##' fitMeasures(fit1c, output = "text",
##'             # lavTestLRT.mi() arguments:
##'             pool.method = "D2", pool.robust = TRUE)
##'
##' @export
lavTestLRT.mi <- function(object, ..., modnames = NULL, asANOVA = TRUE,
                          pool.method = c("D4","D3","D2"),
                          omit.imps = c("no.conv","no.se"),
                          asymptotic = FALSE, pool.robust = FALSE) {
  ## save model names
  objname <- deparse(substitute(object))
  dotnames <- as.character(sapply(substitute(list(...))[-1], deparse))

  ## check class
  if (!inherits(object, "lavaan.mi")) stop("object is not class 'lavaan.mi'")
  ## check for additional arguments
  dots <- list(...)
  if (length(dots)) {
    ## separate lavaan.mi objects from lavTestLRT() arguments
    idx.mi <- which(sapply(dots, inherits, what = "lavaan.mi"))
    if (length(idx.mi)) {
      mods <- dots[idx.mi]
      dots <- dots[-idx.mi]
      if (!is.null(modnames)) {
        objname  <- modnames[1] # remove name for "object"
        modnames <- modnames[-1]
        stopifnot(length(modnames) == length(mods))
        names(mods) <- modnames
      } else {
        modnames <- dotnames[idx.mi]
        nonames <- which(names(mods) == "")
        names(mods)[nonames] <- modnames[nonames]
      }

    } else {
      mods <- NULL
      modnames <- NULL
    }

    ## any tests named for D2.LRT()?
    if ("standard.test" %in% names(dots)) {
      standard.test <- dots$standard.test
    } else standard.test <- NULL
    if ("scaled.test" %in% names(dots)) {
      scaled.test <- dots$scaled.test
    } else scaled.test <- NULL

    ## otherwise, only save lavTestLRT() arguments
    LRT.names <- intersect( names(dots) , names(formals(lavTestLRT)) )
    dots <- if (length(LRT.names)) dots[LRT.names] else NULL

  } else {
    mods          <- NULL
    standard.test <- NULL
    scaled.test   <- NULL
  }

  ## compare 1 or multiple pairs?
  if (length(mods) < 2L) {
    argList <- c(list(object = object), dots,
                 list(pool.method = pool.method, omit.imps = omit.imps,
                      asymptotic = asymptotic, pool.robust = pool.robust))
    if (length(mods) == 1L)      argList$h1            <- mods[[1]]
    if (!is.null(standard.test)) argList$standard.test <- standard.test
    if (!is.null(  scaled.test)) argList$scaled.test   <- scaled.test
    results <- do.call(pairwiseLRT, argList)

    ## store in an anova-class table?
    if (asANOVA) {

      if (length(mods) == 1L) {
        argList$asymptotic <- TRUE # table includes pooled naive chi-squared per model
        argList$h1 <- NULL
        out0 <- do.call(pairwiseLRT, argList)
        argList$object <- mods[[1]]
        out1 <- do.call(pairwiseLRT, argList)

        out <- data.frame(Df    = sort(c(out0["df"]   , out1["df"]   ) ),
                          Chisq = sort(c(out0["chisq"], out1["chisq"]) ))
        rownames(out) <- c(objname, modnames)[order(c(out0["df"], out1["df"]))]

      } else {
        ## results already contain single-model test, but might be F test
        if ("chisq" %in% names(results)) {
          resChi <- results
        } else {
          argList$asymptotic <- TRUE # table includes pooled naive chi-squared per model
          resChi <- do.call(pairwiseLRT, argList)
        }
        out <- data.frame(Df    = c(0, resChi["df"]),
                          Chisq = c(0, resChi["chisq"]))
        rownames(out) <- c("Saturated", objname)
      }
      class(out) <- c("anova", "data.frame")

      if ("chisq" %in% names(results)) {
        robust <- "chisq.scaled" %in% names(results)
        out$`Chisq diff` <- c(NA, ifelse(robust, results["chisq.scaled"], results["chisq"]))
        out$`Df diff`    <- c(NA, ifelse(robust, results["df.scaled"], results["df"]))
        out$`Pr(>Chisq)` <- c(NA, ifelse(robust, results["pvalue.scaled"], results["pvalue"]))
      } else {
        robust <- "F.scaled" %in% names(results)
        out$`F`      <- c(NA, ifelse(robust, results["F.scaled"]  , results["F"]))
        out$`num Df` <- c(NA, ifelse(robust, results["df1.scaled"], results["df1"]))
        out$`den Df` <- c(NA, ifelse(robust, results["df2.scaled"], results["df2"]))
        out$`Pr(>F)` <- c(NA, ifelse(robust, results["pvalue.scaled"], results["pvalue"]))
      }
      out$RIV <- c(NA, results["ariv"])
      out$FMI <- c(NA, results["fmi" ])

    } else return(results)

  } else {
    modList <- c(list(object), mods)
    names(modList) <- c(objname, modnames)
    argList <- c(modList, argsLRT = list(dots),
                 list(pool.method = pool.method, omit.imps = omit.imps,
                      asymptotic = asymptotic, pool.robust = pool.robust))
    if (!is.null(standard.test)) argList$standard.test <- standard.test
    if (!is.null(  scaled.test)) argList$scaled.test   <- scaled.test
    results <- do.call(multipleLRT, argList)

    if (asANOVA) {
      ## get each model's pooled naive chi-squared
      outList <- lapply(modList, pairwiseLRT, asymptotic = TRUE,
                        pool.method = pool.method, pool.robust = pool.robust,
                        omit.imps = omit.imps)
      ## extract pooled test and degrees of freedom
      x2List <- sapply(outList, "[", i = "chisq")
      dfList <- sapply(outList, "[", i = "df")
      ## sort in order of less-to-more constrained
      out <- data.frame(Df = sort(dfList),
                        Chisq = x2List[order(dfList)])
      rownames(out) <- names(modList)[order(dfList)]
      ## add an empty top row to model-comparison results
      resPlus1 <- rbind(NA, results)
      rownames(resPlus1) <- rownames(out)
      if (grepl("chisq", colnames(resPlus1)[1])) {
        colnames(resPlus1) <- c("Chisq diff", "Df diff", "Pr(>Chisq)", "RIV", "FMI")
      } else {
        colnames(resPlus1) <- c("F", "num Df", "den Df", "Pr(>F)", "RIV", "FMI")
      }
      ## combine and add class
      out <- cbind(out, resPlus1)
      class(out) <- c("anova", "data.frame")
    } else return(results)

  }

  out
}


## ---------------------------
## anova() method is a wrapper
## ---------------------------

##' @name lavTestLRT.mi
##' @aliases anova,lavaan.mi-method
##' @export
setMethod("anova", "lavaan.mi",
          ## TDJ borrowed from lavaan's anova method:
function(object, ...) {
  # NOTE: if we add additional arguments, it is not the same generic
  # anova() function anymore, and match.call will be screwed up

  # NOTE: we need to extract the names of the models from match.call here,
  #       otherwise, we loose them in the call stack

  mcall <- match.call(expand.dots = TRUE)
  dots <- list(...)

  # catch SB.classic and SB.H0
  #SB.classic <- TRUE; SB.H0 <- FALSE

  #arg.names <- names(dots)
  #arg.idx <- which(nchar(arg.names) > 0L)
  #if(length(arg.idx) > 0L) {
  #    if(!is.null(dots$SB.classic))
  #        SB.classic <- dots$SB.classic
  #    if(!is.null(dots$SB.H0))
  #        SB.H0 <- dots$SB.H0
  #    dots <- dots[-arg.idx]
  #}

  modp <- if (length(dots)) sapply(dots, inherits, "lavaan.mi") else logical(0)
  mods <- c(list(object), dots[modp])
  NAMES <- sapply(as.list(mcall)[c(FALSE, TRUE, modp)], deparse)

  # use do.call to handle changed dots
  #ans <- do.call("lavTestLRT", c(list(object = object,
  #               SB.classic = SB.classic, SB.H0 = SB.H0,
  #               model.names = NAMES), dots))

  #ans
  lavTestLRT.mi(object = object, ..., modnames = NAMES)
})



#FIXME: delete anova_lavaan_mi (not needed)
##' @importFrom stats anova
##' @importFrom lavaan lavListInspect lavTestLRT
anova_lavaan_mi <- function(object, ...) {
  ## save model names
  objname <- deparse(substitute(object))
  dotnames <- as.character(sapply(substitute(list(...))[-1], deparse))

  ## check class
  if (!inherits(object, "lavaan.mi")) stop("object is not class 'lavaan.mi'")
  ## check for additional arguments
  dots <- list(...)
  if (length(dots)) {
    ## separate lavaan.mi objects from other lavTestLRT.mi() arguments
    idx.mi <- which(sapply(dots, inherits, what = "lavaan.mi"))
    if (length(idx.mi)) {
      mods <- dots[idx.mi]
      dots <- dots[-idx.mi]
      ## save names for mods, so compareFit() doesn't break
      modnames <- dotnames[idx.mi]
      nonames <- which(names(mods) == "")
      names(mods)[nonames] <- modnames[nonames]
    } else {
      mods <- NULL
      modnames <- NULL
    }
    LRT.names <- intersect(names(dots),
                           union(names(formals(lavTestLRT)),
                                 names(formals(lavTestLRT.mi))))
    dots <- if (length(LRT.names)) dots[LRT.names] else NULL
    if (!is.null(dots$h1)) {
      #FIXME: this shouldn't be necessary: mods <- c(mods, list(h1 = dots$h1))
      dots$h1 <- NULL
    }
  } else mods <- NULL

  ## run semTools::compareFit if length(idx.mi) > 1L
  if (length(mods) == 0L) {
    argList <- c(list(object = object), dots)
    results <- do.call(lavTestLRT.mi, argList)
  } else if (length(mods) == 1L) {
    argList <- c(list(object = object, h1 = mods[[1]]), dots)
    results <- do.call(lavTestLRT.mi, argList)
  } else if (length(mods) > 1L) {
    ## is semTools::compareFit available?
    if (!requireNamespace("semTools", quietly = TRUE)) {
      stop('The semTools package is required to simultaneously compare ',
           'multiple lavaan.mi models using a the anova() method. Without ',
           'semTools, a single pair of models can be compared using ',
           'lavTestLRT.mi()')
    }
    if (!"package:semTools" %in% search()) attachNamespace("semTools")
    modList <- c(list(object), mods)
    names(modList) <- c(objname, modnames)
    argList <- c(modList, list(argsLRT = dots, indices = FALSE))
    results <- do.call(semTools::compareFit, argList)@nested
    class(results) <- c("lavaan.data.frame","data.frame")
    attr(results, "header") <- "Nested Model Comparisons:"
  }

  results
}



## ----------------
## Hidden Functions
## ----------------


##' @importFrom lavaan lavListInspect parTable lavTestLRT
D2.LRT <- function(object, h1 = NULL, useImps, asymptotic = FALSE,
                   harm = TRUE, # harmonic or arithmetic mean of scaling.factor
                   return.mean.chisq  = FALSE, # TRUE to ease D4
                   return.scale.shift = FALSE, # TRUE for shift in robustify()
                   standard.test = "standard", scaled.test = "default",
                   pool.robust = FALSE, LRTargs = list()) {
  warn <- lavListInspect(object, "options")$warn

  if (pool.robust && !is.null(h1)) {
    PT1 <- parTable(h1)
    op1 <- lavListInspect(h1, "options")
    oldCall <- object@lavListCall #re-run lavaanList() and save DIFFTEST
    if (!is.null(oldCall$parallel)) {
      if (oldCall$parallel == "snow") {
        oldCall$parallel <- "no"
        oldCall$ncpus <- 1L
        if (warn) warning("Unable to pass lavaan::lavTestLRT() arguments when ",
                          "parallel = 'snow'. Switching to parallel = 'no'. ",
                          "Unless using Windows, parallel = 'multicore' works.")
      }
    }

    ## call lavaanList() again to run lavTestLRT() on each imputation
    oldCall$FUN <- function(obj) {
      fit1 <- try(lavaan::lavaan(PT1, slotOptions = op1, slotData = obj@Data),
                  silent = TRUE)
      if (inherits(fit1, "try-error")) {
        return("fit failed")
      } else {
        argList <- c(list(object = obj, fit1), LRTargs)
      }
      out <- try(do.call(lavTestLRT, argList),
                 silent = TRUE)
      if (inherits(out, "try-error")) return("lavTestLRT() failed")
      c(chisq = out[2, "Chisq diff"], df = out[2, "Df diff"],
        scale = attr(out, "scale")[2], shift = attr(out, "shift")[2])
    }
    FIT <- eval(as.call(oldCall))
    ## check if there are any results
    noFit <- sapply(FIT@funList, function(x) x[1] == "fit failed")
    noLRT <- sapply(FIT@funList, function(x) x[1] == "lavTestLRT() failed")
    if (all(noFit | noLRT)) stop("No success using lavTestScore() on any imputations.")

    if (return.scale.shift) {
      SCALE <- sapply(FIT@funList[ intersect(which(!(noFit | noLRT)), useImps) ],
                      "[[", i = "scale")
      SHIFT <- sapply(FIT@funList[ intersect(which(!(noFit | noLRT)), useImps) ],
                      "[[", i = "shift")
      return(data.frame(scale = SCALE, shift = SHIFT))
    }

    chiList <- sapply(FIT@funList[ intersect(which(!(noFit | noLRT)), useImps) ],
                      "[[", i = "chisq")
    if (return.mean.chisq) return(mean(chiList, na.rm = TRUE))

    dfList <- sapply(FIT@funList[ intersect(which(!(noFit | noLRT)), useImps) ],
                     "[[", i = "df")

    out <- calculate.D2(chiList, DF = mean(dfList), asymptotic)
    names(out) <- paste0(names(out), ".scaled")
    class(out) <- c("lavaan.vector","numeric")
    return(out)
  }
  ## else, return model fit OR naive difference test to be robustified


  test.names <- lavListInspect(object, "options")$test
  if (standard.test[1] %in% c("standard", "browne.residual.nt",
                              "browne.residual.nt.model",
                              "browne.residual.adf",
                              "browne.residual.adf.model")) {
    standard.test <- standard.test[1]
  } else {
    standard.test <- intersect(test.names, c("standard", "browne.residual.nt",
                                             "browne.residual.nt.model",
                                             "browne.residual.adf",
                                             "browne.residual.adf.model"))[1]
  }

  # any non-standard test stats?
  if (length(test.names) > 1L) {
    ## remove standard and any bootstrapped tests
    rm.idx <- which(test.names %in% c("none","default","standard",
                                      "browne.residual.nt",
                                      "browne.residual.nt.model",
                                      "browne.residual.adf",
                                      "browne.residual.adf.model"))
    if (length(rm.idx) > 0L) test.names <- test.names[-rm.idx]
    ## only acknowledge one scaled test statistic
    if (length(test.names) > 1L) {
      if (scaled.test[1] %in% c("satorra.bentler",
                                "yuan.bentler", "yuan.bentler.mplus",
                                "mean.var.adjusted", "scaled.shifted")) {
        scaled.test <- scaled.test[1]
      } else {
        scaled.test <- intersect(test.names,
                                 c("satorra.bentler",
                                   "yuan.bentler", "yuan.bentler.mplus",
                                   "mean.var.adjusted", "scaled.shifted"))[1]
      }

    }
  }

  ## pool Wald tests
  if (is.null(h1)) {
    test <- if (pool.robust) scaled.test else standard.test
    DF <- mean(sapply(object@testList[useImps], function(x) x[[test]]$df  ))
    w  <-      sapply(object@testList[useImps], function(x) x[[test]]$stat)
  } else {
    ## this will not get run if pool.robust because logic catches that first
    DF0 <- mean(sapply(object@testList[useImps], function(x) x[[standard.test]]$df))
    DF1 <- mean(sapply(    h1@testList[useImps], function(x) x[[standard.test]]$df))
    DF <- DF0 - DF1
    w0 <- sapply(object@testList[useImps], function(x) x[[standard.test]]$stat)
    w1 <- sapply(    h1@testList[useImps], function(x) x[[standard.test]]$stat)
    w <- w0 - w1
    ## check DF
    if (DF < 0) {
      w <- -1*w
      DF <- -1*DF
    }
  }
  if (return.mean.chisq) return(mean(w, na.rm = TRUE))

  out <- calculate.D2(w, DF, asymptotic)
  ## add .scaled suffix
  if (pool.robust) names(out) <- paste0(names(out), ".scaled")
  ## for 1 model, add extra info
  if (is.null(h1)) {
    if (pool.robust) {
      ## for consistent fitMeasures() output, add average scaling factor
      scaleList <- sapply(object@testList[useImps],
                          function(x) x[[test]][["scaling.factor"]] )
      if (test == "scaled.shifted") {
        shiftList <- sapply(object@testList[useImps],
                            function(x) x[[test]][["shift.parameter"]])
      }

      if (harm) {
        sc <- c(chisq.scaling.factor = 1/mean(1/scaleList))
        if (test == "scaled.shifted") {
          sc["chisq.shift.parameter"] <- 1/mean(1/shiftList)
        }

      } else {
        sc <- c(chisq.scaling.factor = mean(scaleList))
        if (test == "scaled.shifted") {
          sc["chisq.shift.parameter"] <- mean(shiftList)
        }
      }
      out <- append(out, values = sc,
                    after = which(names(out) == "pvalue.scaled"))
    } else {
      ## redundant if pool.robust, because concatenated with pooled naive stat
      PT <- parTable(object)
      out <- c(out, npar = max(PT$free) - sum(PT$op == "=="),
               ntotal = lavListInspect(object, "ntotal"))
    }
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}

##' @importFrom lavaan parTable lavaan lavListInspect
getLLs <- function(object, useImps, saturated = FALSE, # default h1 model
                   omit.imps = c("no.conv","no.se")) {
  ## FIXME: lavaanList does not return info when fixed because no convergence!
  dataList <- object@DataList[useImps]
  lavoptions <- lavListInspect(object, "options")

  group <- lavListInspect(object, "group")
  if (length(group) == 0L) group <- NULL
  cluster <- lavListInspect(object, "cluster")
  if (length(cluster) == 0L) cluster <- NULL

  if (saturated) {
    #FIXME: below is legacy code, no longer needed?
    # fit <- lavaan(parTable(object), data = dataList[[ useImps[1] ]],
    #               slotOptions = lavoptions, group = group, cluster = cluster)
    # ## use saturated parameter table as new model
    # PT <- lavaan::lav_partable_unrestricted(fit)
    # ## fit saturated parameter table to each imputation, return estimates
    # satParams <- lapply(object@DataList[useImps], function(d) {
    #   parTable(lavaan(model = PT, data = d, slotOptions = lavoptions,
    #                   group = group, cluster = cluster))$est
    # })
    # ## fix them to pooled estimates
    # PT$ustart <- colMeans(do.call(rbind, satParams))

    PT <- object@h1List[[ useImps[1] ]]$PT
    coefList <- lapply(object@h1List[useImps], function(x) x$PT$ustart)
    PT$ustart <- colMeans(do.call(rbind, coefList))

    ## set all parameters fixed
    PT$free <- 0L
    PT$user <- 1L
    PT$start <- NULL
    PT$est <- NULL
    PT$se <- NULL
  } else {
    ## save parameter table as new model
    PT <- parTable(object)
    ## set all parameters fixed
    PT$free <- 0L
    PT$user <- 1L
    ## fix them to pooled estimates
    fixedValues <- coef_lavaan_mi(object, type = "user", omit.imps = omit.imps)
    PT$ustart <- fixedValues
    PT$start <- NULL
    PT$est <- NULL
    PT$se <- NULL
    ## omit (in)equality constraints and user-defined parameters
    params <- !(PT$op %in% c("==","<",">",":="))
    PT <- PT[params, ]
  }
  ## return log-likelihoods
  sapply(dataList, function(d) {
    lavaan::logLik(lavaan(PT, data = d, slotOptions = lavoptions,
                          group = group, cluster = cluster))
  })
  #TODO: use "dry-run" trick from blavaan:::postpred() to save computing time
}

##' @importFrom stats pf pchisq
##' @importFrom lavaan lavListInspect parTable
D3.LRT <- function(object, h1 = NULL, useImps, asymptotic = FALSE,
                   omit.imps = c("no.conv","no.se")) {
  ## NOTE: Need to pass omit.imps= to getLLs(), which calls coef() method

  N <- lavListInspect(object, "ntotal")
  m <- length(useImps)

  if (is.null(h1)) {
    DF <- object@testList[[ useImps[1] ]]$standard[["df"]]
  } else {
    DF1 <- h1@testList[[ useImps[1] ]]$standard[["df"]]
    DF0 <- object@testList[[ useImps[1] ]]$standard[["df"]]
    DF <- DF0 - DF1
    if (DF < 0) stop('The "object" model must be nested within (i.e., have ',
                     'fewer degrees of freedom than) the "h1" model.')
  }

  ## calculate m log-likelihoods under pooled H0 estimates
  LL0 <- getLLs(object, useImps, omit.imps = omit.imps)
  ## calculate m log-likelihoods under pooled H1 estimates
  if (is.null(h1)) {
    LL1 <- getLLs(object, useImps, saturated = TRUE, omit.imps = omit.imps)
  } else {
    LL1 <- getLLs(h1, useImps, omit.imps = omit.imps)
  }
  #FIXME: check whether LL1 or LL0 returned errors?  add try()?

  ## calculate average of m LRTs
  LRT_con <- mean(-2*(LL0 - LL1)) # getLLs() already applies [useImps]
  ## average chisq across imputations
  if (is.null(h1)) {
    LRT_bar <- mean(sapply(object@testList[useImps], function(x) x$standard$stat))
  } else {
    LRT_bar <- mean(sapply(object@testList[useImps], function(x) x$standard$stat) -
                      sapply(h1@testList[useImps], function(x) x$standard$stat))
  }
  ## calculate average relative increase in variance
  a <- DF*(m - 1)
  ariv <- ((m + 1) / a) * (LRT_bar - LRT_con)
  test.stat <- ifelse(DF == 0, yes = 0, no = LRT_con / (DF*(1 + ariv)) )
  if (is.na(test.stat)) stop('D3 test statistic could not be calculated. ',
                             'Try the D2 pooling method.') #FIXME: check whether model-implied Sigma is NPD
  if (test.stat < 0) {
    message('Negative test statistic set to zero \n')
    test.stat <- 0
  }
  if (asymptotic) {
    out <- c("chisq" = test.stat * DF, df = DF,
             pvalue = pchisq(test.stat * DF, df = DF, lower.tail = FALSE),
             ariv = ariv, fmi = ariv / (1 + ariv))
  } else {
    ## F statistic
    if (a > 4) {
      v4 <- 4 + (a - 4) * (1 + (1 - (2 / a))*(1 / ariv))^2 # Enders (eq. 8.34)
    } else {
      v4 <- a*(1 + 1/DF)*(1 + 1/ariv)^2 / 2 # Enders (eq. 8.35)
      # v4 <- (DF + 1)*(m - 1)*(1 + (1 / ariv))^2 / 2 # Grund et al. (eq. 9)
    }
    out <- c("F" = test.stat, df1 = DF, df2 = v4,
             pvalue = pf(test.stat, df1 = DF, df2 = v4, lower.tail = FALSE),
             ariv = ariv, fmi = ariv / (1 + ariv))
  }
  ## add log-likelihood and AIC/BIC for target model
  if (is.null(h1)) {
    PT <- parTable(object)
    npar <- max(PT$free) - sum(PT$op == "==")
    out <- c(out, npar = npar, ntotal = N,
             #FIXME: omit? (also info criteria?)
             logl = mean(LL0), unrestricted.logl = mean(LL1),
             aic = -2*mean(LL0) + 2*npar, bic = -2*mean(LL0) + npar*log(N),
             bic2 = -2*mean(LL0) + npar*log((N + 2) / 24))
    ## NOTE: Mplus reports the average of m likelihoods evaluated at the
    ##       m point estimates, not evaluated at the pooled point estimates.
    ##       Mplus also uses those to calcluate AIC and BIC.
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}


##' @importFrom stats pf pchisq
##' @importFrom lavaan lavListInspect lavTestLRT parTable
D4.LRT <- function(object, h1 = NULL, useImps, asymptotic = FALSE,
                   omit.imps = c("no.conv","no.se")) {

  m <- length(useImps)

  if (is.null(h1)) {
    DF <- object@testList[[ useImps[1] ]]$standard[["df"]]
  } else {
    DF1 <- h1@testList[[ useImps[1] ]]$standard[["df"]]
    DF0 <- object@testList[[ useImps[1] ]]$standard[["df"]]
    DF <- DF0 - DF1
    if (DF < 0) stop('The "object" model must be nested within (i.e., have ',
                     'fewer degrees of freedom than) the "h1" model.')
  }

  ## Stack imputed data
  stackImps <- do.call(rbind, object@DataList[useImps])

  ## isolate @lavListCall arguments passed to lavaan()
  lavListCall <- object@lavListCall[-1]
  lavListArgNames <- names(formals(lavaan::lavaanList))
  lavArgNames <- setdiff(names(lavListCall), lavListArgNames)
  noNames <- which(nchar(lavArgNames) == 0L)
  if (length(noNames)) lavArgNames <- lavArgNames[-noNames]
  lavArgs <- c(list(model = parTable(object), data = stackImps),
               lavListCall[lavArgNames])
  ## simplify call: no SEs, standard test= (to use lavTestLRT's chisq)
  lavArgs$se <- "none"
  lavArgs$test <- "standard"

  ## fit model(s) to stacked data
  fit0 <- do.call(lavaan::lavaan, lavArgs)
  LL0_stacked <- fitMeasures(fit0, "logl")[[1]] / m
  if (!is.null(h1)) {
    lavArgs$model <- parTable(h1)
    fit1 <- do.call(lavaan::lavaan, lavArgs)
    LL1_stacked <- fitMeasures(fit1, "logl")[[1]] / m
  } else LL1_stacked <- fitMeasures(fit0, "unrestricted.logl")[[1]] / m
  LRT_stacked <- 2*(LL1_stacked - LL0_stacked)

  ## get average LRT across imputations for ARIV
  LRT_mean <- D2.LRT(object, h1 = h1, useImps = useImps,
                     return.mean.chisq = TRUE)

  ## calculate average relative increase in variance (Grund eq. 9)
  a <- DF*(m - 1)
  ariv <- max(0, ((m + 1) / a) * (LRT_mean - LRT_stacked))
  test.stat <- ifelse(DF == 0, yes = 0, no = LRT_stacked / (DF*(1 + ariv)) ) # Grund eq. 8

  if (asymptotic) {
    out <- c("chisq" = test.stat * DF, df = DF,
             pvalue = pchisq(test.stat * DF, df = DF, lower.tail = FALSE),
             ariv = ariv, fmi = ariv / (1 + ariv))
  } else {
    v4 <- a * (1 + 1/ariv)^2 # Grund eq. 10
    ## F statistic
    out <- c("F" = test.stat, df1 = DF, df2 = v4,
             pvalue = pf(test.stat, df1 = DF, df2 = v4, lower.tail = FALSE),
             ariv = ariv, fmi = ariv / (1 + ariv))
  }
  ## add log-likelihood and AIC/BIC for target model
  if (is.null(h1)) {
    PT <- parTable(object)
    npar <- max(PT$free) - sum(PT$op == "==") #FIXME: will Yves change this?
    N <- lavListInspect(object, "ntotal")
    LL0 <- getLLs(object, useImps, omit.imps = omit.imps)
    LL1 <- getLLs(object, useImps, saturated = TRUE, omit.imps = omit.imps)

    out <- c(out, npar = npar, ntotal = N,
             #FIXME: omit? (also info criteria?)
             logl = mean(LL0), unrestricted.logl = mean(LL1),
             aic = -2*mean(LL0) + 2*npar, bic = -2*mean(LL0) + npar*log(N),
             bic2 = -2*mean(LL0) + npar*log((N + 2) / 24))
    ## NOTE: Mplus reports the average of m likelihoods evaluated at the
    ##       m point estimates, not evaluated at the pooled point estimates.
    ##       Mplus also uses those to calculate AIC and BIC.
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}

##' @importFrom stats pchisq
##' @importFrom lavaan lavListInspect
robustify <- function(ChiSq, object, h1 = NULL, baseline = FALSE, useImps,
                      harm = TRUE, # harmonic or arithmetic mean of scaling.factor
                      LRTargs = list(), # for D2.LRT() when h1 + shifted
                      standard.test = "standard", scaled.test = "default") {
  scaleshift <- scaled.test == "scaled.shifted"

  if (baseline) {
    TEST <- lapply(object@baselineList[useImps], "[[", i = "test")
  } else TEST <- object@testList[useImps]

  d0 <- mean(sapply(TEST, function(x) x[[scaled.test]][["df"]]))
  if (harm) {
    c0 <- 1/ mean(1/sapply(TEST, function(x) x[[scaled.test]][["scaling.factor"]]))
  } else c0 <- mean(sapply(TEST, function(x) x[[scaled.test]][["scaling.factor"]]))


  if (!is.null(h1)) {
    ## MODEL COMPARISON

    d1 <- mean(sapply(h1@testList[useImps], function(x) x[[scaled.test]][["df"]]))
    if (scaleshift) {
      ScSh <- D2.LRT(object = object, h1 = h1, useImps = useImps, asymptotic = TRUE,
                     return.scale.shift = TRUE, LRTargs = LRTargs,
                     standard.test = standard.test, scaled.test = scaled.test,
                     pool.robust = TRUE) # just so D2.LRT() doesn't skip it
      delta_c <- if (harm) 1/mean(1/ScSh$scale) else mean(ScSh$scale)
      shift   <- mean(ScSh$shift)

    } else {
      if (harm) {
        c1 <- 1/mean(1/sapply(h1@testList[useImps],
                              function(x)  x[[scaled.test]][["scaling.factor"]]))
      } else c1 <- mean(sapply(h1@testList[useImps],
                               function(x) x[[scaled.test]][["scaling.factor"]]))
      delta_c <- ifelse( (d0-d1) == 0, 0, no = (d0*c0 - d1*c1) / (d0 - d1) )
      shift   <- 0
    }

    ChiSq["chisq.scaled"] <- ChiSq[["chisq"]] / delta_c + shift
    ChiSq["df.scaled"] <- d0 - d1
    ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                     df = ChiSq[["df.scaled"]],
                                     lower.tail = FALSE)
    ChiSq["chisq.scaling.factor"] <- delta_c
    if (scaleshift) ChiSq["chisq.shift.parameter"] <- shift

  } else {
    ## TEST SINGLE MODEL (against saturated h1)

    ChiSq["chisq.scaled"] <- ChiSq[["chisq"]] / c0
    ChiSq["df.scaled"] <- d0
    if (scaleshift) {
      ## add average shift parameter
      shift <- mean(sapply(TEST, function(x) {
        #FIXME: sum() no longer necessary, since lavaan returns 1 sum for MG-SEM
        sum(x[[scaled.test]][["shift.parameter"]])
      }))
      ChiSq["chisq.scaled"] <- ChiSq[["chisq.scaled"]] + shift
      ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                       df = ChiSq[["df.scaled"]],
                                       lower.tail = FALSE)
      ChiSq["chisq.scaling.factor"] <- c0
      ChiSq["chisq.shift.parameter"] <- shift
    } else {
      ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                       df = ChiSq[["df.scaled"]],
                                       lower.tail = FALSE)
      ChiSq["chisq.scaling.factor"] <- c0
    }
  }
  ChiSq
}

## Pools the LRT for 1 model or when comparing a single pair of models.
## (formerly lavTestLRT.mi)
##' @importFrom lavaan lavListInspect lavTestLRT
pairwiseLRT <- function(object, h1 = NULL, pool.method = c("D4","D3","D2"),
                        omit.imps = c("no.conv","no.se"), asymptotic = FALSE,
                        standard.test = "standard", scaled.test = "default",
                        pool.robust = FALSE, ...) {
  ## this also checks the class
  useImps <- imps2use(object = object, omit.imps = omit.imps)
  m <- length(useImps)

  DF0 <- object@testList[[ useImps[1] ]]$standard[["df"]]

  ## model comparison?
  if (!is.null(h1)) {
    if (!inherits(h1, "lavaan.mi")) stop("h1 is not class 'lavaan.mi'")
    if (!all(lavListInspect(object, "options")$test == lavListInspect(h1, "options")$test)) {
      stop('Different (sets of) test statistics were requested for the 2 models.')
    }

    useImps1 <- imps2use(object = h1, omit.imps = omit.imps)

    h0.not.h1 <- setdiff(useImps, useImps1)
    if (length(h0.not.h1)) warn0 <- paste('\n\nFor the following imputations, ',
                                          '"object" converged but not "h1":',
                                          paste(h0.not.h1, collapse = ", "),
                                          '\n\n')
    h1.not.h0 <- setdiff(useImps1, useImps)
    if (length(h1.not.h0)) warn1 <- paste('\n\nFor the following imputations, ',
                                          '"h1" converged but not "object":',
                                          paste(h1.not.h0, collapse = ", "),
                                          '\n\n')
    if (length(c(h0.not.h1, h1.not.h0)))
      warning('The models being compared did not converge on the same set of ',
              'imputations. ',
              if (length(h0.not.h1)) warn0 else '',
              if (length(h1.not.h0)) warn1 else '',
              'Likelihood ratio test conducted using only the ',
              'imputations for which both models converged.')
    useImps <- intersect(useImps, useImps1)
    m <- length(useImps) # yes, redefine if necessary
    if (m == 0L) stop('For no imputations did both models converge.')

    ## check DF
    DF1 <- h1@testList[[ useImps[1] ]]$standard[["df"]]
    if (DF0 == DF1) stop("models have equal degrees of freedom")
    if (DF0 < DF1) {
      H0 <- h1
      h1 <- object
      object <- H0
      H0 <- DF1
      DF1 <- DF0
      DF0 <- H0
    }
    DF <- DF0 - DF1
  } else DF <- DF0

  ## arguments passed to lavTestLRT (if D2)
  dots <- list(...)
  ## only keep relevant arguments
  keepArgs <- intersect(names(dots), names(formals(lavTestLRT)))
  if (length(keepArgs)) dots <- dots[keepArgs]

  ## check for lavTestLRT() arguments that imply standard.test= or scaled.test=
  if ("type" %in% keepArgs) {
    if (any(grepl(pattern = "browne", dots$type))) {
      standard.test <- dots$type
      pool.method <- "D2" #TODO: calculate from pooled moments?
    }
  }
  if ("test" %in% keepArgs) {
    if (dots$test[1] %in% c("satorra.bentler",
                            "yuan.bentler", "yuan.bentler.mplus",
                            "mean.var.adjusted", "scaled.shifted")) {
      scaled.test <- dots$test[1]
    }
  }

  ## check test-pooling options, for backward compatibility?
  pool.method <- tolower(pool.method[1])
  if (pool.method == "mplus") {
    pool.method <- "D3"
    asymptotic <- TRUE
  }
  if (tolower(pool.method) %in% c("cm","chan.meng","new.lrt","d4"))      pool.method <- "D4"
  if (tolower(pool.method) %in% c("mr","meng.rubin","old.lrt","d3"))     pool.method <- "D3"
  if (tolower(pool.method) %in% c("lmrr","li.et.al","pooled.wald","d2")) pool.method <- "D2"
  if (pool.method %in% c("D3","D4") && !lavListInspect(object, "options")$estimator %in% c("ML","PML","FML")) {
    message('"D3" and "D4" only available using maximum likelihood estimation. ',
            'Changed to pool.method = "D2".')
    pool.method <- "D2"
  }

  ## check for robust
  TEST_slot <- object@testList[[ useImps[1] ]]
  test.names <- unname(sapply(TEST_slot, "[[", "test"))
  if (standard.test[1] %in% c("standard", "browne.residual.nt",
                              "browne.residual.nt.model",
                              "browne.residual.adf",
                              "browne.residual.adf.model")) {
    standard.test <- standard.test[1]
  } else {
    standard.test <- intersect(test.names, c("standard", "browne.residual.nt",
                                             "browne.residual.nt.model",
                                             "browne.residual.adf",
                                             "browne.residual.adf.model"))[1]
  }

  # any non-standard test stats?
  robust <- FALSE
  if (length(test.names) > 1L) {
    checkMVadj <- scaled.test == "default"
    ## remove standard and any bootstrapped tests
    rm.idx <- which(test.names %in% c("none","default","standard",
                                      "browne.residual.nt",
                                      "browne.residual.nt.model",
                                      "browne.residual.adf",
                                      "browne.residual.adf.model"))
    if (length(rm.idx) > 0L) test.names <- test.names[-rm.idx]
    ## only acknowledge one scaled test statistic
    if (length(test.names) > 1L) {
      if (scaled.test[1] %in% c("satorra.bentler",
                                "yuan.bentler", "yuan.bentler.mplus",
                                "mean.var.adjusted", "scaled.shifted")) {
        scaled.test <- scaled.test[1]
      } else {
        scaled.test <- intersect(test.names,
                                 c("satorra.bentler",
                                   "yuan.bentler", "yuan.bentler.mplus",
                                   "mean.var.adjusted", "scaled.shifted"))[1]
      }
      ## did user request lavTestLRT(scaled.shifted = FALSE)?
      if (scaled.test == "scaled.shifted" && checkMVadj &&
          "scaled.shifted" %in% keepArgs) {
        if (!dots$scaled.shifted) scaled.test <- "mean.var.adjusted"
      }

    } else if (length(test.names) == 1L) {
      scaled.test <- test.names
    } else scaled.test <- character(0)

    robust <- length(scaled.test) > 0L
  }
  ## check alternative tests (method= or type=), which would imply !robust
  if ("method" %in% keepArgs) {
    if (dots$method == "standard") robust <- FALSE
  }
  if ("type" %in% keepArgs) {
    if (any(grepl(pattern = "browne", dots$type))) robust <- FALSE
  }

  if (!robust) pool.robust <- FALSE

  if (robust && !pool.robust && !asymptotic) {
    message('Robust correction can only be applied to pooled chi-squared ',
            'statistic, not F statistic. "asymptotic" was switched to TRUE.')
    asymptotic <- TRUE
  }

  if (pool.robust && pool.method %in% c("D3","D4")) {
    message('pool.robust = TRUE is only applicable when pool.method = "D2".\n',
            'Changed to pool.method = "D2".')
    pool.method <- "D2"
  }

  ## calculate pooled test:
  if (robust && pool.robust) {
    ## pool both the naive and robust test statistics, return both to
    ## make output consistent across options
    out.naive <- D2.LRT(object, h1 = h1, useImps = useImps,
                        asymptotic = asymptotic, pool.robust = FALSE,
                        standard.test = standard.test, scaled.test = scaled.test)
    out.robust <- D2.LRT(object, h1 = h1, useImps = useImps, LRTargs = dots,
                         asymptotic = asymptotic, pool.robust = TRUE,
                         standard.test = standard.test, scaled.test = scaled.test)
    out <- c(out.naive, out.robust)
  } else if (pool.method == "D2") {
    out <- D2.LRT(object, h1 = h1, useImps = useImps, asymptotic = asymptotic,
                  pool.robust = pool.robust,
                  standard.test = standard.test, scaled.test = scaled.test)
  } else if (pool.method == "D3") {
    out <- D3.LRT(object, h1 = h1, useImps = useImps, asymptotic = asymptotic,
                  omit.imps = omit.imps)
  } else if (pool.method == "D4") {
    out <- D4.LRT(object, h1 = h1, useImps = useImps, asymptotic = asymptotic,
                  omit.imps = omit.imps)
  }

  ## If test statistic is negative, return without any indices or robustness
  # if (asymptotic) {
  #   if (out[["chisq"]] == 0) {
  #     message('Negative pooled test statistic was set to zero, so fit will ',
  #             'appear to be arbitrarily perfect. ',
  #             if (robust) 'Robust corrections uninformative, not returned.',
  #             '\n')
  #     class(out) <- c("lavaan.vector","numeric")
  #     return(out)
  #   }
  # } else {
  #   if (out[["F"]] == 0) {
  #     message('Negative pooled test statistic was set to zero, so fit will ',
  #             'appear to be arbitrarily perfect.\n')
  #     class(out) <- c("lavaan.vector","numeric")
  #     return(out)
  #   }
  # }

  ## If robust statistics were not pooled above, robustify naive statistics
  if (robust & !pool.robust) {
    out <- robustify(ChiSq = out, object = object, h1 = h1, useImps = useImps,
                     LRTargs = dots, # in case scaled.shifted (same for below)
                     standard.test = standard.test, scaled.test = scaled.test)
    #TODO: default looks smarter, remove unnecessary distracting warning below
    # if (scaleshift) {
    #   extraWarn <- ' and shift parameter'
    # } else if (any(test.names == "mean.var.adjusted")) {
    #   extraWarn <- ' and degrees of freedom'
    # } else extraWarn <- ''
    # message('Robust corrections are made by pooling the standard chi-squared ',
    #         'statistic across ', m, ' imputations for which the model ',
    #         'converged, then applying the average (across imputations) scaling',
    #         ' factor', extraWarn, ' to that pooled value. \n',
    #         'To instead pool the robust test statistics, set pool.method = "D2"',
    #         ' and pool.robust = TRUE. \n')
  }

  class(out) <- c("lavaan.vector","numeric")
  ## add header
  attr(out, "header") <- paste("Test statistic(s) pooled using the",
                               pool.method, "pooling method")
  out
}

## Iterates over > 2 models to apply pairwiseLRT() sequentially
## (borrowed from semTools::compareFit, code for @nested)
multipleLRT <- function(..., argsLRT = list(), pool.method = c("D4","D3","D2"),
                        omit.imps = c("no.conv","no.se"), asymptotic = FALSE,
                        standard.test = "standard", scaled.test = "default",
                        pool.robust = FALSE) {
  ## separate models from lists of models
  mods <- list(...)

  ## check for lavaan models
  not.lavaan <- !sapply(mods, inherits, what = "lavaan.mi")
  if (any(not.lavaan)) stop("The following are not fitted lavaan.mi models:\n",
                            paste0(names(which(not.lavaan)), collapse = ", "))
  ## check whether any models failed to converge on any imputations
  nonConv <- !sapply(mods, function(fit) {
    any(sapply(fit@convergence, "[", i = "converged"))
  })
  if (all(nonConv)) {
    stop('No models converged')
  } else if (any(nonConv)) {
    message('The following models did not converge on any imputations, ',
            'so they are ignored:\n',
            paste(names(nonConv)[nonConv], collapse = ",\t"))
    mods <- mods[which(!nonConv)]
  }


  ## order models by increasing df (least-to-most constrained)
  DF <- sapply(mods, function(x) x@testList[[1]]$standard$df)
  mods <- mods[order(DF)]

  ## compare each sequential pair
  modsA <- mods[-1]
  modsB <- mods[-length(mods)]
  fitDiff <- list()

  for (i in seq_along(modsA)) {
    fitA <- modsA[[i]]
    fitB <- modsB[[i]]

    tempDiff <- do.call(pairwiseLRT,
                        c(list(object        = fitA,
                               h1            = fitB,
                               pool.method   = pool.method,
                               omit.imps     = omit.imps,
                               asymptotic    = asymptotic,
                               pool.robust   = pool.robust,
                               standard.test = standard.test,
                               scaled.test   = scaled.test),
                          argsLRT))

    if (names(tempDiff)[1] == "F") {
      statNames <-  c("F", "df1", "df2", "pvalue")
    } else statNames <- c("chisq", "df", "pvalue")
    ## check for scaled
    if (any(grepl(pattern = "scaled", x = names(tempDiff)))) {
      statNames <- paste0(statNames, ".scaled")
    }

    diffName <- paste(names(modsA)[i], "-", names(modsB)[i])
    fitDiff[[diffName]] <- tempDiff[ c(statNames, "ariv", "fmi") ]
  }

  out <- as.data.frame(do.call(rbind, fitDiff))
  class(out) <- c("lavaan.data.frame", "data.frame")
  out
}




