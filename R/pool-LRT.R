### Terrence D. Jorgensen & Yves Rosseel
### Last updated: 20 December 2022
### Pooled likelihood ratio test for multiple imputations
### Borrowed source code from lavaan/R/lav_test_LRT.R


#TODO:  whenever Yves adds scale/shift as attribute to anova output...


## -------------
## Main function
## -------------

##' Likelihood Ratio Test for Multiple Imputations
##'
##' Likelihood ratio test (LRT) for lavaan models fitted to multiple imputed
##' data sets.
##'
##' The \code{"D2"} method is available using any estimator and test statistic.
##' When using a likelihood-based estimator, 2 additional methods are available
##' to pool the LRT.
##' \itemize{
##'   \item The Meng & Rubin (1992) method, commonly referred to as \code{"D3"}.
##'         This method has many problems, discussed in Chan & Meng (2022).
##'   \item The Chan & Meng (2022) method, referred to as \code{"D4"} by
##'         Grund et al. (2021 preprint), resolves problems with \code{"D3"}.
##' }
##' When \code{"D2"} is not explicitly requested in situations it is the only
##' applicable method, (e.g., DWLS for categorical outcomes), users are notified
##' that \code{test} was set to \code{"D2"}.
##'
##' \code{test = "Mplus"} implies \code{"D3"} and \code{asymptotic = TRUE}
##' (see Asparouhov & Muthen, 2010).
##'
##' Note that unlike \code{\link[lavaan]{lavTestLRT}}, \code{lavTestLRT.mi} can
##' only be used to compare a single pair of models, not a longer list of
##' models.  To compare several nested models fitted to multiple imputations,
##' see examples on the \code{\link[semTools]{compareFit}} help page, which can
##' be called via the \code{anova} method (see \code{\linkS4class{lavaan.mi}}).
##'
##' @aliases lavTestLRT.mi
##' @importFrom lavaan lavTestLRT
##'
##' @param object,h1 An object of class \code{\linkS4class{lavaan.mi}}.
##'   \code{object} should be nested within (more constrained than) \code{h1}.
##' @param ... Additional arguments passed to \code{\link[lavaan]{lavTestLRT}},
##'   only if \code{test = "D2"} and \code{pool.robust = TRUE}
##' @param modnames Optional \code{character} of model names to use as row names
##'   in the resulting matrix of results (when more than 2 models are compared)
##' @param asANOVA \code{logical} indicating whether to return an object of
##'   class \code{"anova"}.  If \code{FALSE}, a numeric vector is returned for
##'   one (pair of) model(s), or a \code{data.frame} is returned for multiple
##'   pairs of models.
##' @param test \code{character} indicating which pooling method to use.
##'   \itemize{
##'     \item \code{"D4"}, \code{"new.LRT"}, \code{"cm"}, or \code{"chan.meng"}
##'           requests the method described by Chan & Meng (2022).
##'           This is currently the default.
##'     \item \code{"D3"}, \code{"old.LRT"}, \code{"mr"}, or \code{"meng.rubin"}
##'           requests the method described by Meng & Rubin (1992).
##'     \item \code{"D2"}, \code{"LMRR"}, or \code{"Li.et.al"} requests the
##'           complete-data LRT statistic should be calculated using each
##'           imputed data set, which will then be pooled across imputations, as
##'           described in Li, Meng, Raghunathan, & Rubin (1991).
##'   }
##'   Find additional details in Enders (2010, chapter 8).
##'
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'   imputations from pooled results.  Can include any of
##'   \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (\code{"no.npd"}) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases. Specific imputation numbers can also be included in this
##'   argument, in case users want to  apply their own custom omission criteria
##'   (or simulations can use different numbers of imputations without
##'   redundantly refitting the model).
##' @param asymptotic \code{logical}. If \code{FALSE} (default), the pooled test
##'   will be returned as an \emph{F}-distributed statistic with numerator
##'   (\code{df1}) and denominator (\code{df2}) degrees of freedom.
##'   If \code{TRUE}, the pooled \emph{F} statistic will be multiplied by its
##'   \code{df1} on the assumption that its \code{df2} is sufficiently large
##'   enough that the statistic will be asymptotically \eqn{\chi^2} distributed
##'   with \code{df1}.
##' @param pool.robust \code{logical}. Ignored unless \code{test = "D2"} and a
##'   robust test was requested. If \code{pool.robust = TRUE}, the robust test
##'   statistic is pooled, whereas \code{pool.robust = FALSE} will pool
##'   the naive test statistic (or difference statistic) and apply the average
##'   scale/shift parameter to it (unavailable for mean- and variance-adjusted
##'   difference statistics, so \code{pool.robust} will be set \code{TRUE}).
##'
##' @return
##'   A vector containing the LRT statistic (either an \code{F} or \eqn{\chi^2}
##'   statistic, depending on the \code{asymptotic} argument), its degrees of
##'   freedom (numerator and denominator, if \code{asymptotic = FALSE}), its
##'   \emph{p} value, and 2 missing-data diagnostics: the relative increase
##'   in variance (RIV, or average for multiparameter tests: ARIV) and the
##'   fraction missing information (FMI = ARIV / (1 + ARIV)). Robust statistics
##'   will also include the average (across imputations) scaling factor and
##'   (if relevant) shift parameter(s), unless \code{pool.robust = TRUE}.
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##'   Based on source code for \code{\link[lavaan]{lavTestLRT}} by Yves Rosseel
##'
##' @references
##'   Asparouhov, T., & Muthen, B. (2010). \emph{Chi-square statistics
##'   with multiple imputation}. Technical Report. Retrieved from
##'   \url{http://www.statmodel.com/}
##'
##'   Chan, K. W., & Meng, X. L. (2022). Multiple improvements of multiple
##'   imputation likelihood ratio tests. \emph{Statistica Sinica, 32},
##'   1489--1514. \doi{10.5705/ss.202019.0314}
##'
##'   Enders, C. K. (2010). \emph{Applied missing data analysis}.
##'   New York, NY: Guilford.
##'
##'   Grund, S., Lüdtke, O., & Robitzsch, A. (2021). Pooling methods for
##'   likelihood ratio tests in multiply imputed data sets. \emph{PsyArXiv}
##'   preprint: \doi{10.31234/osf.io/d459g}
##'
##'   Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
##'   Significance levels from repeated \emph{p}-values with multiply-imputed
##'   data. \emph{Statistica Sinica, 1}(1), 65--92. Retrieved from
##'   \url{https://www.jstor.org/stable/24303994}
##'
##'   Meng, X.-L., & Rubin, D. B. (1992). Performing likelihood ratio tests with
##'   multiply-imputed data sets. \emph{Biometrika, 79}(1), 103--111.
##'   \doi{10.2307/2337151}
##'
##'   Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}.
##'   New York, NY: Wiley.
##'
##' @seealso \code{\link[lavaan]{lavTestLRT}}, \code{\link[semTools]{compareFit}}
##'
##' @examples
##  \dontrun{
##' ## impose missing data for example
##' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
##'                                       "ageyr","agemo","school")]
##' set.seed(12345)
##' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
##' age <- HSMiss$ageyr + HSMiss$agemo/12
##' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
##'
##' ## impute missing data
##' library(Amelia)
##' set.seed(12345)
##' HS.amelia <- amelia(HSMiss, m = 20, noms = "school", p2s = FALSE)
##' imps <- HS.amelia$imputations
##'
##' ## specify CFA model from lavaan's ?cfa help page
##' HS.model <- '
##'   visual  =~ x1 + b1*x2 + x3
##'   textual =~ x4 + b2*x5 + x6
##'   speed   =~ x7 + b3*x8 + x9
##' '
##'
##' fit1 <- cfa.mi(HS.model, data = imps, estimator = "mlm")
##'
##' ## more constrained model: parallel indicators
##' HS.parallel <- '
##'   visual  =~ x1 + 1*x2 + 1*x3
##'   textual =~ x4 + 1*x5 + 1*x6
##'   speed   =~ x7 + 1*x8 + 1*x9
##' '
##'
##' fitp <- cfa.mi(HS.parallel, data = imps, estimator = "mlm")
##'
##' ## Even more constrained model: orthogonal factors
##' fit0 <- cfa.mi(HS.parallel, data = imps, estimator = "mlm",
##'                orthogonal = TRUE)
##'
##' ## By default, use D4.
##' ## Must request an asymptotic chi-squared statistic
##' ## in order to accommodate a robust correction.
##' lavTestLRT.mi(fit1, fit0, fitp, asymptotic = TRUE)
##'
##' ## Using D2, you can either robustify the pooled naive statistic ...
##' lavTestLRT.mi(fit1, fit0, asymptotic = TRUE, test = "D2")
##' ## ... or pool the robust chi-squared statistic
##' lavTestLRT.mi(fit1, fit0, asymptotic = TRUE, test = "D2",
##'               pool.robust = TRUE)
## }
##'
##' @export
lavTestLRT.mi <- function(object, ..., modnames = NULL, asANOVA = TRUE,
                          test = c("D4","D3","D2"),
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
    LRT.names <- intersect( names(dots) , names(formals(lavTestLRT)) )
    dots <- if (length(LRT.names)) dots[LRT.names] else NULL

  } else mods <- NULL

  ## compare 1 or multiple pairs?
  if (length(mods) < 2L) {
    argList <- c(list(object = object), dots,
                 list(test = test, omit.imps = omit.imps,
                      asymptotic = asymptotic, pool.robust = pool.robust))
    if (length(mods) == 1L) argList$h1 <- mods[[1]]
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
                 list(test = test, omit.imps = omit.imps,
                      asymptotic = asymptotic, pool.robust = pool.robust))
    results <- do.call(multipleLRT, argList)

    if (asANOVA) {
      ## get each model's pooled naive chi-squared
      outList <- lapply(modList, pairwiseLRT, asymptotic = TRUE, test = test,
                        omit.imps = omit.imps, pool.robust = pool.robust)
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


## ----------------
## Hidden Functions
## ----------------


##' @importFrom lavaan lavListInspect parTable lavTestLRT
D2.LRT <- function(object, h1 = NULL, useImps, asymptotic = FALSE,
                   return.mean.chisq = FALSE, # TRUE to ease D4
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
      c(chisq = out[2, "Chisq diff"], df = out[2, "Df diff"])
    }
    FIT <- eval(as.call(oldCall))
    ## check if there are any results
    noFit <- sapply(FIT@funList, function(x) x[1] == "fit failed")
    noLRT <- sapply(FIT@funList, function(x) x[1] == "lavTestLRT() failed")
    if (all(noFit | noLRT)) stop("No success using lavTestScore() on any imputations.")

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
  # lavaan 0.6-5: for now, we only acknowledge the first non-standard @test
  if (length(test.names) > 1L) {
    ## remove standard and any bootstrapped tests
    rm.idx <- which(test.names %in% c("standard","bootstrap","bollen.stine"))
    if (length(rm.idx) > 0L) {
      test.names <- test.names[-rm.idx]
    }
    ## only acknowledge the first scaled test statistic
    if (length(test.names) > 1L) {
      test.names <- test.names[1]
    }
  }

  ## pool Wald tests
  if (is.null(h1)) {
    test <- if (pool.robust) test.names[1] else "standard"
    DF <- mean(sapply(object@testList[useImps], function(x) x[[test]][["df"]]  ))
    w  <-      sapply(object@testList[useImps], function(x) x[[test]][["stat"]])
  } else {
    ## this will not get run if !pool.robust because logic catches that first
    DF0 <- mean(sapply(object@testList[useImps], function(x) x$standard[["df"]]))
    DF1 <- mean(sapply(    h1@testList[useImps], function(x) x$standard[["df"]]))
    DF <- DF0 - DF1
    w0 <- sapply(object@testList[useImps], function(x) x$standard[["stat"]])
    w1 <- sapply(    h1@testList[useImps], function(x) x$standard[["stat"]])
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
  ## for 1 model, add extra info (redundant if pool.robust)
  if (is.null(h1) & !pool.robust) {
    PT <- parTable(object)
    out <- c(out, npar = max(PT$free) - sum(PT$op == "=="),
             ntotal = lavListInspect(object, "ntotal"))
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}

##' @importFrom lavaan parTable lavaan lavListInspect
##' @importFrom methods getMethod
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
    fixedValues <- getMethod("coef","lavaan.mi")(object, type = "user",
                                                 omit.imps = omit.imps)
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
  test.stat <- LRT_con / (DF*(1 + ariv))
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
  test.stat <- LRT_stacked / (DF*(1 + ariv)) # Grund eq. 8

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
robustify <- function(ChiSq, object, h1 = NULL, baseline = FALSE, useImps) {
  test.names <- lavListInspect(object, "options")$test
  # lavaan 0.6-5: for now, we only acknowledge the first non-standard @test
  if (length(test.names) > 1L) {
    ## remove standard and any bootstrapped tests
    rm.idx <- which(test.names %in% c("standard","bootstrap","bollen.stine"))
    if (length(rm.idx) > 0L) {
      test.names <- test.names[-rm.idx]
    }
    ## only acknowledge the first scaled test statistic
    if (length(test.names) > 1L) {
      test.names <- test.names[1]
    }
  }
  scaleshift <- any(test.names %in% "scaled.shifted")

  if (baseline) {
    TEST <- lapply(object@baselineList[useImps], "[[", i = "test")
  } else TEST <- object@testList[useImps]

  d0 <- mean(sapply(TEST, function(x) x[[ test.names[1] ]][["df"]]))
  c0 <- mean(sapply(TEST, function(x) x[[ test.names[1] ]][["scaling.factor"]]))
  if (!is.null(h1)) {
    d1 <- mean(sapply(h1@testList[useImps], function(x) x[[ test.names[1] ]][["df"]]))
    c1 <- mean(sapply(h1@testList[useImps],
                      function(x) x[[ test.names[1] ]][["scaling.factor"]]))
    delta_c <- (d0*c0 - d1*c1) / (d0 - d1)
    ChiSq["chisq.scaled"] <- ChiSq[["chisq"]] / delta_c
    ChiSq["df.scaled"] <- d0 - d1
    ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                     df = ChiSq[["df.scaled"]],
                                     lower.tail = FALSE)
    ChiSq["chisq.scaling.factor"] <- delta_c
  } else {
    ChiSq["chisq.scaled"] <- ChiSq[["chisq"]] / c0
    ChiSq["df.scaled"] <- d0
    if (scaleshift) {
      ## add average shift parameter (or average of sums, if nG > 1)
      shift <- mean(sapply(TEST, function(x) sum(x[[ test.names[1] ]][["shift.parameter"]]) ))
      ChiSq["chisq.scaled"] <- ChiSq[["chisq.scaled"]] + shift
      ChiSq["pvalue.scaled"] <- pchisq(ChiSq[["chisq.scaled"]],
                                       df = ChiSq[["df.scaled"]],
                                       lower.tail = FALSE)
      ChiSq["chisq.scaling.factor"] <- c0
      ChiSq["chisq.shift.parameters"] <- shift
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
pairwiseLRT <- function(object, h1 = NULL, test = c("D4","D3","D2"),
                        omit.imps = c("no.conv","no.se"),
                        asymptotic = FALSE, pool.robust = FALSE, ...) {
  ## check class
  if (!inherits(object, "lavaan.mi")) stop("object is not class 'lavaan.mi'")

  useImps <- rep(TRUE, length(object@DataList))
  if ("no.conv" %in% omit.imps) useImps <- sapply(object@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(object@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.imps) {
    Heywood.lv <- sapply(object@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(object@convergence, "[[", i = "Heywood.ov")
    useImps <- useImps & !(Heywood.lv | Heywood.ov)
  }
  ## custom removal by imputation number
  rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
  if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
  ## whatever is left
  m <- sum(useImps)
  if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
  useImps <- which(useImps)

  DF0 <- object@testList[[ useImps[1] ]]$standard[["df"]]

  ## model comparison?
  if (!is.null(h1)) {
    if (!inherits(h1, "lavaan.mi")) stop("h1 is not class 'lavaan.mi'")
    if (!all(lavListInspect(object, "options")$test == lavListInspect(h1, "options")$test)) {
      stop('Different (sets of) test statistics were requested for the 2 models.')
    }

    useImps1 <- rep(TRUE, length(h1@DataList))
    if ("no.conv" %in% omit.imps) useImps1 <- sapply(h1@convergence, "[[", i = "converged")
    if ("no.se" %in% omit.imps) useImps1 <- useImps1 & sapply(h1@convergence, "[[", i = "SE")
    if ("no.npd" %in% omit.imps) {
      Heywood.lv <- sapply(h1@convergence, "[[", i = "Heywood.lv")
      Heywood.ov <- sapply(h1@convergence, "[[", i = "Heywood.ov")
      useImps1 <- useImps1 & !(Heywood.lv | Heywood.ov)
    }
    m <- sum(useImps1)
    if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
    useImps1 <- which(useImps1)

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

  ## only keep arguments relevant to pass to lavTestLRT (if D2)
  dots <- list(...)[names(formals(lavTestLRT))]

  ## check test-pooling options, for backward compatibility?
  test <- tolower(test[1])
  if (test == "mplus") {
    test <- "D3"
    asymptotic <- TRUE
  }
  if (tolower(test) %in% c("cm","chan.meng","new.lrt","d4")) test <- "D4"
  if (tolower(test) %in% c("mr","meng.rubin","old.lrt","d3")) test <- "D3"
  if (tolower(test) %in% c("lmrr","li.et.al","pooled.wald","d2")) test <- "D2"
  if (test %in% c("D3","D4") && !lavListInspect(object, "options")$estimator %in% c("ML","PML","FML")) {
    message('"D3" and "D4" only available using maximum likelihood estimation. ',
            'Changed test to "D2".')
    test <- "D2"
  }

  ## check for robust
  test.names <- lavListInspect(object, "options")$test
  # lavaan 0.6-5: for now, we only acknowledge the first non-standard @test
  if (length(test.names) > 1L) {
    ## remove standard and any bootstrapped tests
    rm.idx <- which(test.names %in% c("standard","bootstrap","bollen.stine"))
    if (length(rm.idx) > 0L) {
      test.names <- test.names[-rm.idx]
    }
    ## only acknowledge the first scaled test statistic
    #FIXME: loop over all tests? (e.g., "browne.residual.adf")
    if (length(test.names) > 1L) {
      test.names <- test.names[1]
    }
  }

  robust <- any(test.names %in% c("satorra.bentler","yuan.bentler",
                                  "yuan.bentler.mplus","scaled.shifted",
                                  "mean.var.adjusted","satterthwaite"))
  if (!robust) pool.robust <- FALSE

  scaleshift <- any(test.names == "scaled.shifted")
  if (scaleshift && !is.null(h1)) {
    if (test %in% c("D3","D4") | !pool.robust)
      #TODO: find scale/shift in lavaan's anova-table attributes
      message("If test = 'scaled.shifted' (estimator = 'WLSMV' or 'MLMV'), ",
              "model comparison is only available by (re)setting test = 'D2' ",
              "and pool.robust = TRUE.\n",
              "Control more options by passing arguments to lavTestLRT() via ",
              "the '...' argument.\n")
    pool.robust <- TRUE
    test <- 'D2'
  }

  if (robust && !pool.robust && !asymptotic) {
    message('Robust correction can only be applied to pooled chi-squared ',
            'statistic, not F statistic. "asymptotic" was switched to TRUE.')
    asymptotic <- TRUE
  }

  if (pool.robust && test %in% c("D3","D4")) {
    message('pool.robust = TRUE is only applicable when test = "D2". ',
            'Changed test to "D2".')
    test <- "D2"
  }

  ## calculate pooled test:
  if (robust && pool.robust) {
    ## pool both the naive and robust test statistics, return both to
    ## make output consistent across options
    out.naive <- D2.LRT(object, h1 = h1, useImps = useImps,
                        asymptotic = asymptotic, pool.robust = FALSE)
    out.robust <- D2.LRT(object, h1 = h1, useImps = useImps, LRTargs = dots,
                         asymptotic = asymptotic, pool.robust = TRUE)
    out <- c(out.naive, out.robust)
  } else if (test == "D2") {
    out <- D2.LRT(object, h1 = h1, useImps = useImps, asymptotic = asymptotic,
                  pool.robust = pool.robust)
  } else if (test == "D3") {
    out <- D3.LRT(object, h1 = h1, useImps = useImps, asymptotic = asymptotic,
                  omit.imps = omit.imps)
  } else if (test == "D4") {
    out <- D4.LRT(object, h1 = h1, useImps = useImps, asymptotic = asymptotic,
                  omit.imps = omit.imps)
  }

  ## If test statistic is negative, return without any indices or robustness
  if (asymptotic) {
    if (out[["chisq"]] == 0) {
      message('Negative pooled test statistic was set to zero, so fit will ',
              'appear to be arbitrarily perfect. ',
              if (robust) 'Robust corrections uninformative, not returned.',
              '\n')
      class(out) <- c("lavaan.vector","numeric")
      return(out)
    }
  } else {
    if (out[["F"]] == 0) {
      message('Negative pooled test statistic was set to zero, so fit will ',
              'appear to be arbitrarily perfect.\n')
      class(out) <- c("lavaan.vector","numeric")
      return(out)
    }
  }

  ## If robust statistics were not pooled above, robustify naive statistics
  if (robust & !pool.robust) {
    out <- robustify(ChiSq = out, object = object, h1 = h1, useImps = useImps)
    # if (scaleshift) {
    #   extraWarn <- ' and shift parameter'
    # } else if (any(test.names == "mean.var.adjusted")) {
    #   extraWarn <- ' and degrees of freedom'
    # } else extraWarn <- ''
    # message('Robust corrections are made by pooling the naive chi-squared ',
    #         'statistic across ', m, ' imputations for which the model ',
    #         'converged, then applying the average (across imputations) scaling',
    #         ' factor', extraWarn, ' to that pooled value. \n',
    #         'To instead pool the robust test statistics, set test = "D2" and ',
    #         'pool.robust = TRUE. \n')
  }

  class(out) <- c("lavaan.vector","numeric")
  out
}

## Iterates over > 2 models to apply pairwiseLRT() sequentially
## (borrowed from semTools::compareFit, code for @nested)
multipleLRT <- function(..., argsLRT = list(), test = c("D4","D3","D2"),
                        omit.imps = c("no.conv","no.se"),
                        asymptotic = FALSE, pool.robust = FALSE) {
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
                        c(list(object      = fitA,
                               h1          = fitB,
                               test        = test,
                               omit.imps   = omit.imps,
                               asymptotic  = asymptotic,
                               pool.robust = pool.robust),
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





