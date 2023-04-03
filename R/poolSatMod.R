### Terrence D. Jorgensen
### Last updated: 3 April 2023
### pool saturated moments across imputations to fit SEM in "single" step:
###    Normal data: https://doi.org/10.3102/1076998612458320
###    Categorical: https://doi.org/10.1080/00273171.2018.1523000


##' Fit a Saturated lavaan Model to Multiple Imputed Data Sets
##'
##' This function fits a saturated model to a list of imputed data sets, and
##' returns a list of pooled summary statistics to treat as data.
##'
##'
##' @importFrom lavaan lavNames lavListInspect lavInspect parTable
##'
##' @param data A a \code{list} of imputed data sets, or an object class from
##'   which imputed data can be extracted. Recognized classes are
##'   \code{lavaan.mi} (stored in the \code{@@DataList} slot),
##'   \code{amelia} (created by the Amelia package), or
##'   \code{mids} (created by the mice package).
##' @param \dots
##'   Additional arguments passed to \code{\link[lavaan]{lavCor}} or to
##'   \code{\link{lavaan.mi}}.
##' @param return.fit
##'   \code{logical} indicating whether to return a \code{\linkS4class{lavaan.mi}}
##'   object containing the results of fitting the saturated model to multiple
##'   imputed \code{data}.  Could be useful for diagnostic purposes.
##' @param scale.W
##'   \code{logical}. If \code{TRUE} (default), the within- and
##'   between-imputation components will be pooled by scaling the
##'   within-imputation component by the ARIV (see Enders, 2010, p. 235,
##'   for definition and formula). Otherwise, the pooled matrix is
##'   calculated as the weighted sum of the within-imputation and
##'   between-imputation components (see Enders, 2010, ch. 8, for details).
##' @param omit.imps
##'   \code{character} vector specifying criteria for omitting
##'   imputations from pooled results of saturated model.  Can include any of
##'   \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (\code{"no.npd"}) would exclude any imputations which
##'   yielded a nonpositive definite covariance matrix for observed or
##'   latent variables, which would include any "improper solutions" such
##'   as Heywood cases.  NPD solutions are not excluded by default because
##'   they are likely to occur due to sampling error, especially in small
##'   samples.  However, gross model misspecification could also cause
##'   NPD solutions, users can compare pooled results with and without
##'   this setting as a sensitivity analysis to see whether some
##'   imputations warrant further investigation. Specific imputation
##'   numbers can also be included in this argument, in case users want to
##'   apply their own custom omission criteria (or simulation studies can use
##'   different numbers of imputations without redundantly refitting the model).
##'
##' @return If \code{return.fit=TRUE}, a \code{\linkS4class{lavaan.mi}} object.
##'   Otherwise, an object of class \code{lavMoments}, which is a \code{list}
##'   that contains at least \code{$sample.cov} and \code{$sample.nobs},
##'   potentially also \code{$sample.mean}, \code{$sample.th}, \code{$NACOV},
##'   and \code{$WLS.V}.  Also contains \code{$lavOptions} that will be passed
##'   to \code{\link[lavaan]{lavaan}(...)}.
##'
##' @note The \code{$lavOptions} list will always set \code{fixed.x=FALSE} and
##'   \code{conditional.x=FALSE}.  Users should not override those options when
##'   calling \code{\link[lavaan]{lavaan}} because doing so would yield
##'   incorrect \emph{SE}s and test statistics. Computing the correct
##'   \code{$NACOV} argument would depend on which specific variables are
##'   treated as fixed, which would require an argument to \code{poolSat()} for
##'   users to declare names of exogenous variables. This has not yet been
##'   programmed, but that feature may be added in the future in order to reduce
##'   the number of parameters to estimate.
##'   However, if "exogenous" predictors were incomplete and imputed, then they
##'   are not truly fixed (i.e., unvarying across samples), so treating them as
##'   fixed would be illogical and yield biased \emph{SE}s and test statistics.
##'
##'   The information returned by \code{poolSat()} must assume that any fitted
##'   SEM will include all the variables in \code{$sample.cov} and (more
##'   importantly) in \code{$NACOV}.  Although \code{lavaan} can drop unused
##'   rows/columns from \code{$sample.cov}, it cannot be expected to drop the
##'   corresponding sampling variances of those eliminated (co)variances from
##'   \code{$NACOV}.  Thus, it is necessary to use \code{poolSat()} to obtain
##'   the appropriate summary statistics for any particular SEM (see Examples).
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Lee, T., & Cai, L. (2012). Alternative multiple imputation inference for
##'   mean and covariance structure modeling.
##'   \emph{Journal of Educational and Behavioral Statistics, 37}(6), 675--702.
##'   \doi{10.3102/1076998612458320}
##'
##'   Chung, S., & Cai, L. (2019). Alternative multiple imputation inference
##'   for categorical structural equation modeling,
##'   \emph{Multivariate Behavioral Research, 54}(3), 323--337.
##'   \doi{10.1080/00273171.2018.1523000}
##'
##' @examples
##'
##' ## impose missing data for example
##' HSMiss <- HolzingerSwineford1939[ , c(paste0("x", 1:9),
##'                                       "ageyr","agemo","school")]
##' set.seed(123)
##' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
##' age <- HSMiss$ageyr + HSMiss$agemo/12
##' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
##'
##' ## impute missing data with Amelia
##' library(Amelia)
##' set.seed(456)
##' HS.amelia <- amelia(HSMiss, m = 20, noms = "school", p2s = FALSE)
##' imps <- HS.amelia$imputations
##'
##' ## fit saturated model to imputations, pool those results
##' impSubset1 <- lapply(imps, "[", i = paste0("x", 1:9)) # only modeled variables
##' (prePooledData <- poolSat(impSubset1))
##'
##' ## Note: no means were returned (default lavOption() is meanstructure=FALSE)
##' (prePooledData <- poolSat(impSubset1, meanstructure = TRUE))
##'
##' ## specify CFA model from lavaan's ?cfa help page
##' HS.model <- '
##'   visual  =~ x1 + x2 + x3
##'   textual =~ x4 + x5 + x6
##'   speed   =~ x7 + x8 + x9
##' '
##'
##' ## fit model to summary statistics in "prePooledData"
##' fit <- cfa(HS.model, data = prePooledData, std.lv = TRUE)
##'
##' ## In the summary() output, IGNORE the standard test statistic
##' ## reported directly under the heading "Model Test User Model:".
##' summary(fit)
##' ## Instead, refer to "Browne's residual-based (ADF) test".
##' ## Browne's (1984) test statistic better maintains the chosen Type I error rate.
##'
##' ## Browne's (1984) test statistic is also available via lavTest():
##' lavTest(fit, test = "browne.residual.adf", output = "text")
##'
##' ## Browne's test statistic will not be reported by fitMeasures(), which only
##' ## calculates fit indices using the standard (not Browne's) statistic.
##' fitMeasures(fit, output = "text") # missing Browne, interpret with caution!
##'
##'
##' ## optionally, save the saturated-model lavaan.mi object
##' satFit <- poolSat(impSubset1, return.fit = TRUE)
##' ## this could be helpful for diagnosing convergence problems per imputation
##'
##'
##' ## FITTING MODELS TO DIFFERENT (SUBSETS OF) VARIABLES
##'
##' ## If you only want to analyze a subset of these variables, you will get
##' ## an error:
##' \dontrun{
##' mod.vis <- 'visual  =~ x1 + x2 + x3'
##' fit.vis <- cfa(mod.vis, data = prePooledData) # error
##' }
##'
##' ## As explained in the "Note" section, you must use poolSat() again for
##' ## this subset of variables
##' impSubset3 <- lapply(imps, "[", i = paste0("x", 1:3)) # only modeled variables
##' visData <- poolSat(impSubset3)
##' fit.vis <- cfa(mod.vis, data = visData) # no problem
##'
##'
##' ## OTHER lavaan OPTIONS
##'
##' ## fit saturated MULIPLE-GROUP model to imputations
##' impSubset2 <- lapply(imps, "[", i = c(paste0("x", 1:9), "school"))
##' (prePooledData2 <- poolSat(impSubset2, group = "school", se = "robust.sem"))
##' ## Nonnormality accounted for in this step by requesting robust
##' ## standard errors.  These are used as weights in the next step:
##' fit.config <- cfa(HS.model, data = prePooledData2, group = "school",
##'                   std.lv = TRUE)
##' ## standard errors and Browne's chi-squared test both robust to nonnormality
##' summary(fit.config)
##'
##'
##' ## COMPARING NESTED MODELS
##'
##' ## Browne's (1984) residual-based statistics have not been generalized to
##' ## the case of nested-model comparisons.  Thus, the lavTestLRT() function
##' ## does not provide such an option, but would instead provide a standard
##' ## test that IGNORES between-imputation variance.
##'
##' ## fit more-constrained levels of measurement equivalence
##' fit.metric <- cfa(HS.model, data = prePooledData2, group = "school",
##'                   std.lv = TRUE, group.equal = "loadings")
##' fit.scalar <- cfa(HS.model, data = prePooledData2, group = "school",
##'                   std.lv = TRUE, group.equal = c("loadings","intercepts"))
##'
##' ## Cannot trust Type I error rate of standard LRT statistic
##' ## because it IGNORES between-imputation variance
##' lavTestLRT(fit.config, fit.metric, fit.scalar)
##'
##' ## Instead, request the chi-squared difference test using Browne's statistic
##' lavTestLRT(fit.config, fit.metric, fit.scalar, type = "browne.residual.adf")
##'
##'
##' ## CATEGORICAL OUTCOMES
##'
##' ## discretize imputed data
##' HS3cat <- lapply(impSubset1, function(x) {
##'   as.data.frame( lapply(x, cut, breaks = 3, labels = FALSE) )
##' })
##' ## pool polychoric correlations and thresholds
##' (prePooledData3 <- poolSat(HS3cat, ordered = paste0("x", 1:9)))
##'
##' fitc <- cfa(HS.model, data = prePooledData3, std.lv = TRUE)
##' summary(fitc) # again, only pay attention to Browne's ADF test
##'
##' ## Optionally, use unweighted least-squares estimation.
##' ## Note that the robust SEs and Browne's test will still be used.
##' fitcu <- cfa(HS.model, data = prePooledData3, std.lv = TRUE, estimator = "ULS")
##' summary(fitcu)
##'
##' @export
poolSat <- function(data, ..., return.fit = FALSE, scale.W = TRUE,
                    omit.imps = c("no.conv","no.se")) {
  ## Validate list of imputed data sets
  imputedData <- NULL
  if (is.data.frame(data)) {
    stop('lavaan.mi(data=) cannot be a single data.frame')

  } else if (is.list(data)) {

    ## check whether it is a mids or amelia object (both inherit from list)
    if (inherits(data, "mids")) {
      requireNamespace("mice")
      m <- data$m
      imputedData <- vector("list", m)
      for (i in 1:m) {
        imputedData[[i]] <- mice::complete(data, action = i, include = FALSE)
      }
    } else if (inherits(data, "amelia")) {
      imputedData <- data$imputations
      m <- length(imputedData)
      class(imputedData) <- "list" # override "mi" inheritance

    } else {
      ## if not mids or amelia, check it is a list of data.frames
      stopifnot(all(sapply(data, inherits, what = "data.frame")))
      imputedData <- data
      m <- length(imputedData)
      class(imputedData) <- "list" # override inheritance (e.g., "mi" if Amelia)
    }

  } else if (inherits(data, "lavaan.mi")) {
    imputedData <- data@DataList
    m <- length(imputedData)
  } else stop("data= must be a list of imputed data.frames or an object of ",
              "class 'lavaan.mi', 'mids', or 'amelia'")

  ## prepare dots to pass to lavCor()
  dots <- list(...)
  dots$output <- "fit"
  dots$object <- imputedData[[1]]
  ## fit template model to first imputation
  satTemp <- do.call(lavaan::lavCor, dots)
  ## prepare dots to pass to lavaan()
  dots$object   <- NULL
  dots$model    <- parTable(satTemp)
  dots$data     <- imputedData
  dots$h1       <- FALSE # redundant
  dots$baseline <- FALSE # unnecessary
  ## save NACOV per imputation
  dots$FUN <- function(obj) list(gamma   = lavaan::lavInspect(obj, "gamma"),
                                 wls.obs = lavaan::lavInspect(obj, "wls.obs") )
  ## get rid of estimates in parameter table
  dots$model$est   <- NULL
  dots$model$se    <- NULL
  dots$model$start <- NULL

  ## remove any lavCor()-specific arguments
  ## FIXME: use setdiff() of formals() to automate?
  if (!is.null(dots$ov.names.x    )) dots$ov.names.x     <- NULL
  if (!is.null(dots$cor.smooth    )) dots$cor.smooth     <- NULL
  if (!is.null(dots$cor.smooth.tol)) dots$cor.smooth.tol <- NULL
  dots$output <- NULL
  ## fit saturated model to imputations
  satFit <- do.call(lavaan.mi, dots)
  if (return.fit) return(satFit)

  ## apply omit.imps= criteria
  useImps <- rep(TRUE, length(satFit@DataList))
  if ("no.conv" %in% omit.imps) useImps <- sapply(satFit@convergence, "[[", i = "converged")
  if ("no.se" %in% omit.imps) useImps <- useImps & sapply(satFit@convergence, "[[", i = "SE")
  if ("no.npd" %in% omit.imps) {
    Heywood.lv <- sapply(satFit@convergence, "[[", i = "Heywood.lv")
    Heywood.ov <- sapply(satFit@convergence, "[[", i = "Heywood.ov")
    useImps <- useImps & !(Heywood.lv | Heywood.ov)
  }
  ## custom removal by imputation number
  rm.imps <- omit.imps[ which(omit.imps %in% 1:length(useImps)) ]
  if (length(rm.imps)) useImps[as.numeric(rm.imps)] <- FALSE
  ## whatever is left
  m <- sum(useImps)
  if (m == 0L) stop('No imputations meet "omit.imps" criteria.')
  useImps <- which(useImps)


  N             <- lavListInspect(satFit, "nobs")
  ngroups       <- lavListInspect(satFit, "ngroups")
  meanstructure <- lavListInspect(satFit, "meanstructure")
  categorical   <- lavListInspect(satFit, "categorical")

  ## average coefficients
  satMoments  <- fitted_lavaan_mi(satFit, omit.imps = omit.imps)
  if (ngroups > 1L) {
    sample.cov <- sapply(satMoments, "[[", i = "cov", simplify = FALSE)
    if (meanstructure) {
      sample.mean <- sapply(satMoments, "[[", i = "mean", simplify = FALSE)
    }
    if (categorical) {
      sample.th <- sapply(satMoments, "[[", i = "th", simplify = FALSE)
      th.idx <- lavInspect(satTemp, "th.idx") # names not in lavaanList
    }

  } else {
    ## store single group as list
    sample.cov <- list(satMoments$cov)
    if (meanstructure) sample.mean <- list(satMoments$mean)
    if (categorical) {
      sample.th <- list(satMoments$th)
      th.idx <- list(lavInspect(satTemp, "th.idx")) # names not in lavaanList
    }
  }

  ## extract and average NACOV (W component), pool with B = var(wls.obs)
  NACOV <- vector("list", length = ngroups)
  if (ngroups > 1L) {

    for (g in 1:ngroups) {
      B <- (N[g] - 1L)*cov(do.call(rbind, lapply(satFit@funList[useImps], function(x) x$wls.obs[[g]] )))
      W <- Reduce("+", lapply(satFit@funList[useImps], function(x) x$gamma[[g]])) / m #FIXME: don't use omit.imps= argument?
      if (scale.W) {
        inv.W <- try(solve(W), silent = TRUE)
        if (inherits(inv.W, "try-error")) {
          warning("Could not invert within-imputation covariance matrix. ",
                  "Generalized inverse used instead.\nIt may be ",
                  "safer to set `scale.W = FALSE' (and `asymptotic = TRUE').")
          inv.W <- MASS::ginv(W)
        }
        ## relative increase in variance due to missing data
        r <- (1 + 1/m) * mean(diag(B %*% inv.W)) # Enders (2010, p. 235) eqs. 8.20-21
        NACOV[[g]] <- (1 + r) * W

      } else {
        ## less reliable, but constraints prevent inversion of W
        NACOV[[g]] <- W + (1 + 1/m)*B
      }

    }

  } else {
    B <- (N - 1L)*cov(do.call(rbind, lapply(satFit@funList[useImps], "[[", i = "wls.obs")))
    W <- Reduce("+", lapply(satFit@funList[useImps], "[[", i = "gamma")) / m #FIXME: don't use omit.imps= argument?
    if (scale.W) {
      inv.W <- try(solve(W), silent = TRUE)
      if (inherits(inv.W, "try-error")) {
        warning("Could not invert within-imputation covariance matrix. ",
                "Generalized inverse used instead.\nIt may be ",
                "safer to set `scale.W = FALSE' (and `asymptotic = TRUE').")
        inv.W <- MASS::ginv(W)
      }
      ## relative increase in variance due to missing data
      r <- (1 + 1/m) * mean(diag(B %*% inv.W)) # Enders (2010, p. 235) eqs. 8.20-21
      NACOV[[1]] <- (1 + r) * W

    } else {
      ## less reliable, but constraints prevent inversion of W
      NACOV[[1]] <- W + (1 + 1/m)*B
    }

  }

  ## Invert NACOV first, then take diagonal for DWLS weights
  #FIXME: if (categorical) {}
  WLS.V <-  vector("list", length = ngroups)
  for (g in 1:ngroups) {
    WLS.V[[g]] <- diag(diag(solve(NACOV[[g]])))
    dimnames(WLS.V[[g]]) <- dimnames(NACOV[[g]])
    class(WLS.V[[g]]) <- c("lavaan.matrix.symmetric","matrix")
  }

  ## return list of moments
  if (ngroups > 1L) {
    out <- list(sample.cov = sample.cov, sample.nobs = N,
                NACOV = NACOV, WLS.V = WLS.V)
    if (meanstructure) out$sample.mean <- sample.mean
    if (categorical) {
      out$sample.th <- sample.th
      attr(out$sample.th, "th.idx") <- th.idx
    }

  } else {
    ## drop list for single group
    out <- list(sample.cov = sample.cov[[1]], sample.nobs = N,
                NACOV = NACOV[[1]], WLS.V = WLS.V[[1]])
    if (meanstructure) out$sample.mean <- sample.mean[[1]]
    if (categorical) {
      out$sample.th <- sample.th[[1]]
      attr(out$sample.th, "th.idx") <- th.idx[[1]]
    }
  }

  ## necessary?  lavaan sets this already when NACOV= or WLS.V= are provided
  out$ov.order <- "data"

  ## set recommended arguments (estimator, se, test)
  out$lavOptions <- list(sample.cov.rescale = FALSE,
                         fixed.x            = FALSE, #TODO: model-specific NACOV based on ov.x
                         conditional.x      = FALSE,
                         estimator          = ifelse(categorical, "DWLS","ML"),
                         se                 = "robust.sem",
                         test               = "Browne.residual.adf")
  ## set class and return
  class(out) <- c("lavMoments", "list")
  out
}


##' @exportS3Method print lavMoments
print.lavMoments <- function(x, ...) {
  nameX <- substitute(x)

  cat('This lavMoments-class object contains summary statistics in a list ',
      'consisting of the following elements:\n  ',
      sep = '')
  cat(names(x), sep = ', ')
  cat('\nYou can view list elements with str(', nameX, ') or standard ',
      'extraction methods (e.g., ', nameX, '$sample.cov).\n\n', sep = '')

  cat('This object can be passed to lavaan() as the data= argument, in which ',
      'case each element in this object (e.g., $sample.cov) will be internally ',
      'passed to the corresponding lavaan() argument.\n', sep = '')

  if ("lavOptions" %in% names(x)) {
    cat('The following recommended lavOptions() will also be passed to lavaan():\n  ',
        sep = '')
    cat(paste0(names(x$lavOptions),
               ifelse(sapply(x$lavOptions, is.character), ' = "', ' = '),
               x$lavOptions,
               ifelse(sapply(x$lavOptions, is.character), '"', '')),
        sep = ',\n  ')
    cat('\nYou can override these settings in your lavaan() call, but you ',
        'should NOT set fixed.x=TRUE, which would yield incorrect standard ',
        'errors and test statistics.\n\n',
        sep = '')
  }

  return(invisible(x))
}



