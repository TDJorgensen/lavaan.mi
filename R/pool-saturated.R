### Terrence D. Jorgensen
### Last updated: 7 March 2025
### pool saturated moments across imputations to fit SEM in "single" step:
###    Normal data: https://doi.org/10.3102/1076998612458320
###    Categorical: https://doi.org/10.1080/00273171.2018.1523000


##' Fit a Saturated `lavaan` Model to Multiple Imputed Data Sets
##'
##' This function fits a saturated model to a list of imputed data sets, and
##' returns a list of pooled summary statistics to treat as data.
##'
##'
##' @importFrom lavaan lavNames lavListInspect lavInspect parTable
##'
##' @param data A `list` of imputed data sets, or an object class from
##'   which imputed data can be extracted. Recognized classes are
##'   `lavaan.mi` (list of imputations stored in the `@@DataList` slot),
##'   `amelia` (created by the `Amelia` package), or
##'   `mids` (created by the `mice` package).
##' @param \dots
##'   Additional arguments passed to [lavaan::lavCor()] or to
##'   [lavaan.mi()].
##' @param return.fit
##'   `logical` indicating whether to return a [lavaan.mi-class]
##'   object containing the results of fitting the saturated model to multiple
##'   imputed `data`.  Could be useful for diagnostic purposes.
##' @param scale.W
##'   `logical`. If `TRUE` (default), the within- and
##'   between-imputation components will be pooled by scaling the
##'   within-imputation component by the ARIV (see Enders, 2010, p. 235,
##'   for definition and formula). Otherwise, the pooled matrix is
##'   calculated as the weighted sum of the within-imputation and
##'   between-imputation components (see Enders, 2010, ch. 8, for details).
##' @param omit.imps
##'   `character` vector specifying criteria for omitting
##'   imputations from pooled results of saturated model.  Can include any of
##'   `c("no.conv", "no.se", "no.npd")`, the first 2 of which are the
##'   default setting, which excludes any imputations that did not
##'   converge or for which standard errors could not be computed.  The
##'   last option (`"no.npd"`) would exclude any imputations which
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
##' @return If `return.fit=TRUE`, a [lavaan.mi-class] object.
##'   Otherwise, an object of class `lavMoments`, which is a `list`
##'   that contains at least `$sample.cov` and `$sample.nobs`,
##'   potentially also `$sample.mean`, `$sample.th`, `$NACOV`,
##'   and `$WLS.V`.  Also contains `$lavOptions` that will be passed
##'   to `lavaan(...)`.
##'
##' @note The `$lavOptions` list will always set `fixed.x=FALSE` and
##'   `conditional.x=FALSE`.  Users should not override those options when
##'   calling [lavaan::lavaan()] because doing so would yield
##'   incorrect *SE*s and test statistics. Computing the correct
##'   `$NACOV` argument would depend on which specific variables are
##'   treated as fixed, which would require an argument to `poolSat()` for
##'   users to declare names of exogenous variables. This has not yet been
##'   programmed, but that feature may be added in the future in order to reduce
##'   the number of parameters to estimate.
##'   However, if "exogenous" predictors were incomplete and imputed, then they
##'   are not truly fixed (i.e., unvarying across samples), so treating them as
##'   fixed would be illogical and yield biased *SE*s and test statistics.
##'
##'   The information returned by `poolSat()` must assume that any fitted
##'   SEM will include all the variables in `$sample.cov` and (more
##'   importantly) in `$NACOV`.  Although `lavaan` can drop unused
##'   rows/columns from `$sample.cov`, it cannot be expected to drop the
##'   corresponding sampling variances of those eliminated (co)variances from
##'   `$NACOV`.  Thus, it is necessary to use `poolSat()` to obtain
##'   the appropriate summary statistics for any particular SEM (see **Examples**).
##'
##' @seealso [lavaan.mi()] for traditional method (fit SEM to each imputation,
##'   pool results afterward).
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Lee, T., & Cai, L. (2012). Alternative multiple imputation inference for
##'   mean and covariance structure modeling.
##'   *Journal of Educational and Behavioral Statistics, 37*(6), 675--702.
##'   \doi{10.3102/1076998612458320}
##'
##'   Chung, S., & Cai, L. (2019). Alternative multiple imputation inference
##'   for categorical structural equation modeling,
##'   *Multivariate Behavioral Research, 54*(3), 323--337.
##'   \doi{10.1080/00273171.2018.1523000}
##'
##' @examples
##'
##' data(HS20imps) # import a list of 20 imputed data sets
##'
##' ## fit saturated model to imputations, pool those results
##' impSubset1 <- lapply(HS20imps, "[", i = paste0("x", 1:9)) # only modeled variables
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
##' ## By default, the "Scaled" column provides a "scaled.shifted" test
##' ## statistic that maintains an approximately nominal Type I error rate.
##' summary(fit, fit.measures = TRUE, standardized = "std.all")
##' ## Note that this scaled statistic does NOT account for deviations from
##' ## normality, because the default normal-theory standard errors were
##' ## requested when running poolSat().  See below about non-normality.
##'
##' ## Alternatively, "Browne's residual-based (ADF) test" is also available:
##' lavTest(fit, test = "browne.residual.adf", output = "text")
##'
##' ## Optionally, save the saturated-model lavaan.mi object, which
##' ## could be helpful for diagnosing convergence problems per imputation.
##' satFit <- poolSat(impSubset1, return.fit = TRUE)
##'
##'
##' ## FITTING MODELS TO DIFFERENT (SUBSETS OF) VARIABLES
##'
##' ## If you only want to analyze a subset of these variables,
##' mod.vis <- 'visual  =~ x1 + x2 + x3'
##' ## you will get an error:
##' try(
##'   fit.vis <- cfa(mod.vis, data = prePooledData) # error
##' )
##'
##' ## As explained in the "Note" section, you must use poolSat() again for
##' ## this subset of variables
##' impSubset3 <- lapply(HS20imps, "[", i = paste0("x", 1:3)) # only modeled variables
##' visData <- poolSat(impSubset3)
##' fit.vis <- cfa(mod.vis, data = visData) # no problem
##'
##'
##' ## OTHER lavaan OPTIONS
##'
##' \donttest{
##' ## fit saturated MULIPLE-GROUP model to imputations
##' impSubset2 <- lapply(HS20imps, "[", i = c(paste0("x", 1:9), "school"))
##' (prePooledData2 <- poolSat(impSubset2, group = "school",
##'                            ## request standard errors that are ROBUST
##'                            ## to violations of the normality assumption:
##'                            se = "robust.sem"))
##' ## Nonnormality-robust standard errors are implicitly incorporated into the
##' ## pooled weight matrix (NACOV= argument), so they are
##' ## AUTOMATICALLY applied when fitting the model:
##' fit.config <- cfa(HS.model, data = prePooledData2, group = "school",
##'                   std.lv = TRUE)
##' ## standard errors and chi-squared test of fit both robust to nonnormality
##' summary(fit.config)
##'
##'
##' ## CATEGORICAL OUTCOMES
##'
##' ## discretize the imputed data, for an example of 3-category data
##' HS3cat <- lapply(impSubset1, function(x) {
##'   as.data.frame( lapply(x, cut, breaks = 3, labels = FALSE) )
##' })
##' ## pool polychoric correlations and thresholds
##' (prePooledData3 <- poolSat(HS3cat, ordered = paste0("x", 1:9)))
##'
##' fitc <- cfa(HS.model, data = prePooledData3, std.lv = TRUE)
##' summary(fitc)
##'
##' ## Optionally, use unweighted least-squares estimation.  However,
##' ## you must first REMOVE the pooled weight matrix (WLS.V= argument)
##' ## or replace it with an identity matrix of the same dimensions:
##' prePooledData4 <- prePooledData3
##' prePooledData4$WLS.V <- NULL
##' ## or prePooledData4$WLS.V <- diag(nrow(prePooledData3$WLS.V))
##' fitcu <- cfa(HS.model, data = prePooledData4, std.lv = TRUE, estimator = "ULS")
##' ## Note that the SEs and test were still appropriately corrected:
##' summary(fitcu)
##' }
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
      loadNamespace("mice")
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

  ## apply omit.imps= criteria (this also checks the class is lavaan.mi)
  useImps <- imps2use(object = satFit, omit.imps = omit.imps)
  m <- length(useImps)


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
    WLS.V[[g]] <- solve(diag(diag(NACOV[[g]])))
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
                         test               = "scaled.shifted")
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



