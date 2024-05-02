### Terrence D. Jorgensen
### Last updated: 3 May 2024
### Class and Methods for lavaan.mi object


##' Class for a lavaan Model Fitted to Multiple Imputations
##'
##' This class extends the [lavaanList-class] class, created by
##' fitting a lavaan model to a list of data sets. In this case, the list of
##' data sets are multiple imputations of missing data.
##'
##'
##' @name lavaan.mi-class
##' @importClassesFrom lavaan lavaanList
##' @aliases lavaan.mi-class  show,lavaan.mi-method  summary,lavaan.mi-method
##'   fitMeasures,lavaan.mi-method  fitmeasures,lavaan.mi-method
##'   nobs,lavaan.mi-method  coef,lavaan.mi-method  vcov,lavaan.mi-method
##'   fitted,lavaan.mi-method  fitted.values,lavaan.mi-method
##' @docType class
##'
##' @slot coefList `list` of estimated coefficients in matrix format (one
##'   per imputation) as output by `lavInspect(fit, "est")`
##' @slot phiList `list` of model-implied latent-variable covariance
##'   matrices (one per imputation) as output by
##'   `lavInspect(fit, "cov.lv")`
##' @slot miList `list` of modification indices output by
##'   [lavaan::modindices()]
##' @slot lavListCall call to [lavaan::lavaanList()] used to fit the
##'   model to the list of imputed data sets in `@DataList`, stored as a
##'   `list` of arguments
##' @slot convergence `list` of `logical` vectors indicating whether,
##'   for each imputed data set, (1) the model converged on a solution, (2)
##'   *SE*s could be calculated, (3) the (residual) covariance matrix of
##'   latent variables (\eqn{\Psi}) is non-positive-definite, and (4) the
##'   residual covariance matrix of observed variables (\eqn{\Theta}) is
##'   non-positive-definite.
##' @slot lavaanList_slots All remaining slots are from
##'   [lavaanList-class], but [lavaan.mi()] only populates a
##'   subset of the `list` slots, two of them with custom information:
##' @slot version Named `character` vector indicating the `lavaan` and
##'   `lavaan.mi` version numbers.
##' @slot DataList The `list` of imputed data sets
##' @slot SampleStatsList List of output from
##'   `lavInspect(fit, "sampstat")` applied to each fitted model.
##' @slot ParTableList See [lavaanList-class]
##' @slot vcovList See [lavaanList-class]
##' @slot testList See [lavaanList-class]
##' @slot h1List See [lavaanList-class]. An additional element is
##'   added to the `list`: `$PT` is the "saturated" model's parameter
##'   table, returned by [lavaan::lav_partable_unrestricted()].
##' @slot baselineList See [lavaanList-class]
##'
##' @param object An object of class [lavaan.mi-class]
##' @param header,fit.measures,fm.args,estimates,ci,standardized,std,cov.std,rsquare,remove.unused,modindices,nd,output
##'        See descriptions of `summary()` arguments in the help page for
##'        [lavaan-class] class. Also see [lavaan::fitMeasures()] for arguments
##'        `fit.measures` and `fm.args`.
##' @param baseline.model,h1.model See [lavaan::fitMeasures()].
##' @param ... Additional arguments passed to [lavTestLRT.mi()], or
##'        subsequently to [lavaan::lavTestLRT()].
##' @param fmi `logical` indicating whether to add the Fraction Missing
##'        Information (FMI) and (average) relative increase in variance (ARIV)
##'        to the output.
##' @param asymptotic `logical`. If `FALSE` (typically a default, but
##'        see **Value** section for details using various methods), pooled
##'        tests (of fit or pooled estimates) will be *F* or *t*
##'        statistics with associated degrees of freedom (*df*). If
##'        `TRUE`, the (denominator) *df* are assumed to be
##'        sufficiently large for a *t* statistic to follow a normal
##'        distribution, so it is printed as a *z* statistic; likewise,
##'        *F* times its numerator *df* is printed, assumed to follow
##'        a \eqn{\chi^2} distribution.
##' @param scale.W `logical`. If `TRUE` (default), the `vcov`
##'        method will calculate the pooled covariance matrix by scaling the
##'        within-imputation component by the ARIV (see Enders, 2010, p. 235,
##'        for definition and formula). Otherwise, the pooled matrix is
##'        calculated as the weighted sum of the within-imputation and
##'        between-imputation components (see Enders, 2010, ch. 8, for details).
##'        This in turn affects how the `summary` method calculates its
##'        pooled standard errors, as well as the Wald test
##'        ([lavTestWald.mi()]).
##' @param omit.imps `character` vector specifying criteria for omitting
##'        imputations from pooled results.  Can include any of
##'        `c("no.conv", "no.se", "no.npd")`, the first 2 of which are the
##'        default setting, which excludes any imputations that did not
##'        converge or for which standard errors could not be computed.  The
##'        last option (`"no.npd"`) would exclude any imputations which
##'        yielded a nonpositive definite covariance matrix for observed or
##'        latent variables, which would include any "improper solutions" such
##'        as Heywood cases.  NPD solutions are not excluded by default because
##'        they are likely to occur due to sampling error, especially in small
##'        samples.  However, gross model misspecification could also cause
##'        NPD solutions, users can compare pooled results with and without
##'        this setting as a sensitivity analysis to see whether some
##'        imputations warrant further investigation. Specific imputation
##'        numbers can also be included in this argument, in case users want to
##'        apply their own custom omission criteria (or simulations can use
##'        different numbers of imputations without redundantly refitting the
##'        model).
##' @param labels `logical` indicating whether the `coef()` output
##'        should include parameter labels. Default is `TRUE`.
##' @param total `logical` (default: `TRUE`) indicating whether the
##'        `nobs()` method should return the total sample size or (if
##'        `FALSE`) a vector of group sample sizes.
##' @param type The meaning of this argument varies depending on which method it
##'        it used for. Find detailed descriptions in the **Value** section
##'        under `coef()` and `vcov()`.
##'
##' @return
##'
##' \item{coef}{`signature(object = "lavaan.mi", type = "free",
##'   labels = TRUE, omit.imps = c("no.conv","no.se"))`:
##'   See argument description on the help page for [lavaan-class] class.
##'   Returns the pooled point estimates (i.e., averaged across imputed data
##'   sets; see Rubin, 1987).}
##'
##' \item{vcov}{`signature(object = "lavaan.mi", scale.W = TRUE,
##'   omit.imps = c("no.conv","no.se"),
##'   type = c("pooled","between","within","ariv"))`:  By default, returns the
##'   pooled covariance matrix of parameter estimates (`type = "pooled"`),
##'   the within-imputations covariance matrix (`type = "within"`), the
##'   between-imputations covariance matrix (`type = "between"`), or the
##'   average relative increase in variance (`type = "ariv"`) due to
##'   missing data.}
##'
##' \item{fitted.values}{`signature(object = "lavaan.mi",
##'   omit.imps = c("no.conv","no.se"))`: See corresponding [lavaan-class] method.
##'   Returns model-implied moments, evaluated at the pooled point estimates.}
##' \item{fitted}{alias for `fitted.values`}
##'
##' \item{nobs}{`signature(object = "lavaan.mi", total = TRUE)`: either
##'   the total (default) sample size or a vector of group sample sizes
##'   (`total = FALSE`).}
##'
##' \item{fitMeasures}{`signature(object = "lavaan.mi",
##'     fit.measures = "all", baseline.model = NULL, h1.model = NULL,
##'     fm.args = list(standard.test = "default", scaled.test = "default",
##'     rmsea.ci.level = 0.90, rmsea.h0.closefit = 0.05,
##'     rmsea.h0.notclosefit = 0.08, robust = TRUE, cat.check.pd = TRUE),
##'     output = "vector", omit.imps = c("no.conv","no.se"), ...)`:
##'   See [lavaan::fitMeasures()] for details.
##'   Pass additional arguments to [lavTestLRT.mi()] via `...`.}
##' \item{fitmeasures}{alias for `fitMeasures`.}
##'
##' \item{show}{`signature(object = "lavaan.mi")`: returns a message about
##'   convergence rates and estimation problems (if applicable) across imputed
##'   data sets.}
##'
##' \item{summary}{`signature(object = "lavaan.mi", header = TRUE,
##'    fit.measures = FALSE,fm.args = list(standard.test = "default",
##'    scaled.test = "default", rmsea.ci.level = 0.90, rmsea.h0.closefit = 0.05,
##'    rmsea.h0.notclosefit = 0.08, robust = TRUE, cat.check.pd = TRUE),
##'    estimates = TRUE, ci = FALSE, standardized = FALSE, std = standardized,
##'    cov.std = TRUE, rsquare = FALSE, fmi = FALSE, asymptotic = FALSE,
##'    scale.W = !asymptotic, omit.imps = c("no.conv","no.se"),
##'    remove.unused = TRUE, modindices = FALSE, nd = 3L, ...)`:
##'  Analogous to `summary()` for `lavaan-class` objects.
##'  By default, `summary` returns output from [parameterEstimates.mi()],
##'  with some cursory information in the header.
##'  Setting `fit.measures=TRUE` will additionally run `fitMeasures()`, and
##'  setting `modindices=TRUE` will additionally run [modificationIndices.mi()].}
##'
##' @section Objects from the Class: See the [lavaan.mi()] function
##'   for details. Wrapper functions include [cfa.mi()],
##'   [sem.mi()], and [growth.mi()].
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Enders, C. K. (2010). *Applied missing data analysis*. New York, NY:
##'   Guilford.
##'
##'   Rubin, D. B. (1987). *Multiple imputation for nonresponse in surveys*.
##'   New York, NY: Wiley. \doi{10.1002/9780470316696}
##'
##'
##' @examples
##'
##' ## See ?lavaan.mi help page
##'
##' @export
setClass("lavaan.mi", contains = "lavaanList",
         slots = c(coefList = "list",     # coefficients in matrix format
                   phiList = "list",      # list of model-implied latent covariance matrices
                   miList = "list",       # modification indices
                   lavListCall = "list",  # store actual call to lavaanList
                   convergence = "list")) # also check SEs and Heywood cases


## show() basic info at the top of summary()
lavaan_mi_short_summary <- function(object, return.string = FALSE) {
  nData <- object@meta$ndat

  useImps <- sapply(object@convergence, "[[", i = "converged")
  nConverged <- sum(useImps)

  SE <- sapply(object@convergence, "[[", "SE")
  SE[is.na(SE)] <- FALSE

  Heywood.ov <- sapply(object@convergence, "[[", "Heywood.ov")
  Heywood.ov[is.na(Heywood.ov)] <- FALSE

  Heywood.lv <- sapply(object@convergence, "[[", "Heywood.lv")
  Heywood.lv[is.na(Heywood.lv)] <- FALSE

  ## assemble a message to print
  MESSAGE <- paste0('lavaan.mi object fit to ', nData,
                    ' imputed data sets using:\n',
                    ' - lavaan    (', object@version[1],
                    ')\n - lavaan.mi (', object@version[2],')\n',
                    'See class?lavaan.mi help page for available methods. \n\n',
                    'Convergence information:\n', 'The model converged on ',
                    nConverged, ' imputed data sets.\n')

  if (all(SE)) {
    MESSAGE <- c(MESSAGE,
                 'Standard errors were available for all imputations.\n')
  } else {
    MESSAGE <- c(MESSAGE,
                 paste('Standard errors could not be computed for data set(s)',
                       paste(which(!SE), collapse = ", "), '\nTry fitting the',
                       'model to the individual data set(s) to diagnose',
                       'problems. If they cannot be fixed, try inspecting the',
                       'imputations. It may be necessary to reimpute the data',
                       'with some restrictions imposed.\n'))
  }

  if (any(Heywood.ov | Heywood.lv)) {
    MESSAGE <- c(MESSAGE,
                 paste('Heywood cases detected for data set(s)',
                       paste(which(Heywood.ov | Heywood.lv), collapse = ", "),
                       '\nThese are not necessarily a cause for concern, unless',
                       'a pooled estimate is also a Heywood case.\n\n'))
  }

  #TODO: add test stat? (like lavaan-class prints)

  if (return.string) return(paste(MESSAGE, collapse = ""))
  ## else
  object
}
##' @name lavaan.mi-class
##' @aliases show,lavaan.mi-method
##' @importFrom methods show
##' @export
setMethod("show", "lavaan.mi", function(object) {
  MESSAGE <- lavaan_mi_short_summary(object, return.string = TRUE)
  cat(MESSAGE)
  object
})



## analog to lavaan:::lav_object_summary(), which creates a list of
## components to appear in the summary() output.  Creating a similar
## object allows lavaan.mi to capitalize on lavaan:::print.lavaan.summary()
lavaan_mi_object_summary <- function(object, omit.imps = c("no.conv", "no.se"),
                                     asymptotic = FALSE, scale.W = !asymptotic,
                                     #TODO: add pool.method=
                                     header = TRUE, fit.measures = FALSE,
                                     fm.args = list(standard.test   = "default",
                                                    scaled.test     = "default",
                                                    rmsea.ci.level       = 0.90,
                                                    rmsea.h0.closefit    = 0.05,
                                                    rmsea.h0.notclosefit = 0.08,
                                                    robust               = TRUE,
                                                    cat.check.pd         = TRUE),
                                     estimates = TRUE,
                                     ci = FALSE,
                                     fmi = FALSE,
                                     standardized = FALSE,
                                     std = standardized,
                                     remove.unused = TRUE,
                                     cov.std = TRUE,
                                     rsquare = FALSE,
                                     modindices = FALSE) {
  #TODO: Keep a close eye on changes in lavaan:::lav_object_summary().
  #      The code below is basically copy/pasted, with
  #      a few changes appropriate for lavaanList objects.

  # return a list with the main ingredients
  res <- list()

  # this is to avoid partial matching of 'std' with std.nox
  if (is.logical(std) && is.logical(standardized)) {
    standardized <- std || standardized
  } else {
    # At least 1 is not logical. Retain only valid standardization options.
    standardized <- intersect(union(tolower(std), tolower(standardized)),
                              c("std.lv","std.all","std.nox"))
  }

  if (header) {
    ## custom "top" for lavaan.mi objects:
    res$top_of_lavaanmi <- lavaan_mi_short_summary(object, return.string = TRUE)
    ## the rest below is copied from lavaan:::lav_object_summary()

    # 2. summarize optim info (including estimator)
    res$optim <- list(
      estimator = object@Options$estimator,
      estimator.args = object@Options$estimator.args,
      optim.method = object@Options$optim.method,
      npar = object@Model@nx.free,
      eq.constraints = object@Model@eq.constraints,
      nrow.ceq.jac = nrow(object@Model@ceq.JAC),
      nrow.cin.jac = nrow(object@Model@cin.JAC),
      nrow.con.jac = nrow(object@Model@con.jac),
      con.jac.rank = qr(object@Model@con.jac)$rank
    )

    # 4. summarize lavdata (copy lavaan's internal function)
    ##' @importFrom methods .hasSlot
    lav_data_summary_short <- function(lavdata) {

      # two or three columns (depends on nobs/norig)
      threecolumn <- FALSE
      for (g in 1:lavdata@ngroups) {
        if (lavdata@nobs[[g]] != lavdata@norig[[g]]) {
          threecolumn <- TRUE
          break
        }
      }

      # clustered data?
      clustered <- FALSE
      if (.hasSlot(lavdata, "cluster") && # in case we have an old obj
          length(lavdata@cluster) > 0L) {
        clustered <- TRUE
      }

      # multilevel data?
      multilevel <- FALSE
      if (.hasSlot(lavdata, "nlevels") && # in case we have an old obj
          lavdata@nlevels > 1L) {
        multilevel <- TRUE
      }

      # extract summary information
      datasummary <- list(
        ngroups = lavdata@ngroups,
        nobs = unlist(lavdata@nobs)
      )

      # norig?
      if (threecolumn) {
        datasummary$norig <- unlist(lavdata@norig)
      }

      # multiple groups?
      if (lavdata@ngroups > 1L) {
        datasummary$group.label <- lavdata@group.label
      }

      # sampling weights?
      if ((.hasSlot(lavdata, "weights")) && # in case we have an old object
          (!is.null(lavdata@weights[[1L]]))) {
        datasummary$sampling.weights <- lavdata@sampling.weights
      }

      # clustered/multilevel data?
      if (clustered) {
        if (multilevel) {
          datasummary$nlevels <- lavdata@nlevels
        }
        datasummary$cluster <- lavdata@cluster

        if (lavdata@ngroups == 1L) {
          datasummary$nclusters <- unlist(lavdata@Lp[[1]]$nclusters)
        } else {
          tmp <- vector("list", length = lavdata@ngroups)
          for (g in seq_len(lavdata@ngroups)) {
            tmp[[g]] <- unlist(lavdata@Lp[[g]]$nclusters)
          }
          datasummary$nclusters <- tmp
        }
      }

      # missing data?
      if (!is.null(lavdata@Mp[[1L]])) {
        datasummary$npatterns <- sapply(lavdata@Mp, "[[", "npatterns")
        if (multilevel && !is.null(lavdata@Mp[[1L]]$Zp)) {
          datasummary$npatterns2 <- sapply(lapply(
            lavdata@Mp,
            "[[", "Zp"
          ), "[[", "npatterns")
        }
      }

      datasummary
    }

    res$data <- lav_data_summary_short(lavdata = object@Data)

    ## lavaan prints the test here, but pooling takes time.
    ## Instead, conditionally add pooled test with(out) fit.measures
  }
  if (fit.measures) {
    res$fit <- mi_fit_indices_via_lavaan(object, fm.args = fm.args,
                                         omit.imps = omit.imps,
                                         fit.measures = "default")
  } else {
    ## just the pooled test(s)
    chiSqTests <- c("chisq","df","pvalue",   "chisq.scaling.factor",
                    "chisq.scaled","df.scaled","pvalue.scaled")
    res$fit <- mi_fit_indices_via_lavaan(object, fm.args = fm.args,
                                         fit.measures = chiSqTests)
    attr(res$fit, "add.h0") <- TRUE
  }

  if (estimates) {
    PE <- parameterEstimates.mi(object, omit.imps = omit.imps,
                                asymptotic = asymptotic, scale.W = scale.W,
                                ci = ci, standardized = standardized,
                                rsquare = rsquare, fmi = fmi, cov.std = cov.std,
                                remove.eq = FALSE, remove.system.eq = TRUE,
                                remove.ineq = FALSE, remove.def = FALSE,
                                remove.nonfree = FALSE,
                                remove.unused = remove.unused,
                                output = "text", header = TRUE)
    res$pe <- as.data.frame(PE)
  }

  if (modindices) {
    MI <- modificationIndices.mi(object, omit.imps = omit.imps,
                                 #TODO: add pool.method=
                                 standardized = TRUE, cov.std = cov.std)
    res$mi <- MI
  }

  ## return lavaan.summary S3-class object for lavaan's print method
  class(res) <- c("lavaan.summary", "list")
  res
}
##' @importFrom stats pt qt pnorm qnorm
##' @importFrom lavaan lavListInspect parTable lavNames
##' @importFrom methods getMethod
summary_lavaan_mi <- function(object, header = TRUE,
                              fit.measures = FALSE,
                              fm.args = list(standard.test        = "default",
                                             scaled.test          = "default",
                                             rmsea.ci.level       = 0.90,
                                             rmsea.h0.closefit    = 0.05,
                                             rmsea.h0.notclosefit = 0.08,
                                             robust               = TRUE,
                                             cat.check.pd         = TRUE),
                              estimates = TRUE,
                              ci = FALSE, #level = .95,
                              ## standardization
                              standardized = FALSE, std = standardized,
                              cov.std = TRUE, rsquare = FALSE,
                              ## control over pooling
                              fmi = FALSE, asymptotic = FALSE,
                              scale.W = !asymptotic,
                              omit.imps = c("no.conv","no.se"),
                              ## remove rows?
                              remove.unused = TRUE,
                              modindices = FALSE, nd = 3L,
                              ...) {

  SUM <- lavaan_mi_object_summary(object = object, omit.imps = omit.imps,
                                  asymptotic = asymptotic, scale.W = scale.W,
                                  header = header,
                                  fit.measures = fit.measures, fm.args = fm.args,
                                  estimates = estimates, ci = ci, fmi = fmi,
                                  std = std, standardized = standardized,
                                  cov.std = cov.std, rsquare = rsquare,
                                  remove.unused = remove.unused,
                                  modindices = modindices)
  attr(SUM, "nd") <- nd # save as attribute for lavaan's print method
  SUM
}
##' @name lavaan.mi-class
##' @aliases summary,lavaan.mi-method
##' @export
setMethod("summary", "lavaan.mi", summary_lavaan_mi)


##' @name lavaan.mi-class
##' @aliases nobs,lavaan.mi-method
##' @importFrom stats nobs
##' @importFrom lavaan lavListInspect
##' @export
setMethod("nobs", "lavaan.mi", function(object, total = TRUE) {
  if (total) return(lavListInspect(object, "ntotal"))
  #FIXME: cluster N for multilevel?
  N <- lavListInspect(object, "norig")
  if (length(N) > 1L) names(N) <- lavListInspect(object, "group.label")
  N
})



##' @importFrom stats coef
##' @importFrom lavaan parTable
coef_lavaan_mi <- function(object, type = "free", labels = TRUE,
                           omit.imps = c("no.conv","no.se")) {
  useImps <- imps2use(object = object, omit.imps = omit.imps)

  PT <- parTable(object)
  if (type == "user" || type == "all") {
    type <- "user"
    idx <- 1:length(PT$lhs)
  } else if (type == "free") {
    ## FIXME: duplicated leftover from old way of handling EQ constraints?
    idx <- which(PT$free > 0L & !duplicated(PT$free))
  }
  ## extract coefficients for converged models
  coefList <- lapply(object@ParTableList[useImps], "[[", i = "est")
  out <- colMeans(do.call(rbind, coefList))[idx]
  ## attach names, set class
  if (labels) names(out) <- lavaan::lav_partable_labels(PT, type = type)
  class(out) <- c("lavaan.vector","numeric")
  out
}
##' @name lavaan.mi-class
##' @aliases coef,lavaan.mi-method
##' @export
setMethod("coef", "lavaan.mi", coef_lavaan_mi)



##' @importFrom stats cov vcov
##' @importFrom lavaan lavListInspect parTable
vcov_lavaan_mi <- function(object, type = c("pooled","between","within","ariv"),
                           scale.W = TRUE, omit.imps = c("no.conv","no.se")) {
  useImps <- imps2use(object = object, omit.imps = omit.imps)
  m <- length(useImps)

  if (lavListInspect(object, "options")$se == "none") {
    warning('requested se="none", so only between-imputation (co)variance can',
            ' be computed')
    type <- "between"
  }
  type <- tolower(type[1])
  if (!(type %in% c("pooled","between","within","ariv")))
    stop("'", type, "' is not a valid option for 'type'")

  PT <- parTable(object)
  ncon <- sum(PT$op == "==")
  npar <- max(PT$free) - ncon

  coefList <- lapply(object@ParTableList[useImps], "[[", i = "est")
  B <- cov(do.call(rbind, coefList)[ , PT$free > 0L & !duplicated(PT$free)])
  class(B) <- c("lavaan.matrix.symmetric","matrix")
  rownames(B) <- colnames(B) <- lavaan::lav_partable_labels(PT, type = "free")
  if (type == "between") return(B)

  W <- Reduce("+", lapply(object@vcovList[useImps], function(x) x$vcov)) / m
  class(W) <- c("lavaan.matrix.symmetric","matrix")
  dimnames(W) <- dimnames(B)
  if (type == "within") return(W)

  ## check whether equality constraints prevent inversion of W
  if (scale.W || type == "ariv") {
    inv.W <- if (ncon == 0) try(solve(W), silent = TRUE) else MASS::ginv(W)
    if (inherits(inv.W, "try-error")) {
      if (ncon == 0) {
        warning("Could not invert within-imputation covariance matrix. ",
                "Generalized inverse used instead.\nIt may be ",
                "safer to set `scale.W = FALSE' (and `asymptotic = TRUE').")
      }
      inv.W <- MASS::ginv(W)
    }
    ## relative increase in variance due to missing data
    r <- (1 + 1/m)/npar * sum(diag(B %*% inv.W)) # Enders (2010, p. 235) eqs. 8.20-21
    if (type == "ariv") return(r)
    Total <- (1 + r) * W # FIXME: asked Yves for a hack, says it can't be inverted back to infoMat
  } else {
    ## less reliable, but constraints prevent inversion of W
    Total <- W + B + (1/m)*B ## Enders (2010, p. 235) eq. 8.19
  }
  ## return pooled variance
  Total
}
##' @name lavaan.mi-class
##' @aliases vcov,lavaan.mi-method
##' @export
setMethod("vcov", "lavaan.mi", vcov_lavaan_mi)



##' @importFrom lavaan lavListInspect lavNames
##' @importFrom stats fitted fitted.values
fitted_lavaan_mi <- function(object, momentsNblocks = TRUE, # the way users see it
                             ## momentsNblocks = FALSE is how lavaan stores it
                             omit.imps = c("no.conv","no.se")) {
  useImps <- imps2use(object = object, omit.imps = omit.imps)

  ## how many blocks to loop over
  nG <- lavListInspect(object, "ngroups")
  nlevels <- lavListInspect(object, "nlevels")
  nBlocks <- nG * nlevels #FIXME: always?
  group.label <- if (nG > 1L) lavListInspect(object, "group.label") else NULL
  clus.label <- if (nlevels > 1L) c("within", lavListInspect(object, "cluster")) else NULL
  if (nBlocks > 1L) {
    block.label <- paste(rep(group.label, each = nlevels), clus.label,
                         sep = if (nG > 1L && nlevels > 1L) "_" else "")
  }

  est <- coef_lavaan_mi(object, omit.imps = omit.imps)
  setpar <- lavaan::lav_model_set_parameters(object@Model, x = est)
  impMats <- lavaan::lav_model_implied(setpar)
  # if (lavListInspect(object, "categorical")) {
  #   th.idx <- lavListInspect(object, "th.idx") # to select $(res.)th
  #   if (nBlocks == 1L) th.idx <- list(th.idx)  # to loop over
      #FIXME when multilevel accepts categorical
  # }

  ## blocks nested in moments, for use in mi2lavaan(rmr=TRUE)
  if (!momentsNblocks) return(impMats)


  #TODO: adapt to multilevel, multigroup, or both
  ## loop over (blocks and) moments
  Implied <- vector("list", nBlocks)
  for (b in 1:nBlocks) {
    for (nm in names(impMats)) {

      ## skip any empty objects
      if (is.null(impMats[[nm]][[b]])) next

      Implied[[b]][[nm]] <- impMats[[nm]][[b]]

      ## assign names and classes
      if (nm %in% c("cov","res.cov")) {
        NAMES <- lavNames(object, type = "ov.model", block = b)
        dimnames(Implied[[b]][[nm]]) <- list(NAMES, NAMES)
        class(Implied[[b]][[nm]]) <- c("lavaan.matrix.symmetric","matrix")

      } else if (nm %in% c("mean","res.int")) {
        Implied[[b]][[nm]] <- as.numeric(Implied[[b]][[nm]]) # remove matrix
        names(Implied[[b]][[nm]]) <- lavNames(object, type = "ov.model", block = b)
        class(Implied[[b]][[nm]]) <- c("lavaan.vector","numeric")

      } else if (nm %in% c("th","res.th")) {
        #FIXME: When lavaan allows multilevel categorical, thresholds only
        ##      apply once (not to each level, like for all groups).
        ##      Will lavaan return a vector of zeros for all but "within"?
        ##      If not, it will not exist for each block, so count over groups.
        Implied[[b]][[nm]] <- as.numeric(Implied[[b]][[nm]]) #[ th.idx[[b]] ] # remove matrix & numeric -means
        names(Implied[[b]][[nm]]) <- lavNames(object, type = "th",
                                              block = b) #FIXME?
        class(Implied[[b]][[nm]]) <- c("lavaan.vector","numeric")

      } else if (nm == "group.w") {
        ## Only for (D)WLS estimation, but when is it relevant?
        ## For now, assign no names/class


      ## The remaining only exist when conditional.x
      } else if (nm %in% c("slopes","res.slopes")) {
        dimnames(Implied[[b]][[nm]]) <- list(lavNames(object, type = "ov.nox", block = b),
                                             lavNames(object, type = "ov.x", block = b))
        class(Implied[[b]][[nm]]) <- c("lavaan.matrix","matrix")

      } else if (nm == "cov.x") {
        NAMES <- lavNames(object, type = "ov.x", block = b)
        dimnames(Implied[[b]][[nm]]) <- list(NAMES, NAMES)
        class(Implied[[b]][[nm]]) <- c("lavaan.matrix.symmetric","matrix")

      } else if (nm == "mean.x") {
        Implied[[b]][[nm]] <- as.numeric(Implied[[b]][[nm]]) # remove matrix
        names(Implied[[b]][[nm]]) <- lavNames(object, type = "ov.x", block = b)
        class(Implied[[b]][[nm]]) <- c("lavaan.vector","numeric")
      }

    ## end loops
    }
  }

  ## drop list for 1 block, or add labels for multiple
  if (nBlocks == 1L) {
    Implied <- Implied[[1]]
  } else names(Implied) <- block.label

  Implied
}
##' @name lavaan.mi-class
##' @aliases fitted,lavaan.mi-method
##' @export
setMethod("fitted", "lavaan.mi",
          function(object, omit.imps = c("no.conv","no.se")) {
  fitted_lavaan_mi(object, momentsNblocks = TRUE, omit.imps = omit.imps)
})
##' @name lavaan.mi-class
##' @aliases fitted.values,lavaan.mi-method
##' @export
setMethod("fitted.values", "lavaan.mi",
          function(object, omit.imps = c("no.conv","no.se")) {
  fitted_lavaan_mi(object, momentsNblocks = TRUE, omit.imps = omit.imps)
})


## utility function called within mi2lavaan()
## and formerly within resid_lavaan_mi()
pool_h1 <- function(object, momentsNblocks = TRUE, # the way users see it
                    ## momentsNblocks = FALSE is how lavaan stores it
                    omit.imps = c("no.conv","no.se")) {
  useImps <- imps2use(object = object, omit.imps = omit.imps)
  m <- length(useImps)

  nBlocks <- object@Model@nblocks
  momentNames <- names(object@h1List[[ useImps[1] ]]$implied)

  if (lavListInspect(object, "categorical")) {
    th.idx <- lavListInspect(object, "th.idx") # to select $(res.)th
    if (nBlocks == 1L) th.idx <- list(th.idx)  # to loop over
    #FIXME when multilevel accepts categorical
  }

  ## template to store saturated moments
  OBS <- vector("list", ifelse(momentsNblocks, nBlocks, length(momentNames)))

  ## loop over (blocks and) moments
  for (b in 1:nBlocks) {
    for (nm in momentNames) {

      ## skip if Implied element is not part of the saturated list
      if (is.null(object@h1List[[ useImps[1] ]]$implied[[nm]][[b]])) next

      ## H1 (saturated model) implied moments
      ## (block-list nested in moments-list)
      momentList <- lapply(object@h1List[useImps],
                           function(x) x$implied[[nm]][[b]])
      target <- Reduce("+", momentList) / m
      #TODO: unnecessary calculation if standardized and nm %in% c("th","slopes")

      ## only for fitted() method
      if (momentsNblocks  &&  nm %in% c("th","res.th")) {
        ## remove numeric -means from thresholds
        target <- as.numeric(target)[ th.idx[[b]] ]
      }

      if (momentsNblocks) {
        OBS[[b]][[nm]] <- target # formerly used by (old_)resid_lavaan_mi()
      } else {
        OBS[[nm]][[b]] <- target # used by mi2lavaan()
      }
      ## end loop over moments
    }
    ## end loop over blocks
  }

  OBS
}



