### Terrence D. Jorgensen
### Last updated: 1 April 2024
### Class and Methods for lavaan.mi object


##' Class for a lavaan Model Fitted to Multiple Imputations
##'
##' This class extends the \code{\linkS4class{lavaanList}} class, created by
##' fitting a lavaan model to a list of data sets. In this case, the list of
##' data sets are multiple imputations of missing data.
##'
##'
##' @name lavaan.mi-class
##' @importClassesFrom lavaan lavaanList
##' @aliases lavaan.mi-class show,lavaan.mi-method summary,lavaan.mi-method
##'   fitMeasures,lavaan.mi-method fitmeasures,lavaan.mi-method
##'   anova,lavaan.mi-method nobs,lavaan.mi-method coef,lavaan.mi-method
##'   vcov,lavaan.mi-method fitted,lavaan.mi-method fitted.values,lavaan.mi-method
##'   residuals,lavaan.mi-method resid,lavaan.mi-method
##' @docType class
##'
##' @slot coefList \code{list} of estimated coefficients in matrix format (one
##'   per imputation) as output by \code{\link[lavaan]{lavInspect}(fit, "est")}
##' @slot phiList \code{list} of model-implied latent-variable covariance
##'   matrices (one per imputation) as output by
##'   \code{\link[lavaan]{lavInspect}(fit, "cov.lv")}
##' @slot miList \code{list} of modification indices output by
##'   \code{\link[lavaan]{modindices}}
##' @slot lavListCall call to \code{\link[lavaan]{lavaanList}} used to fit the
##'   model to the list of imputed data sets in \code{@@DataList}, stored as a
##'   \code{list} of arguments
##' @slot convergence \code{list} of \code{logical} vectors indicating whether,
##'   for each imputed data set, (1) the model converged on a solution, (2)
##'   \emph{SE}s could be calculated, (3) the (residual) covariance matrix of
##'   latent variables (\eqn{\Psi}) is non-positive-definite, and (4) the
##'   residual covariance matrix of observed variables (\eqn{\Theta}) is
##'   non-positive-definite.
##' @slot lavaanList_slots All remaining slots are from
##'   \code{\linkS4class{lavaanList}}, but \code{\link{lavaan.mi}} only populates a
##'   subset of the \code{list} slots, two of them with custom information:
##' @slot DataList The \code{list} of imputed data sets
##' @slot SampleStatsList List of output from
##'   \code{\link[lavaan]{lavInspect}(fit, "sampstat")} applied to each fitted
##'   model
##' @slot ParTableList See \code{\linkS4class{lavaanList}}
##' @slot vcovList See \code{\linkS4class{lavaanList}}
##' @slot testList See \code{\linkS4class{lavaanList}}
##' @slot h1List See \code{\linkS4class{lavaanList}}. An additional element is
##'   added to the \code{list}: \code{$PT} is the "saturated" model's parameter
##'   table, returned by \code{\link[lavaan]{lav_partable_unrestricted}}.
##' @slot baselineList See \code{\linkS4class{lavaanList}}
##'
##' @param object An object of class \code{lavaan.mi}
##' @param se,ci,level,standardized,rsquare,header,output See
##'        \code{\link[lavaan]{parameterEstimates}}. \code{output}
##'        can also be passed to \code{\link[lavaan]{fitMeasures}}.
##' @param fmi \code{logical} indicating whether to include the Fraction Missing
##'        Information (FMI) for parameter estimates in the \code{summary}
##'        output (see \bold{Value} section).
##' @param asymptotic \code{logical}. If \code{FALSE} (typically a default, but
##'        see \bold{Value} section for details using various methods), pooled
##'        tests (of fit or pooled estimates) will be \emph{F} or \emph{t}
##'        statistics with associated degrees of freedom (\emph{df}). If
##'        \code{TRUE}, the (denominator) \emph{df} are assumed to be
##'        sufficiently large for a \emph{t} statistic to follow a normal
##'        distribution, so it is printed as a \emph{z} statistic; likewise,
##'        \emph{F} times its numerator \emph{df} is printed, assumed to follow
##'        a \eqn{\chi^2} distribution.
##' @param scale.W \code{logical}. If \code{TRUE} (default), the \code{vcov}
##'        method will calculate the pooled covariance matrix by scaling the
##'        within-imputation component by the ARIV (see Enders, 2010, p. 235,
##'        for definition and formula). Otherwise, the pooled matrix is
##'        calculated as the weighted sum of the within-imputation and
##'        between-imputation components (see Enders, 2010, ch. 8, for details).
##'        This in turn affects how the \code{summary} method calculates its
##'        pooled standard errors, as well as the Wald test
##'        (\code{\link{lavTestWald.mi}}).
##' @param omit.imps \code{character} vector specifying criteria for omitting
##'        imputations from pooled results.  Can include any of
##'        \code{c("no.conv", "no.se", "no.npd")}, the first 2 of which are the
##'        default setting, which excludes any imputations that did not
##'        converge or for which standard errors could not be computed.  The
##'        last option (\code{"no.npd"}) would exclude any imputations which
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
##' @param labels \code{logical} indicating whether the \code{coef} output
##'        should include parameter labels. Default is \code{TRUE}.
##' @param total \code{logical} (default: \code{TRUE}) indicating whether the
##'        \code{nobs} method should return the total sample size or (if
##'        \code{FALSE}) a vector of group sample sizes.
##' @param type The meaning of this argument varies depending on which method it
##'        it used for. Find detailed descriptions in the \bold{Value} section
##'        under \code{coef}, \code{vcov}, and \code{residuals}.
##' @param fit.measures,baseline.model,fm.args See \code{\link[lavaan]{fitMeasures}}.
##'        \code{summary(object, fit.measures = TRUE)} will print (but not
##'        return) a table of fit measures to the console.
##' @param ... Additional arguments passed to \code{\link{lavTestLRT.mi}}, or
##'        subsequently to \code{\link[lavaan]{lavTestLRT}}.
##'
##' @return
##'
##' \item{coef}{\code{signature(object = "lavaan.mi", type = "free",
##'   labels = TRUE, omit.imps = c("no.conv","no.se"))}:
##'   See \code{\linkS4class{lavaan}}. Returns the pooled point estimates (i.e.,
##'   averaged across imputed data sets; see Rubin, 1987).}
##'
##' \item{vcov}{\code{signature(object = "lavaan.mi", scale.W = TRUE,
##'   omit.imps = c("no.conv","no.se"),
##'   type = c("pooled","between","within","ariv"))}:  By default, returns the
##'   pooled covariance matrix of parameter estimates (\code{type = "pooled"}),
##'   the within-imputations covariance matrix (\code{type = "within"}), the
##'   between-imputations covariance matrix (\code{type = "between"}), or the
##'   average relative increase in variance (\code{type = "ariv"}) due to
##'   missing data.}
##'
##' \item{fitted.values}{\code{signature(object = "lavaan.mi",
##'   omit.imps = c("no.conv","no.se"))}: See \code{\linkS4class{lavaan}}.
##'   Returns model-implied moments, evaluated at the pooled point estimates.}
##' \item{fitted}{alias for \code{fitted.values}}
##'
##' \item{residuals}{\code{signature(object = "lavaan.mi",
##'   type = c("raw","cor"), omit.imps = c("no.conv","no.se"))}:
##'   See \code{\linkS4class{lavaan}}. By default (\code{type = "raw"}), returns
##'   the difference between the model-implied moments from \code{fitted.values}
##'   and the pooled observed moments (i.e., averaged across imputed data sets).
##'   Standardized residuals are also available, using Bollen's
##'   (\code{type = "cor"} or \code{"cor.bollen"}) or Bentler's
##'   (\code{type = "cor.bentler"}) formulas.}
##' \item{resid}{alias for \code{residuals}}
##'
##' \item{nobs}{\code{signature(object = "lavaan.mi", total = TRUE)}: either
##'   the total (default) sample size or a vector of group sample sizes
##'   (\code{total = FALSE}).}
##'
##' \item{anova}{\code{signature(object = "lavaan.mi", ...)}:
##'   Returns a test of model fit for a single model (\code{object}) or test(s)
##'   of the difference(s) in fit between nested models passed via \code{...}.
##'   This is just a wrapper around \code{\link{lavTestLRT.mi}}, where you can
##'   find details about additional arguments.}
##'
##' \item{fitMeasures}{\code{signature(object = "lavaan.mi",
##'   fit.measures = "all", baseline.model = NULL,
##'   fm.args = list(standard.test = "default", scaled.test = "default",
##'   rmsea.ci.level = 0.90, rmsea.close.h0 = 0.05, rmsea.notclose.h0 = 0.08,
##'   cat.check.pd = TRUE), output = "vector", omit.imps = c("no.conv","no.se"),
##'   ...)}: See lavaan's  \code{\link[lavaan]{fitMeasures}} for details.
##'   Pass additional arguments to \code{\link{lavTestLRT.mi}} via \code{...}.}
##' \item{fitmeasures}{alias for \code{fitMeasures}.}
##'
##' \item{show}{\code{signature(object = "lavaan.mi")}: returns a message about
##'   convergence rates and estimation problems (if applicable) across imputed
##'   data sets.}
##'
##' \item{summary}{\code{signature(object = "lavaan.mi", se = TRUE, ci = FALSE,
##'   level = .95, standardized = FALSE, rsquare = FALSE, fmi = FALSE,
##'   scale.W = !asymptotic, omit.imps = c("no.conv","no.se"),
##'   asymptotic = FALSE, header = TRUE, output = "text", fit.measures = FALSE,
##'   fm.args = list(standard.test = "default", scaled.test = "default",
##'                  rmsea.ci.level = 0.90, rmsea.h0.closefit = 0.05,
##'                  rmsea.h0.notclosefit = 0.08), ...)}:
##'  see \code{\link[lavaan]{parameterEstimates}} for details.
##'  By default, \code{summary} returns pooled point and \emph{SE}
##'  estimates, along with \emph{t} test statistics and their associated
##'  \emph{df} and \emph{p} values. If \code{ci = TRUE}, confidence intervals
##'  are returned with the specified confidence \code{level} (default 95\% CI).
##'  If \code{asymptotic = TRUE}, \emph{z} instead of \emph{t} tests are
##'  returned. \code{standardized} solution(s) can also be requested by name
##'  (\code{"std.lv"} or \code{"std.all"}) or both are returned with \code{TRUE}.
##'  \emph{R}-squared for endogenous variables can be requested, as well as the
##'  Fraction Missing Information (FMI) for parameter estimates. By default, the
##'  output will appear like \code{lavaan}'s \code{summary} output, but if
##'  \code{output == "data.frame"}, the returned \code{data.frame} will resemble
##'  the \code{parameterEstimates} output. The \code{scale.W} argument is
##'  passed to \code{vcov} (see description above).
##'  Setting \code{fit.measures=TRUE} will additionally print fit measures to
##'  the console, but they will not be returned; additional arguments may be
##'  passed via \code{...} to \code{\link[lavaan]{fitMeasures}} and
##'  subsequently to \code{\link{lavTestLRT.mi}}.}
##'
##' @section Objects from the Class: See the \code{\link{lavaan.mi}} function
##'   for details. Wrapper functions include \code{\link{cfa.mi}},
##'   \code{\link{sem.mi}}, and \code{\link{growth.mi}}.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Enders, C. K. (2010). \emph{Applied missing data analysis}. New York, NY:
##'   Guilford.
##'
##'   Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}.
##'   New York, NY: Wiley.
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



##' @name lavaan.mi-class
##' @aliases show,lavaan.mi-method
##' @importFrom methods show
##' @export
setMethod("show", "lavaan.mi", function(object) {
  nData <- object@meta$ndat

  useImps <- sapply(object@convergence, "[[", i = "converged")
  nConverged <- sum(useImps)

  SE <- sapply(object@convergence, "[[", "SE")
  SE[is.na(SE)] <- FALSE

  Heywood.ov <- sapply(object@convergence, "[[", "Heywood.ov")
  Heywood.ov[is.na(Heywood.ov)] <- FALSE

  Heywood.lv <- sapply(object@convergence, "[[", "Heywood.lv")
  Heywood.lv[is.na(Heywood.lv)] <- FALSE

  cat('lavaan.mi object based on ', nData, ' imputed data sets. \n',
      'See class?lavaan.mi help page for available methods. \n\n',
      'Convergence information:\n', 'The model converged on ',
      nConverged, ' imputed data sets \n\n', sep = "")

  if (!all(SE)) cat('Standard errors could not be computed for data set(s)',
                    paste(which(!SE), collapse = ", "), '\nTry fitting the',
                    'model to the individual data set(s) to diagnose',
                    'problems. If they cannot be fixed, try inspecting the',
                    'imputations. It may be necessary to reimpute the data',
                    'with some restrictions imposed. \n\n')

  if (any(Heywood.ov | Heywood.lv))
    cat('Heywood cases detected for data set(s)',
        paste(which(Heywood.ov | Heywood.lv), collapse = ", "),
        '\nThese are not necessarily a cause for concern, unless a pooled',
        'estimate is also a Heywood case. \n\n')

  object
})
#TODO: add test stat, like lavaan-class prints



## analog to lavaan:::lav_object_summary(), which creates a list of
## components to appear in the summary() output.  Creating a similar
## object allows lavaan.mi to capitalize on lavaan:::print.lavaan.summary()
lavaan_mi_object_summary <- function(object) {

}
##' @importFrom stats pt qt pnorm qnorm
##' @importFrom lavaan lavListInspect parTable lavNames
##' @importFrom methods getMethod
summary_lavaan_mi <- function(object, header = TRUE,
                              # fit.measures = FALSE,
                              # fm.args =
                              #   list(
                              #     standard.test = "default",
                              #     scaled.test = "default",
                              #     rmsea.ci.level = 0.90,
                              #     rmsea.h0.closefit = 0.05,
                              #     rmsea.h0.notclosefit = 0.08
                              #   ),
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
                              remove.step1 = TRUE,
                              remove.unused = TRUE,
                              # modindices = FALSE,
                              ...) {
  PE <- parameterEstimates.mi(object = object,
                              ci = ci, standardized = standardized,
                              rsquare = rsquare, fmi = fmi,
                              cov.std = cov.std,
                              remove.eq = FALSE, remove.system.eq = TRUE,
                              remove.ineq = FALSE, remove.def = FALSE,
                              remove.nonfree = FALSE,
                              remove.step1 = remove.step1,
                              remove.unused = remove.unused,
                              output = "text",
                              header = TRUE)

  #TODO: move this to a print method for an object class this function returns
  #      (like lav_object_summary() does)
  # if (output == "text") {
  #   getMethod("show", "lavaan.mi")(object)
  #   cat(messPool)
  # }


  # if (fit.measures) {
  #   indices <- c("chisq","df","pvalue","cfi","tli","rmsea","srmr")
  #   FITS <- suppressWarnings(fitMeasures_mi(object, output = "text",
  #                                           fit.measures = indices,
  #                                           fm.args = fm.args, ...))
  #   try(print(FITS, add.h0 = TRUE), silent = TRUE)
  # }

  #TODO: add score-test option: if (modindices) {}

  PE
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

##' @name lavaan.mi-class
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
setMethod("fitted", "lavaan.mi", fitted_lavaan_mi)
##' @name lavaan.mi-class
##' @aliases fitted.values,lavaan.mi-method
##' @export
setMethod("fitted.values", "lavaan.mi", fitted_lavaan_mi)


## utility function called within resid_lavaan_mi() and mi2lavaan()
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

      ## only for fitted() and resid() methods
      if (momentsNblocks  &&  nm %in% c("th","res.th")) {
        ## remove numeric -means from thresholds
        target <- as.numeric(target)[ th.idx[[b]] ]
      }

      if (momentsNblocks) {
        OBS[[b]][[nm]] <- target # used by resid_lavaan_mi()
      } else {
        OBS[[nm]][[b]] <- target # used by mi2lavaan()
      }
      ## end loop over moments
    }
    ## end loop over blocks
  }

  OBS
}



