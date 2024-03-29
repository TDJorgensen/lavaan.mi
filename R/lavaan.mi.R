### Terrence D. Jorgensen
### Last updated: 1 November 2022
### function that creates lavaan.mi object, inherits from lavaanList class


##' Fit a lavaan Model to Multiple Imputed Data Sets
##'
##' This function fits a lavaan model to a list of imputed data sets.
##'
##'
##' @aliases lavaan.mi cfa.mi sem.mi growth.mi
##'
##' @param model The analysis model can be specified using lavaan
##'   \code{\link[lavaan]{model.syntax}} or a parameter table (as returned by
##'   \code{\link[lavaan]{parTable}}).
##' @param data A a \code{list} of imputed data sets, or an object class from
##'   which imputed data can be extracted. Recognized classes are
##'   \code{lavaan.mi} (stored in the \code{@@DataList} slot),
##'   \code{amelia} (created by the Amelia package), or
##'   \code{mids} (created by the mice package).
##' @param \dots additional arguments to pass to \code{\link[lavaan]{lavaan}} or
##'   \code{\link[lavaan]{lavaanList}}. See also \code{\link[lavaan]{lavOptions}}.
##'   Note that \code{lavaanList} provides parallel computing options, as well as
##'   a \code{FUN=} argument so the user can extract custom output after the model
##'   is fitted to each imputed data set (see \strong{Examples}).  TIP: If a
##'   custom \code{FUN} is used \emph{and} \code{parallel = "snow"} is requested,
##'   the user-supplied function should explicitly call \code{library} or use
##'   \code{\link[base]{::}} for any functions not part of the base distribution.
##'
##' @return A \code{\linkS4class{lavaan.mi}} object
##'
##' @note This functionality was originally provided via \code{runMI()} in the
##'   \code{semTools} package, but there are differences.  See the README file
##'   on the GitHub page for this package (find link in DESCRIPTION).
##'
##' @author
##'   Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Enders, C. K. (2010). \emph{Applied missing data analysis}.
##'   New York, NY: Guilford.
##'
##'   Rubin, D. B. (1987). \emph{Multiple imputation for nonresponse in surveys}.
##'   New York, NY: Wiley.
##'
##' @examples
##'  \dontrun{
##' ## impose missing data for example
##' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
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
##' ## specify CFA model from lavaan's ?cfa help page
##' HS.model <- '
##'   visual  =~ x1 + x2 + x3
##'   textual =~ x4 + x5 + x6
##'   speed   =~ x7 + x8 + x9
##' '
##'
##' ## fit model to imputed data sets
##' fit <- cfa.mi(HS.model, data = imps)
##'
##' summary(fit, fit.measures = TRUE)
##' summary(fit, ci = FALSE, fmi = TRUE, output = "data.frame")
##' summary(fit, ci = FALSE, stand = TRUE, rsq = TRUE)
##'
##' ## model fit. D3 includes information criteria
##' anova(fit)
##' ## equivalently:
##' lavTestLRT.mi(fit)
##' ## request D2
##' anova(fit, test = "D2")
##' ## request fit indices
##' fitMeasures(fit)
##'
##'
##' ## fit multigroup model without invariance constraints
##' mgfit.config <- cfa.mi(HS.model, data = imps, estimator = "mlm",
##'                        group = "school")
##' ## add invariance constraints, and use previous fit as "data"
##' mgfit.metric <- cfa.mi(HS.model, data = mgfit.config, estimator = "mlm",
##'                        group = "school", group.equal = "loadings")
##' mgfit.scalar <- cfa.mi(HS.model, data = mgfit.config, estimator = "mlm",
##'                        group = "school",
##'                        group.equal = c("loadings","intercepts"))
##'
##' ## compare fit of 2 models to test metric invariance
##' ## (scaled likelihood ratio test)
##' lavTestLRT.mi(mgfit.config, mgfit.metric, mgfit.scalar, # or use anova()
##'               method = "satorra.bentler.2010") # pass argument to lavTestLRT()
##'
##' ## correlation residuals to investigate local misfit
##' resid(mgfit.scalar, type = "cor.bentler")
##' ## modification indices for fixed parameters, to investigate local misfit
##' modindices.mi(mgfit.scalar)
##' ## or lavTestScore.mi for modification indices about equality constraints
##' lavTestScore.mi(mgfit.scalar)
##'
##' ## Wald test of whether latent means are == (fix 3 means to zero in group 2)
##' eq.means <- ' .p70. == 0
##'               .p71. == 0
##'               .p72. == 0 '
##' lavTestWald.mi(mgfit.scalar, constraints = eq.means)
##'
##'
##'
##' ## ordered-categorical data
##' HSbinary <- as.data.frame( lapply(HSMiss[ , paste0("x", 1:9)],
##'                                   FUN = cut, breaks = 2, labels = FALSE) )
##' HSbinary$school <- HSMiss$school
##'
##' ## impute binary missing data using mice package
##' library(mice)
##' set.seed(456)
##' miceImps <- mice(HSbinary)
##' ## save imputations in a list of data.frames
##' impList <- list()
##' for (i in 1:miceImps$m) impList[[i]] <- complete(miceImps, action = i)
##'
##' ## fit model
##' catout <- cfa.mi(HS.model, data = impList, # can also pass data = miceImps
##'                  # use lavaanList(FUN=) argument to save zero-cell tables
##'                  # and obsolete "WRMR" fit index per imputation
##'                  FUN = function(fit) {
##'                    list(wrmr = lavaan::fitMeasures(fit, "wrmr"),
##'                         zeroCells = lavaan::lavInspect(fit, "zero.cell.tables"))
##'                  })
##' summary(catout)
##' lavTestLRT.mi(catout, test = "D2", pool.robust = TRUE)
##' fitMeasures(catout, fit.measures = c("rmsea","srmr","cfi"),
##'             test = "D2", pool.robust = TRUE)
##'
##' ## extract custom output
##' sapply(catout@funList, function(x) x$wrmr) # WRMR for each imputation
##' catout@funList[[1]]$zeroCells # zero-cell tables for first imputation
##' catout@funList[[2]]$zeroCells # zero-cell tables for second imputation ...
##'
##' }
##'
##' @importFrom lavaan lavInspect parTable
##' @importFrom methods as
##' @export
lavaan.mi <- function(model, data, ...) {
  CALL <- match.call()
  dots <- list(...)
  if (is.null(dots$cmd)) dots$cmd <- "lavaan"

  ## check for (Bollen-Stine) bootstrap request
  if (all(!is.null(dots$test),
          tolower(dots$test) %in% c("boot","bootstrap","bollen.stine")) ||
      all(!is.null(dots$se), tolower(dots$se) %in% c("boot","bootstrap"))) {
    stop('Bootstraping unavailable (and not recommended) in combination with ',
         'multiple imputations. For robust confidence intervals of indirect',
         ' effects, see the ?semTools::monteCarloCI help page. To bootstrap ',
         'within each imputation, users can pass a custom function to the ',
         'FUN= argument (see ?lavaanList) to save bootstrap distributions in ',
         'the @funList slot, then manually combine afterward.')
  }

  ## Validate list of imputed data sets
  imputedData <- NULL
  if (missing(data)) {
    #TODO: check for summary statistics
    #TODO: make lavaanList() accept lists of summary stats
    #TODO: Add argument to implement Li Cai's pool-polychorics first, pass
    #      to lavaan for DWLS with pooled WLS.V= and NACOV=, return(lavaan).

  } else if (is.data.frame(data)) {
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

  ## Function to get custom output for lavaan.mi object
  ## NOTE: Need "lavaan::" to allow for parallel computations
  .getOutput. <- function(obj) {
    converged <- lavaan::lavInspect(obj, "converged")
    if (converged) {
      se <- lavaan::parTable(obj)$se
      se.test <- all(!is.na(se)) & all(se >= 0) & any(se != 0)
      if (lavaan::lavInspect(obj, "ngroups") == 1L && lavaan::lavInspect(obj, "nlevels") == 1L) {
        Heywood.lv <- det(lavaan::lavInspect(obj, "cov.lv")) <= 0
        Heywood.ov <- det(lavaan::lavInspect(obj, "theta")) <= 0
      } else {
        Heywood.lv <- !all(sapply(lavaan::lavInspect(obj, "cov.lv"), det) > 0)
        Heywood.ov <- !all(sapply(lavaan::lavInspect(obj, "theta"), det) > 0)
      }
      suppressWarnings(MIs <- try(lavaan::modindices(obj), silent = TRUE))

    } else {
      se.test <- Heywood.lv <- Heywood.ov <- NA
      MIs <- NULL
    }

    list(sampstat = lavaan::lavInspect(obj, "sampstat"),
         coefMats = lavaan::lavInspect(obj, "est"),
         satPT = data.frame(lavaan::lav_partable_unrestricted(obj),
                            #FIXME: do starting values ALWAYS == estimates?
                            stringsAsFactors = FALSE),
         modindices = MIs,
         cov.lv = lavaan::lavInspect(obj, "cov.lv"), #TODO: calculate from pooled estimates for reliability()
         converged = converged, SE = se.test,
         Heywood.lv = Heywood.lv, Heywood.ov = Heywood.ov)
  }

  ## fit model using lavaanList
  lavListCall <- c(list(quote(lavaan::lavaanList), model = model,
                        dataList = imputedData), dots)
  lavListCall$store.slots <- c("partable","vcov","test","h1","baseline")
  lavListCall$FUN <- if (is.null(dots$FUN)) .getOutput. else function(obj) {
    temp1 <- .getOutput.(obj)
    temp2 <- dots$FUN(obj)
    if (!is.list(temp2)) temp2 <- list(userFUN1 = temp2)
    if (is.null(names(temp2))) names(temp2) <- paste0("userFUN", 1:length(temp2))
    duplicatedNames <- which(sapply(names(temp2), function(x) {
      x %in% c("sampstat","coefMats","satPT","modindices","converged",
               "SE","Heywood.lv","Heywood.ov","cov.lv")
    }))
    for (i in duplicatedNames) names(temp2)[i] <- paste0("userFUN", i)
    c(temp1, temp2)
  }
  fit <- eval(as.call(lavListCall))
  ## Store custom @DataList and @SampleStatsList
  fit@SampleStatsList <- lapply(fit@funList, "[[", i = "sampstat")
  fit@DataList <- imputedData
  ## add parameter table to @h1List
  for (i in 1:m) fit@h1List[[i]] <- c(fit@h1List[[i]],
                                      list(PT = fit@funList[[i]]$satPT))
  ## assign class and add new slots
  fit <- as(fit, "lavaan.mi")
  fit@coefList <- lapply(fit@funList, "[[", i = "coefMats")
  fit@miList <- lapply(fit@funList, "[[", i = "modindices")
  fit@phiList <- lapply(fit@funList, "[[", i = "cov.lv")
  fit@call <- CALL
  fit@lavListCall <- lavListCall
  convList <- lapply(fit@funList, "[", i = c("converged","SE",
                                             "Heywood.lv","Heywood.ov"))
  nonConv <- which(sapply(convList, is.null))
  if (length(nonConv)) for (i in nonConv) {
    convList[[i]] <- list(converged = FALSE, SE = NA, Heywood.lv = NA, Heywood.ov = NA)
  }
  fit@convergence <- lapply(convList, function(x) do.call(c, x))
  conv <- which(sapply(fit@convergence, "[", i = "converged"))
  if (!length(conv)) warning('The model did not converge for any imputed data sets.')

  ## keep any remaining funList slots (if allowing users to supply custom FUN)
  funNames <- names(fit@funList[[1]])
  keepIndex <- which(!sapply(funNames, function(x) {
    x %in% c("sampstat","coefMats","satPT","modindices","converged",
             "SE","Heywood.lv","Heywood.ov","cov.lv")
  }))
  if (length(keepIndex)) {
    fit@funList <- lapply(fit@funList, "[", i = keepIndex)
    if (length(keepIndex) > 1L) {
      keepNames <- funNames[keepIndex]
      noNames <- which(keepNames == "")
      for (i in seq_along(noNames)) keepNames[ noNames[i] ] <- paste0("userFUN", i)
      fit@funList <- lapply(fit@funList, "names<-", value = keepNames)
    }
  } else fit@funList <- list()

  NewStartVals <- try(coef_lavaan_mi(fit, type = "user", labels = FALSE),
                           silent = TRUE)
  if (!inherits(NewStartVals, "try-error")) fit@ParTable$start <- NewStartVals
  #FIXME: else do what? warn?

  fit
}

##' @rdname lavaan.mi
##' @export
cfa.mi <- function(model, data, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc$cmd <- "cfa"
  mc[[1L]] <- quote(lavaan.mi::lavaan.mi)
  eval(mc, parent.frame())
}

##' @rdname lavaan.mi
##' @export
sem.mi <- function(model, data, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc$cmd <- "sem"
  mc[[1L]] <- quote(lavaan.mi::lavaan.mi)
  eval(mc, parent.frame())
}

##' @rdname lavaan.mi
##' @export
growth.mi <- function(model, data, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc$cmd <- "growth"
  mc[[1L]] <- quote(lavaan.mi::lavaan.mi)
  eval(mc, parent.frame())
}


