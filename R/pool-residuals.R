### Terrence D. Jorgensen
### Last updated: 30 April 2024
### pool covariance/correlation residuals
### define resid() method and lavResiduals.mi()

## ---------------------------
## Standard resid(uals) method
## ---------------------------

##' @importFrom lavaan lavListInspect
##' @importFrom stats cov2cor resid residuals
resid_lavaan_mi <- function(object, type = c("raw","cor"),
                            omit.imps = c("no.conv","no.se")) {
  ## @SampleStatsList is (for each imputation) output from:
  ##    getSampStats <- function(obj) lavInspect(obj, "sampstat")
  ## Code below gets equivalent information from @h1List

  useImps <- imps2use(object = object, omit.imps = omit.imps)
  m <- length(useImps)

  ## check type options
  type <- tolower(type[1])
  if (type %in% c("raw","rmr")) {
    type = "raw"
  } else if (type %in% c("cor","cor.bollen","crmr")) {
    type <- "cor.bollen"
  } else if (type %in% c("cor.bentler","cor.eqs","srmr")) {
    type <- "cor.bentler"
  } else stop('type="', type, '" not supported for lavaan.mi objects')

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

  ## H0 (structured model) implied moments
  Implied <- fitted_lavaan_mi(object, momentsNblocks = TRUE, omit.imps = omit.imps)
  if (nBlocks == 1L) Implied <- list(Implied) # store single block in a block-list

  ## H1 (saturated model) implied moments
  OBS <- pool_h1(object, momentsNblocks = TRUE, omit.imps = omit.imps)

  ## template to store observed moments & residuals
  RES <- vector("list", nBlocks)

  ## loop over (blocks and) moments (moments-list nested in block-list)
  for (b in 1:nBlocks) {
    for (nm in names(Implied[[b]])) {

      ## skip if Implied element is not part of the saturated list
      if (is.null(object@h1List[[ useImps[1] ]]$implied[[nm]][[b]])) next

      ## calculate residuals
      if (type == "raw") {
        RES[[b]][[nm]] <- OBS[[b]][[nm]] - Implied[[b]][[nm]]
        class(RES[[b]][[nm]]) <- class(Implied[[b]][[nm]])


        ## correlation residuals
      } else if (type == "cor.bollen") {

        if (nm %in% c("cov","res.cov")) {
          RES[[b]][[nm]] <- cov2cor(OBS[[b]][[nm]]) - cov2cor(Implied[[b]][[nm]])
          class(RES[[b]][[nm]]) <- c("lavaan.matrix.symmetric","matrix")

          ## mean structure
        } else if (nm == "mean") {
          std.obs.M <- OBS[[b]][[nm]] / sqrt(diag(OBS[[b]]$cov))
          std.mod.M <- Implied[[b]][[nm]] / sqrt(diag(Implied[[b]]$cov))
          RES[[b]][[nm]] <- std.obs.M - std.mod.M
          class(RES[[b]][[nm]]) <- c("lavaan.vector","numeric")
        } else if (nm == "res.int") {
          std.obs.M <- OBS[[b]][[nm]] / sqrt(diag(OBS[[b]]$res.cov))
          std.mod.M <- Implied[[b]][[nm]] / sqrt(diag(Implied[[b]]$res.cov))
          RES[[b]][[nm]] <- std.obs.M - std.mod.M
          class(RES[[b]][[nm]]) <- c("lavaan.vector","numeric")

          ## thresholds, slopes, cov.x, mean.x
        } else {
          #FIXME: lavaan currently (0.6-4.1399) returns nothing
          next
        }


        ## standardized (by observed SDs) residuals
      } else if (type == "cor.bentler") {

        if (nm %in% c("cov","mean")) {
          SDs <- diag(sqrt(diag(OBS[[b]]$cov)))
          dimnames(SDs) <- dimnames(OBS[[b]][[nm]])
        } else if (nm %in% c("res.cov","res.int")) {
          SDs <- diag(sqrt(diag(OBS[[b]]$res.cov)))
          dimnames(SDs) <- dimnames(OBS[[b]][[nm]])
        } else {
          #FIXME: lavaan currently (0.6-4.1399) returns nothing for "th" or "slopes"
          next
        }


        if (nm %in% c("cov","res.cov")) {
          RES[[b]][[nm]] <- solve(SDs) %*% (OBS[[b]][[nm]] - Implied[[b]][[nm]]) %*% solve(SDs)
          class(RES[[b]][[nm]]) <- c("lavaan.matrix.symmetric","matrix")
        } else if (nm %in% c("mean","res.int")) {
          RES[[b]][[nm]] <- (OBS[[b]][[nm]] - Implied[[b]][[nm]]) / diag(SDs)
          class(RES[[b]][[nm]]) <- c("lavaan.vector","numeric")
        }

      }

      ## copy names from fitted() results
      if (is.null(dim(RES[[b]][[nm]]))) {
        names(RES[[b]][[nm]]) <- names(Implied[[b]][[nm]])
      } else dimnames(RES[[b]][[nm]]) <- dimnames(Implied[[b]][[nm]])

      ## end loop over moments
    }

    ## add type to beginning of each block's list
    RES[[b]] <- c(list(type = type), RES[[b]])

    #TODO: Rename (res.)cov to (res.)cor?  lavResiduals() does not

  }

  ## drop list for 1 block
  if (nBlocks == 1L) {
    RES <- RES[[1]]
  } else names(RES) <- block.label #FIXME: will lavaan do this in the future?

  RES
}
##' @name lavaan.mi-class
##' @aliases residuals,lavaan.mi-method
##' @export
setMethod("residuals", "lavaan.mi", resid_lavaan_mi)
##' @name lavaan.mi-class
##' @aliases resid,lavaan.mi-method
##' @export
setMethod("resid", "lavaan.mi", resid_lavaan_mi)



## ---------------------------
## Analog for lavResiduals()
## ---------------------------

##' Covariance and Correlation Residuals
##'
##' This function calculates residuals for sample moments (e.g., means and
##' (co)variances, means) from a lavaan model fitted to multiple imputed data
##' sets, along with summary and inferential statistics about the residuals.
##'
##' @importFrom lavaan lavInspect
##'
##' @param object An object of class `lavaan.mi`
##' @param omit.imps `character` indicating criteria for excluding imputations
##'        from pooled results. See [lavaan.mi-class] for argument details.
##' @param \dots Arguments passed to [lavaan::lavResiduals()].
##'
##' @return
##' A `list` of residuals and other information (see [lavaan::lavResiduals()]).
##'
##' @seealso
##' [lavaan::lavResiduals()] for details about other arguments.
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @examples
##'
##' data(HS20imps) # import a list of 20 imputed data sets
##'
##' ## specify CFA model from lavaan's ?cfa help page
##' HS.model <- '
##'   visual  =~ x1 + x2 + x3
##'   textual =~ x4 + x5 + x6
##'   speed   =~ x7 + x8 + x9
##' '
##' ## fit model to 20 imputed data sets
##' fit <- cfa.mi(HS.model, data = HS20imps)
##'
##' ## default type = "cor.bentler" (standardized covariance residuals)
##' lavResiduals.mi(fit, zstat = FALSE)
##' ## SRMR is in the $summary
##'
##' ## correlation residuals
##' lavResiduals.mi(fit, zstat = FALSE, type = "cor")
##' ## CRMR is in the $summary
##'
##' ## raw covariance residuals
##' lavResiduals.mi(fit, type = "raw") # zstat=TRUE by default
##' ## RMR is in the $summary
##' ## "normalized" residuals are in $cov.z
##'
##' @export
lavResiduals.mi <- function(object,
                            omit.imps = c("no.conv","no.se"), ...) {
  ## lavResiduals() doesn't work with PML
  if (lavInspect(object, "options")$estimator == "PML") {
    stop("lavResiduals() is currently unavailable when using PML estimation.")
  }

  ## make fake lavaan object from pooled results
  FIT <- mi2lavaan(object, resid = TRUE, omit.imps = omit.imps)

  ## pass it to lavaan::lavResiduals()
  MC <- match.call()
  MC[[1]] <- quote(lavaan::lavResiduals)
  MC$object <- FIT
  MC$omit.imps <- NULL
  eval(MC)
}



## -----------------------------------------------
## utility function called within fitMeasures.mi()
## -----------------------------------------------

##' @importFrom lavaan lavListInspect lavNames
getSRMR <- function(object, type = "cor.bentler", level = "within",
                    include.mean = TRUE, omit.imps = c("no.conv","no.se")) {
  conditional.x <- lavListInspect(object, "options")$conditional.x
  include.mean <- include.mean && lavListInspect(object, "meanstructure")
  include.diag <- type %in% c("cor.bentler","raw")
  mplus <- type == "mplus"
  if (mplus) type <- "cor.bollen"

  ## how many blocks to loop over
  nG <- lavListInspect(object, "ngroups")
  nlevels <- lavListInspect(object, "nlevels")
  ## save relevant sample sizes
  if (nlevels > 1L && level != "within") {
    n.per.group <- lavListInspect(object, "nclusters") #FIXME: only works for 2 levels
    N <- sum(n.per.group)
  } else {
    n.per.group <- lavListInspect(object, "nobs")
    N <- lavListInspect(object, "ntotal")
  }

  ## grab residuals
  R <- resid_lavaan_mi(object, type = type, omit.imps = omit.imps)
  if (mplus) Rd <- resid_lavaan_mi(object, omit.imps = omit.imps,
                                   type = "cor.bentler")
  ## restructure, if necessary
  if (nG == 1L) {
    loopBlocks <- 1L

    ## extract relevant level
    if (nlevels > 1L) {
      R <- R[[level]]
      if (mplus) Rd <- Rd[[level]]
    }
    ## to loop over blocks
    R <- list(R)
    if (mplus) Rd <- list(Rd)


    ## multiple groups AND multilevel
  } else if (nlevels > 1L) {
    loopBlocks <- 2*(1:nG)
    if (level == "within") loopBlocks <- loopBlocks - 1L
    R <- R[loopBlocks]
    if (mplus) Rd <- Rd[loopBlocks]

  } else loopBlocks <- 1:nG # no restructure necessary for multigroup 1-level models


  ## store vector of squared residuals
  RR <- vector("list", nG)
  for (b in loopBlocks) {
    index <- if (conditional.x) "res.cov" else "cov"

    RR[[b]] <- R[[b]][[index]][lower.tri(R[[b]][[index]], diag = FALSE)]^2
    ## only capture means/variances of numeric modeled variables (not conditional.x)
    vv <- intersect(lavNames(object, type = "ov.num", block = b),
                    lavNames(object, type = "ov.model", block = b))
    if (include.diag)  RR[[b]] <- c(RR[[b]], diag(R[[b]][[index]])[vv]^2)
    if (mplus)  RR[[b]] <- c(RR[[b]], diag(Rd[[b]][[index]])[vv]^2)

    if (include.mean) {
      index <- if (conditional.x) "res.int" else "mean"
      RR[[b]] <- c(RR[[b]], R[[b]][[index]][vv]^2)
    }
  }

  ## take weighted average of group means
  as.numeric( (n.per.group %*% sqrt(sapply(RR, mean))) / N )
}
#TODO: Update according to lav_fit changes?
#      No, develop mi2lavaan() to rely on lavaan's calculations


