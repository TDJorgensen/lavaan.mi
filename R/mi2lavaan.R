### Terrence D. Jorgensen
### Last updated: 1 November 2023
### Create faux lavaan-class object with
### - pooled stat(s) in @test (basline.model in @external)
### - pooled moments in @SampleStats and @implied
      #FIXME: This won't be nearly so easy for ML-SEM
### in order to pass to
### - lavaan::fitMeasures() within fitMeasures.mi()
### - standardizedSolution()
### - lavResiduals()


##' @importFrom lavaan lavaan lavaanList
mi2lavaan <- function(object, omit.imps = c("no.conv","no.se"),
                      chi2 = TRUE, implied = TRUE,
                      ## pass arguments to lavTestLRT.mi()
                      ...) {
  stopifnot(inherits(object, "lavaan.mi"))

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


  DOTS <- list(...)

  ## assemble a lavaan() call
  CALL             <- list(quote(lavaan::lavaan))
  CALL$model       <- data.frame(object@ParTable)
  CALL$model$start <- coef_lavaan_mi(object, type = "all", # pooled estimates
                                     omit.imps = omit.imps)
  CALL$data        <- object@DataList[[ useImps[1] ]]
  ## copy lavaan() arguments and lavOptions()
  lavArgs <- setdiff(names(object@lavListCall)[-1],
                     names(formals(lavaan::lavaanList)))
  if (length(lavArgs)) CALL[lavArgs] <- object@lavListCall[lavArgs]
  # ## don't actually fit, but force it to be recognized as converged
  # CALL$do.fit <- FALSE
  # CALL$optim.force.converged <- TRUE
  ## don't auto-fit baseline.model
  ## (instead, run this function separately on fitted baseline.model)
  CALL$baseline <- FALSE
  FIT <- eval(as.call(CALL))

  ## store pooled test statistic(s)?
  if (chi2) {
    ## default fm.args
    default.fm.args <- list(standard.test     = "default",
                            scaled.test       = "default",
                            rmsea.ci.level    = 0.90,
                            rmsea.close.h0    = 0.05,
                            rmsea.notclose.h0 = 0.08,
                            robust            = TRUE,
                            cat.check.pd      = TRUE)
    if (!is.null(DOTS$fm.args)) {
      fm.args <- modifyList(default.fm.args, DOTS$fm.args)
    } else {
      fm.args <- default.fm.args
    }

    ## check requested test stats
    test.names <- unname(sapply(FIT@test, "[[", "test"))
    if (test.names[1] == "none") {
      stop("lavaan ERROR: fit measures not available if test = \"none\".")
    }
    standard.test <- "standard"
    scaled.test   <- fm.args$scaled.test
    ## do we have a scaled test statistic? if so, which one?
    scaled.flag <- FALSE
    if (scaled.test != "none" &&
        any(test.names %in% c("satorra.bentler",
                              "yuan.bentler", "yuan.bentler.mplus",
                              "mean.var.adjusted", "scaled.shifted"))) {
      scaled.flag <- TRUE
      if (scaled.test %in% c("standard", "default")) {
         tmp.idx <- which(test.names %in% c("satorra.bentler",
                                            "yuan.bentler", "yuan.bentler.mplus",
                                            "mean.var.adjusted", "scaled.shifted"))
        scaled.test <- test.names[tmp.idx[1]]
      }
    }
    ## remove additional (scaled) tests from list
    if (scaled.flag) {
      if (length(object@testList[[ useImps[1] ]]) > 2L)
        for (imp in useImps) {
        object@testList[[imp]] <- object@testList[[imp]][c(standard.test, scaled.test)]
      }
    }


    ## pool model's test stats
    TEST <- lavTestLRT.mi(object, omit.imps = omit.imps, asANOVA = FALSE, ...)
    FIT@test$standard$stat.group <- NULL
    FIT@test$standard$stat       <- TEST[["chisq"]]
    FIT@test$standard$df         <- TEST[["df"]]
    FIT@test$standard$pvalue     <- TEST[["pvalue"]]

    if (scaled.flag) {
      FIT@test[[scaled.test]]$stat.group <- NULL
      FIT@test[[scaled.test]]$stat       <- TEST[["chisq.scaled"]]
      FIT@test[[scaled.test]]$df         <- TEST[["df.scaled"]]
      FIT@test[[scaled.test]]$pvalue     <- TEST[["pvalue.scaled"]]
      FIT@test[[scaled.test]]$scaled.test.stat  <- TEST[["chisq"]]
      FIT@test[[scaled.test]]$scaling.factor    <- TEST[["chisq.scaling.factor"]]
      if (scaled.test == "scaled.shifted") {
        FIT@test$scaled.shifted$shift.parameter <- TEST[["chisq.shift.parameter"]]
      } else if (scaled.test == "mean.var.adjusted") {
        ## only necessary for rmsea.scaled
        FIT@test[[scaled.test]]$trace.UGamma <- mean(sapply(object@testList[[useImps]],
                                                            function(i) i[[scaled.test]]$trace.UGamma),
                                                     na.rm = TRUE)
      }
    }

    ##TODO: pool baseline model's test stats
    #       or call this function again to store in @external?

  }

  ##TODO: store pooled summary statistics?
  if (implied) {

  }

  ##TODO: store anything for standardizedSolution()?

  FIT
}



#TODO: mi_fit_indices_via_lavaan <- function(object) {}

