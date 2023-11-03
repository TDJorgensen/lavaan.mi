### Terrence D. Jorgensen
### Last updated: 3 November 2023
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
  ## this also checks the class
  useImps <- imps2use(object = object, omit.imps = omit.imps)


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
    TEST <- lavTestLRT.mi(object, omit.imps = omit.imps, ...) # asANOVA = FALSE set in mi_fit_indices_via_lavaan()
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
        FIT@test[[scaled.test]]$trace.UGamma <- mean(sapply(object@testList[useImps],
                                                            function(i) i[[scaled.test]]$trace.UGamma),
                                                     na.rm = TRUE)
      }
    }

    ## add log-likelihoods and information criteria?  (only for ML estimators)
    if (!is.na(TEST["logl"])) {
      FIT@loglik$loglik  <- TEST["logl"]
      FIT@h1$logl$loglik <- TEST["unrestricted.logl"]
      FIT@loglik$AIC     <- TEST["aic"]
      FIT@loglik$BIC     <- TEST["bic"]
      FIT@loglik$BIC2    <- TEST["bic2"]
    }

  }

  ##TODO: store pooled summary statistics?
  if (implied) {

  }

  ##TODO: store anything for standardizedSolution()?

  FIT
}


##' @importFrom lavaan lavTestLRT lavListInspect
mi_fit_indices_via_lavaan <- function(object, fit.measures = "all", baseline.model = NULL,
                                      fm.args = list(standard.test     = "default",
                                                     scaled.test       = "default",
                                                     rmsea.ci.level    = 0.90,
                                                     rmsea.close.h0    = 0.05,
                                                     rmsea.notclose.h0 = 0.08,
                                                     cat.check.pd      = TRUE),
                                      output = "vector", omit.imps = c("no.conv","no.se"),
                                      ...) {
  ## this also checks the class
  useImps <- imps2use(object = object, omit.imps = omit.imps)

  lavoptions <- lavListInspect(object, "options")

  ## necessary to pool CHI-SQUARED for fit indices?
  poolChiSq <- TRUE
  if (object@testList[[ useImps[1] ]][[1]]$test == "none") poolChiSq <- FALSE
  if (!poolChiSq) stop('lavaan ERROR: fit measures not available if test = "none".')
  #TODO: check fit.measures for indices calculated from chi-squared


  ## check for additional arguments
  dots <- list(...)
  if (length(dots)) {
    LRT.names <- intersect(names(dots),
                           union(names(formals(lavTestLRT)),
                                 names(formals(lavTestLRT.mi))))
    dots <- if (length(LRT.names)) dots[LRT.names] else list(asymptotic = TRUE)
  } else dots <- list(asymptotic = TRUE)

  ## check test options (adapted from lavTestLRT.mi, limits duplicate warnings)
  test <- dots$test
  if (is.null(test)) {
    test <- "d4"
  } else test <- tolower(test[1])
  if (tolower(test) %in% c("cm","chan.meng","new.lrt","d4")) test <- "D4"
  if (tolower(test) %in% c("mr","meng.rubin","old.lrt","mplus","d3")) test <- "D3"
  if (tolower(test) %in% c("lmrr","li.et.al","pooled.wald","d2")) test <- "D2"
  if (test %in% c("D3","D4") && !lavoptions$estimator %in% c("ML","PML","FML")) {
    message('"D3" and "D4" only available using maximum likelihood estimation. ',
            'Changed test to "D2".')
    test <- "D2"
  }

  ## check for robust
  test.names <- lavoptions$test
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

  robust <- any(test.names %in% c("satorra.bentler","yuan.bentler",
                                  "yuan.bentler.mplus","scaled.shifted",
                                  "mean.var.adjusted","satterthwaite"))
  if (robust) {
    ## assign pool.robust option to object
    if (is.null(dots$pool.robust)) {
      pool.robust <- formals(lavTestLRT.mi)$pool.robust # default value
    } else {
      pool.robust <- dots$pool.robust # user-specified value
    }
  } else dots$pool.robust <- pool.robust <- FALSE

  scaleshift <- any(test.names == "scaled.shifted")

  if (pool.robust && test %in% c("D3","D4")) {
    message('pool.robust = TRUE is only applicable when test = "D2". ',
            'Changed test to "D2".')
    test <- "D2"
  }

  dots$test <- test

  ## mi2lavaan(...) arguments passed to lavTestLRT.mi() to pool the test stat
  argList <- list(object = object, omit.imps = omit.imps)
  argList <- c(argList, dots) # attach lavTestLRT.mi() arguments
  argList$asANOVA    <- FALSE
  argList$standard.test <- fm.args$standard.test # passed to pairwiseLRT() via ...
  argList$scaled.test   <- fm.args$scaled.test   # passed to pairwiseLRT() via ...
  argList$chi2    <- poolChiSq
  argList$implied <- FALSE
  #TODO: argList$implied <- any(grepl(pattern = "rmr", x = tolower(fit.measures)))

  FIT <- do.call(mi2lavaan, argList)


  ## BASELINE model (if necessary)
  incremental <- c("cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni")
  checkEach <- sapply(incremental, function(i) {
    any(grepl(pattern = i, x = tolower(fit.measures)))
  })
  if (any(checkEach) || fit.measures[1] == "all") {

    if (inherits(baseline.model, "lavaan.mi")) {
      baseFit <- baseline.model

    } else if (inherits(object@external$baseline.model, "lavaan.mi")) {
      baseFit <- object@external$baseline.model

      ## fit default baseline.model using PTb from @baselineList
    } else {
      PTb <- object@baselineList[[ useImps[1] ]]$partable
      PTb[c("est","se")] <- NULL
      # FIXME: shouldn't need this line, but lav_partable_merge() fails when
      #        lavaan:::lav_object_extended() returns a NULL slot instead of "plabel"
      PTb$plabel <- paste0(".p", PTb$id, ".")
      group <- lavListInspect(object, "group")
      if (length(group) == 0L) group <- NULL
      cluster <- lavListInspect(object, "cluster")
      if (length(cluster) == 0L) cluster <- NULL

      ## check for estimator= and test= arguments, to avoid unnecessary
      ## concatenation of "standard"
      ## https://github.com/yrosseel/lavaan/issues/306
      if (!is.null(object@call$estimator)) {
        lavoptions$estimator <- object@call$estimator
      }
      if (!is.null(object@call$test)) {
        lavoptions$test <- object@call$test
      }
      baseFit <- lavaan.mi(model = PTb, data = object@DataList[useImps],
                           group = group, cluster = cluster,
                           test = lavoptions$test, estimator = lavoptions$estimator,
                           fixed.x = lavoptions$fixed.x, se = "none", # to save time
                           conditional.x = lavoptions$conditional.x,
                           ordered = lavListInspect(object, "ordered"),
                           parameterization = lavoptions$parameterization)
    }

    #FIXME? Is length(baseImps) always only as long as length(useImps), not the
    #       original m? (even when passed as argument or stored in @external)?
    argList$object <- baseFit
    argList$omit.imps  <- setdiff(omit.imps, "no.se") # se="none" in baseFit

    BASE <- do.call(mi2lavaan, argList)
  } else BASE <- NULL


  lavaan::fitMeasures(FIT, fit.measures = fit.measures,
                      baseline.model = BASE,
                      fm.args = fm.args, output = output)
}

