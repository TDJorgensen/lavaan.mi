### Terrence D. Jorgensen
### Last updated: 23 April 2024
### Create faux lavaan-class object in order to pass to
### - lavaan::fitMeasures() within fitMeasures_mi()
### - standardizedSolution() within standardizedSolution.mi()
### - lavResiduals() within lavResiduals.mi()


##' @importFrom lavaan lavaan lavaanList
mi2lavaan <- function(object, omit.imps = c("no.conv","no.se"),
                      chi2  = FALSE, # store info for chisq-based indices
                      rmr   = FALSE, # store info for residual-based indices
                      resid = FALSE, # store additional info for lavResiduals()
                      std   = FALSE, # store info for standardizedSolution()
                      ## pass arguments to lavTestLRT.mi()
                      ...) {
  ## this also checks the class
  useImps <- imps2use(object = object, omit.imps = omit.imps)

  if (resid) rmr <- TRUE # lavResiduals() requires same info as RMR

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

  ## assemble saturated-model call?
  if (resid) {
    CALL1 <- object@call
    CALL1[[1]] <- poolSat
    CALL1$model <- NULL
    if (!is.null(CALL1$cmd)) CALL1$cmd <- NULL
    CALL1$omit.imps <- omit.imps
    CALL1$data <- lapply(object@DataList,
                         function(x) x[ c(lavNames(object), object@Data@group) ])
    FIT1 <- eval(as.call(CALL1))
  }

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
      fm.args <- utils::modifyList(default.fm.args, DOTS$fm.args)
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
    TEST <- lavTestLRT.mi(object, omit.imps = omit.imps,
                          ...) # asANOVA = FALSE set in mi_fit_indices_via_lavaan()
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

  ## store pooled summary statistics?
  if (rmr) {
    ## @nobs & @ntotal should already be correct

    FIT@SampleStats@missing.flag <- FALSE # so lavaan doesn't look in @missing.h1

    ## pool unstructured sample moments
    SAT <- pool_h1(object, momentsNblocks = FALSE, omit.imps = omit.imps)

    ## Update @SampleStats slot
    sameNames <- intersect(names(SAT), methods::slotNames(FIT@SampleStats))
    for (nn in sameNames) methods::slot(FIT@SampleStats, nn) <- SAT[[nn]]

    ## Update @h1$impled (same info as @SampleStats)
    sameNames1 <- intersect(names(SAT), names(FIT@h1$implied))
    for (nn in sameNames1) FIT@h1$implied[[nn]] <- SAT[[nn]]

    ## Update @implied slot (analogous info from structured model)
    sameNames0 <- intersect(names(SAT), methods::slotNames(FIT@implied))
    for (nn in sameNames0) methods::slot(FIT@implied, nn) <- SAT[[nn]]

    ## group.w is not part of @SampleStats, so add the pooled one
    if (!is.null(FIT@implied$group.w[[1]])) {
      FIT@implied$group.w <- fitted_lavaan_mi(object, momentsNblocks = FALSE)$group.w
    }


    ## Additional information used by lavResiduals()
        # - In @Data slot:
        #     - @weights (average group weights across imputations? or assume constant?)
        #     - @eXo and @X (only for PML)
        #       - these could vary, making average not categorical
        #       - just don't allow lavResiduals.mi() for PML
    if (resid) {
      ## @x.idx should already be correct

      ## N times the ACOV of the SATURATED model (use poolSat)
      if (is.list(FIT1$NACOV)) {
        FIT@SampleStats@NACOV <- FIT1$NACOV
      } else FIT@SampleStats@NACOV <- list(FIT1$NACOV)

      ## invert COV, unless conditional.x
      if (!FIT@Options$conditional.x) {
        FIT@SampleStats@icov <- lapply(FIT@SampleStats@cov, solve)
      }

      ## Weight matrices for least-squares estimators
      if (FIT@Options$estimator == "ULS") {
        FIT@SampleStats@WLS.V <- lapply(FIT@SampleStats@NACOV,
                                        function(x) diag(1, nrow = nrow(x)))
        FIT@SampleStats@WLS.VD <- lapply(FIT@SampleStats@WLS.V, diag)

        ## other least-squares estimators
      } else if (grepl(pattern = "LS", x = FIT@Options$estimator)) {
        FIT@SampleStats@WLS.V <- lapply(FIT@SampleStats@NACOV,
                                        function(x) solve(diag(diag(x))) )
        FIT@SampleStats@WLS.VD <- lapply(FIT@SampleStats@WLS.V, diag)
      }

      ## sampling weights?
      if (length(FIT@Data@sampling.weights)) {
        #TODO: take (harmonic) mean across imputations? Rescale if sum != 1?

        ## For now, assume constant across imputations (first stored in @Data)
      }
      ## end resid
    }
    ## end rmr
  }

  ##TODO: store anything for standardizedSolution()?
  if (std) {
    ## pooled estimates for standardizedSolution()
    est <- coef_lavaan_mi(object, omit.imps = omit.imps)
    ## update @Model@GLIST
    object@Model <- lavaan::lav_model_set_parameters(object@Model, x = est)
    ## store in fake object
    FIT@Model@GLIST <- object@Model@GLIST

    ## start with existing parameter table (without estimates / start values)
    PT <- data.frame(FIT@ParTable)
    if (!is.null(PT$start)) PT$start <- NULL
    if (!is.null(PT$est  )) PT$est   <- NULL
    if (!is.null(PT$se   )) PT$se    <- NULL
    ## get pooled estimates
    FIT@ParTable$est <- coef_lavaan_mi(object, type = "user")
    # PE <- parameterEstimates.mi(object, se = FALSE, omit.imps = omit.imps,
    #                             remove.system.eq = FALSE, remove.eq = FALSE)
    ## merge pooled estimates into parameter table
    # FIT@ParTable$est <- merge(PT, PE, sort = FALSE, all.x = TRUE, all.y = FALSE)$est


    ## store pooled exogenous covariance matrix?
    if (object@Options$conditional.x) {
      FIT@implied$cov.x <- fitted_lavaan_mi(object, omit.imps = omit.imps,
                                            momentsNblocks = FALSE)$cov.x
    }

    ## ACOV of target model, for delta-method SEs
    FIT@vcov$vcov <- vcov_lavaan_mi(object, omit.imps = omit.imps)
  }

  FIT
}


