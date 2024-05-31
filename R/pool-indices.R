### Terrence D. Jorgensen
### Last updated: 28 May 2024
### pool fit indices
### define fitMeasures() method for lavaan.mi


##' @importFrom lavaan lavTestLRT lavListInspect
fitMeasures_mi <- function(object, fit.measures = "all",
                           baseline.model = NULL, h1.model = NULL,
                           fm.args = list(standard.test        = "default",
                                          scaled.test          = "default",
                                          rmsea.ci.level       = 0.90,
                                          rmsea.h0.closefit    = 0.05,
                                          rmsea.h0.notclosefit = 0.08,
                                          robust               = 0.08,
                                          cat.check.pd         = TRUE),
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
    if (length(LRT.names)) dots <- dots[LRT.names]
  }
  dots$asymptotic <- TRUE # ALWAYS for fit indices

  ## check test options (adapted from lavTestLRT.mi, limits duplicate warnings)
  pool.method <- dots$pool.method
  if (is.null(pool.method)) {
    pool.method <- "d4"
  } else pool.method <- tolower(pool.method[1])
  if (tolower(pool.method) %in% c("cm","chan.meng","new.lrt","d4")) pool.method <- "D4"
  if (tolower(pool.method) %in% c("mr","meng.rubin","old.lrt","mplus","d3")) pool.method <- "D3"
  if (tolower(pool.method) %in% c("lmrr","li.et.al","pooled.wald","d2")) pool.method <- "D2"
  if (pool.method %in% c("D3","D4") && !lavoptions$estimator %in% c("ML","PML","FML")) {
    message('"D3" and "D4" only available using maximum likelihood estimation. ',
            'Changed to pool.method = "D2".')
    pool.method <- "D2"
  }

  ## check for robust
  test.names <- lavoptions$test
  #FIXME: lavaan 0.6-5: for now, we only acknowledge the first non-standard @test
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

  if (pool.robust && pool.method %in% c("D3","D4")) {
    message('pool.robust = TRUE is only applicable when pool.method = "D2". ',
            'Changed to pool.method = "D2".')
    pool.method <- "D2"
  }

  dots$pool.method <- pool.method

  ## mi2lavaan(...) arguments passed to lavTestLRT.mi() to pool the test stat
  argList <- list(object = object, omit.imps = omit.imps)
  argList <- c(argList, dots) # attach lavTestLRT.mi() arguments
  argList$asANOVA    <- FALSE
  argList$standard.test <- fm.args$standard.test # passed to pairwiseLRT() via ...
  argList$scaled.test   <- fm.args$scaled.test   # passed to pairwiseLRT() via ...
  argList$chi2    <- poolChiSq
  argList$rmr <- ( fit.measures[1] %in% c("all","default")   ||
                   any(grepl(pattern = "rmr", x = tolower(fit.measures))) )

  FIT <- do.call(mi2lavaan, argList)


  ## BASELINE model (if necessary)
  incremental <- c("cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni")
  checkEach <- sapply(incremental, function(i) {
    any(grepl(pattern = i, x = tolower(fit.measures)))
  })
  if (any(checkEach) || fit.measures[1] %in% c("all","default")) {

    if (inherits(baseline.model, "lavaan.mi")) {
      baseFit <- baseline.model

    } else if (inherits(object@external$baseline.model, "lavaan.mi")) {
      baseFit <- object@external$baseline.model

      ## fit default baseline.model using PTb from @baselineList
    } else {
      if (is.null(object@baselineList[[ useImps[1] ]]$partable)) {
        stop('No baseline.model provided or fitted. Fit a custom baseline.model',
             ' or set baseline=TRUE when fitting the target model.')
      }
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
      ## Yves updated lavaan to avoid this.  Unnecessary now?
      # if (!is.null(object@call$estimator)) {
      #   ESTIMATOR <- try(object@call$estimator, silent = TRUE)
      #   # evalq(meth, envir = parent.frame())
      #   lavoptions$estimator <- object@call$estimator
      # }
      # if (!is.null(object@call$test)) {
      #   lavoptions$test <- object@call$test
      # }
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


  #FIXME: This will NOT yield a pooled chisq-difference test.
  #       It will calculate the difference between 2 pooled chisq stats.
  ## custom h1.model (if necessary)
  # if (inherits(h1.model, "lavaan.mi")) {
  #   h1Fit <- h1.model
  # } else if (inherits(object@external$h1.model, "lavaan.mi")) {
  #   h1Fit <- object@external$h1.model
  # } else h1Fit <- NULL
  #
  # if (is.null(h1Fit)) {
  #   H1 <- NULL
  # } else {
  #   argList$object <- h1Fit
  #   argList$omit.imps  <- omit.imps # in case it was changed for BASE
  #   H1 <- do.call(mi2lavaan, argList)
  # }

  ## extract names of pooled tests
  standard.test <- FIT@external$mi2lavaan$standard.test
  scaled.test   <- FIT@external$mi2lavaan$scaled.test

  ## calculate fit measures
  OUT <- lavaan::fitMeasures(FIT, fit.measures = fit.measures,
                             baseline.model = BASE,
                             h1.model = NULL, #FIXME: Possible? Necessary?
                             fm.args = fm.args, output = output)

  scaled.flag   <- any(grepl("scaled", names(OUT)))

  ## add header
  HEADER <- paste0("Test statistic(s) pooled using the ", pool.method,
                   " pooling method.\n  Pooled statistic: ",
                   ifelse(pool.robust && pool.method == "D2",
                          yes = dQuote(scaled.test), no = dQuote(standard.test)))
  if (scaled.flag) {
    ## indicate whether pool.robust
    HEADER <- paste0(HEADER,
                     ifelse(pool.method == "D2",
                            yes = paste0("  (pool.robust=", pool.robust, ")"),
                            no = ""))
    if (!pool.robust) {
      ## indicate the scaling method
      HEADER <- paste0(HEADER, "\n  Method to robustify pooled statistic:  ",
                       dQuote(scaled.test))
    }
  }
  attr(OUT, "header") <- HEADER

  ## give the header more attributes, ignored by lavaan:::print.lavaan.vector(),
  ## but caught by lavaan:::print.lavaan.fitMeasures() for summary()
  attr(attr(OUT, "header"), "standard.test") <- standard.test
  attr(attr(OUT, "header"),   "scaled.test") <-   scaled.test
  if (isTRUE(scaled.test == "scaled.shifted")) {
    ## not in the vector, save it as an attribute
    attr(attr(OUT, "header"), "shift") <- FIT@external$mi2lavaan$shift
  }

  OUT
}

##' @name lavaan.mi-class
##' @aliases fitMeasures,lavaan.mi-method
##' @importFrom lavaan fitMeasures
##' @export
setMethod("fitMeasures", "lavaan.mi", fitMeasures_mi)
## lowercase 'm'
##' @name lavaan.mi-class
##' @aliases fitmeasures,lavaan.mi-method
##' @importFrom lavaan fitmeasures
##' @export
setMethod("fitmeasures", "lavaan.mi", fitMeasures_mi)



#FIXME: delete when this is no longer necessary
##' @importFrom lavaan lavNames lavListInspect
##' @importFrom stats pchisq uniroot
fitMeas_manual <- function(object, fit.measures = "all",
                           baseline.model = NULL, h1.model = NULL,
                           fm.args = list(standard.test        = "default",
                                          scaled.test          = "default",
                                          rmsea.ci.level       = 0.90,
                                          rmsea.h0.closefit    = 0.05,
                                          rmsea.h0.notclosefit = 0.08,
                                          robust               = 0.08,
                                          cat.check.pd         = TRUE),
                           output = "vector", omit.imps = c("no.conv","no.se"),
                           ...) {

  useImps <- imps2use(object = object, omit.imps = omit.imps)
  m <- length(useImps)

  lavoptions <- lavListInspect(object, "options")

  fit.measures <- tolower(fit.measures)
  if (length(fit.measures) == 0L) fit.measures <- "all"
  ## narrow down fit indices
  incremental <- c("cfi","tli","nnfi","rfi","nfi","pnfi","ifi","rni")
  if ("all" %in% fit.measures) {
    indices <- c("chisq","df","pvalue","scaling", incremental,
                 "rmsea","rmr","mfi","gammahat")
  } else {
    indices <- grep(pattern = paste(c("chisq","df","pvalue","scaling",
                                      incremental, "mfi","rmsea",
                                      "gammahat","rmr"), collapse = "|"),
                    x = fit.measures, ignore.case = TRUE, value = TRUE)
  }

  ## CHI-SQUARED-BASED FIT INDICES
  notest <- length(lavoptions$test) == 1L && lavoptions$test == "none"
  if (notest || any(!grepl(pattern = "rmr", x = indices))) {

    ## check for additional arguments
    dots <- list(...)
    if (length(dots)) {
      LRT.names <- intersect(names(dots),
                             union(names(formals(lavTestLRT)),
                                   names(formals(lavTestLRT.mi))))
      dots <- if (length(LRT.names)) dots[LRT.names] else list(asymptotic = TRUE)
    } else dots <- list(asymptotic = TRUE)

    ## check test options (adapted from lavTestLRT.mi, limits duplicate warnings)
    pool.method <- dots$pool.method
    if (is.null(pool.method)) {
      pool.method <- "d4"
    } else pool.method <- tolower(pool.method[1])
    if (tolower(pool.method) %in% c("cm","chan.meng","new.lrt","d4")) pool.method <- "D4"
    if (tolower(pool.method) %in% c("mr","meng.rubin","old.lrt","mplus","d3")) pool.method <- "D3"
    if (tolower(pool.method) %in% c("lmrr","li.et.al","pooled.wald","d2")) pool.method <- "D2"
    if (pool.method %in% c("D3","D4") && !grepl(pattern = "ML", x = lavoptions$estimator)) {
      message('"D3" and "D4" only available using maximum likelihood estimation. ',
              'Changed to pool.method = "D2".')
      pool.method <- "D2"
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
    # if (scaleshift) {
    #   if (pool.method %in% c("D3","D4")) {
    #     message("If test = 'scaled.shifted' (estimator = 'WLSMV' or 'MLMV'), ",
    #             "model evaluation is only available by (re)setting ",
    #             "pool.method = 'D2'.\nControl more options by passing arguments to ",
    #             "lavTestLRT() via the '...' argument.\n")
    #     pool.method <- 'D2'
    #   }
    # }

    if (pool.robust && pool.method %in% c("D3","D4")) {
      message('pool.robust = TRUE is only applicable when pool.method = "D2". ',
              'Changed to pool.method = "D2".')
      pool.method <- "D2"
    }

    dots$pool.method <- pool.method


    ## pooled test statistic(s)
    argList <- c(list(object = object), dots)
    argList$asymptotic <- TRUE # in case it wasn't set in list(...)
    argList$omit.imps  <- omit.imps
    argList$asANOVA    <- FALSE
    argList$standard.test <- fm.args$standard.test
    argList$scaled.test   <- fm.args$scaled.test
    out <- do.call(lavTestLRT.mi, argList)
    ## check for scaled test statistic (if not, set robust=FALSE)
    if (robust && is.na(out["chisq.scaled"])) robust <- FALSE

    ## fit baseline model if necessary
    if (any(indices %in% incremental)) {
      if (inherits(baseline.model, "lavaan.mi")) {
        baseFit <- baseline.model
      } else if (inherits(object@external$baseline.model, "lavaan.mi")) {
        baseFit <- object@external$baseline.model

        ## MUST fit PTb for D3 or D4 likelihoods, but for D2 use @baselineList
      } else if (pool.method == "D2") {
        ## length(baseImps) == m, not just length(useImps)
        baseImps <- object@meta$baseline.ok
        if (!all(baseImps[useImps])) warning('The default independence model ',
                                             'did not converge for data set(s): ',
                                             which(!baseImps))
        ## only use imputations that converged for both
        baseImps <- intersect(useImps, which(baseImps))

        w <- sapply(object@baselineList[baseImps],
                    function(x) x$test$standard[["stat"]])
        if (is.list(w)) {
          #TODO: figure out why this happens!
          w <- unlist(w)
          DF <- mean(unlist(sapply(object@baselineList[baseImps],
                                   function(x) x$test$standard[["df"]])))
        } else {
          DF <- mean(sapply(object@baselineList[baseImps],
                            function(x) x$test$standard[["df"]]))
        }
        baseOut <- calculate.D2(w, DF, asymptotic = TRUE)
        if (robust) {
          if (pool.robust) {
            w.r <- sapply(object@baselineList[baseImps],
                          function(x) x$test[[ test.names[1] ]][["stat"]])
            if (is.list(w.r)) {
              w.r <- unlist(w.r)
              DF.r <- mean(unlist(sapply(object@baselineList[baseImps],
                                         function(x) x$test[[ test.names[1] ]][["df"]])))
            } else {
              DF.r <- mean(sapply(object@baselineList[baseImps],
                                  function(x) x$test[[ test.names[1] ]][["df"]]))
            }
            base.robust <- calculate.D2(w.r, DF.r, asymptotic = TRUE)
            names(base.robust) <- paste0(names(base.robust), ".scaled")
            baseOut <- c(baseOut, base.robust)
          } else {
            baseOut <- robustify(ChiSq = baseOut, object = object,
                                 baseline = TRUE, useImps = baseImps)
          }
        }
        baseFit <- NULL # for later checking, to avoid unnecessary calls

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
        baseFit <- lavaan.mi(model = PTb, data = object@DataList[useImps],
                             group = group, cluster = cluster,
                             test = lavoptions$test, estimator = lavoptions$estimator,
                             fixed.x = lavoptions$fixed.x, se = "none", # to save time
                             conditional.x = lavoptions$conditional.x,
                             ordered = lavListInspect(object, "ordered"),
                             parameterization = lavoptions$parameterization)
      }

      if (!is.null(baseFit)) {
        #FIXME? Is length(baseImps) always only as long as length(useImps), not the
        #       original m? (even when passed as argument or stored in @external)?
        baseImps <- sapply(baseFit@convergence, "[[", i = "converged")
        if (!all(baseImps)) warning('baseline.model did not converge for data set(s): ',
                                    useImps[!baseImps])
        argList <- c(list(object = baseFit), dots)
        argList$asymptotic <- TRUE # in case it wasn't set in list(...)
        argList$omit.imps  <- setdiff(omit.imps, "no.se") # se="none" in baseFit
        argList$asANOVA    <- FALSE
        argList$standard.test <- fm.args$standard.test
        argList$scaled.test   <- fm.args$scaled.test
        baseOut <- do.call(lavTestLRT.mi, argList)
      }
      # else { already used "D2" with @baselineList info to make baseOut }

    }


    X2 <- out[["chisq"]]
    DF <- out[["df"]]
    if (robust) {
      X2.sc <- out[["chisq.scaled"]]
      DF.sc <- out[["df.scaled"]] ## for mean.var.adjusted, mean DF across imputations
      if (!pool.robust) ch <- out[["chisq.scaling.factor"]] ## mean c_hat across imputations
      if (X2 < .Machine$double.eps && DF == 0) ch <- 0
      ## for RMSEA
      if ("rmsea" %in% indices) {
        d <- mean(sapply(object@testList[useImps],
                         function(x) sum(x[[ test.names[1] ]][["trace.UGamma"]])))
        if (is.na(d) || d == 0) d <- NA # FIXME: only relevant when mean.var.adjusted?
      }
    }

    ## for CFI, TLI, etc.
    if (any(indices %in% incremental)) {
      bX2 <- baseOut[["chisq"]]
      bDF <- baseOut[["df"]]
      out <- c(out, baseline.chisq = bX2, baseline.df = bDF,
               baseline.pvalue = baseOut[["pvalue"]])
      if (robust) {
        out["baseline.chisq.scaled"] <- bX2.sc <- baseOut[["chisq.scaled"]]
        out["baseline.df.scaled"]    <- bDF.sc <- baseOut[["df.scaled"]]
        out["baseline.pvalue.scaled"] <- baseOut[["pvalue.scaled"]]
        if (!pool.robust) {
          cb <- baseOut[["chisq.scaling.factor"]]
          out["baseline.chisq.scaling.factor"] <- cb
          if (scaleshift) out["baseline.chisq.shift.parameter"] <- baseOut[["chisq.shift.parameter"]]
        }
      }
    }

    if ("cfi" %in% indices) {
      t1 <- max(X2 - DF, 0)
      t2 <- max(X2 - DF, bX2 - bDF, 0)
      out["cfi"] <- if(t1 == 0 && t2 == 0) 1 else 1 - t1/t2
      if (robust) {
        ## scaled
        t1 <- max(X2.sc - DF.sc, 0)
        t2 <- max(X2.sc - DF.sc, bX2.sc - bDF.sc, 0)
        if (is.na(t1) || is.na(t2)) {
          out["cfi.scaled"] <- NA
        } else if (t1 == 0 && t2 == 0) {
          out["cfi.scaled"] <- 1
        } else out["cfi.scaled"] <- 1 - t1/t2
        ## Brosseau-Liard & Savalei MBR 2014, equation 15
        if (!pool.robust & test.names[1] %in%
            c("satorra.bentler","yuan.bentler","yuan.bentler.mplus")) {
          t1 <- max(X2 - ch*DF, 0)
          t2 <- max(X2 - ch*DF, bX2 - cb*bDF, 0)
          if (is.na(t1) || is.na(t2)) {
            out["cfi.robust"] <- NA
          } else if (t1 == 0 && t2 == 0) {
            out["cfi.robust"] <- 1
          } else out["cfi.robust"] <- 1 - t1/t2
        }
      }
    }
    if ("rni" %in% indices) {
      t1 <- X2 - DF
      t2 <- bX2 - bDF
      out["rni"] <- if (t2 == 0) NA else 1 - t1/t2
      if (robust) {
        ## scaled
        t1 <- X2.sc - DF.sc
        t2 <- bX2.sc - bDF.sc
        if (is.na(t1) || is.na(t2)) {
          out["rni.scaled"] <- NA
        } else if (t2 == 0) {
          out["rni.scaled"] <- NA
        } else out["rni.scaled"] <- 1 - t1/t2
        ## Brosseau-Liard & Savalei MBR 2014, equation 15
        if (!pool.robust & test.names[1] %in%
            c("satorra.bentler","yuan.bentler","yuan.bentler.mplus")) {
          t1 <- X2 - ch*DF
          t2 <- bX2 - cb*bDF
          if (is.na(t1) || is.na(t2)) {
            out["rni.robust"] <- NA
          } else if (t1 == 0 && t2 == 0) {
            out["rni.robust"] <- NA
          } else out["rni.robust"] <- 1 - t1/t2
        }
      }
    }
    if (any(indices %in% c("tli","nnfi"))) {
      t1 <- (X2 - DF)*bDF
      t2 <- (bX2 - bDF)*DF
      out["tli"] <- out["nnfi"] <- if (DF > 0) 1 - t1/t2 else 1
      if (robust) {
        ## scaled
        t1 <- (X2.sc - DF.sc)*bDF.sc
        t2 <- (bX2.sc - bDF.sc)*DF.sc
        if (is.na(t1) || is.na(t2)) {
          out["tli.scaled"] <- out["nnfi.scaled"] <- NA
        } else if (DF > 0 && t2 != 0) {
          out["tli.scaled"] <- out["nnfi.scaled"] <- 1 - t1/t2
        } else {
          out["tli.scaled"] <- out["nnfi.scaled"] <- 1
        }
        ## Brosseau-Liard & Savalei MBR 2014, equation 15
        if (!pool.robust & test.names[1] %in%
            c("satorra.bentler","yuan.bentler","yuan.bentler.mplus")) {
          t1 <- (X2 - ch*DF)*bDF
          t2 <- (bX2 - cb*bDF)*DF
          if (is.na(t1) || is.na(t2)) {
            out["tli.robust"] <- out["nnfi.robust"] <- NA
          } else if (t1 == 0 && t2 == 0) {
            out["tli.robust"] <- out["nnfi.robust"] <- 1 - t1/t2
          } else out["tli.robust"] <- out["nnfi.robust"] <- 1
        }
      }
    }
    if ("rfi" %in% indices) {
      if (DF > 0) {
        t2 <- bX2 / bDF
        t1 <- t2 - X2/DF
        out["rfi"] <- if (t1 < 0 || t2 < 0) 1 else t1/t2
      } else out["rfi"] <- 1
      if (robust) {
        if (DF > 0) {
          t2 <- bX2.sc / bDF.sc
          t1 <- t2 - X2.sc/DF.sc
          out["rfi.scaled"] <- if (t1 < 0 || t2 < 0) 1 else t1/t2
        } else out["rfi.scaled"] <- 1
      }
    }
    if ("nfi" %in% indices) {
      if (DF > 0) {
        t1 <- bX2 - X2
        t2 <- bX2
        out["nfi"] <- t1 / t2
      } else out["nfi"] <- 1
      if (robust) out["nfi.scaled"] <- (bX2.sc - X2.sc) / bX2.sc
    }
    if ("pnfi" %in% indices) {
      t1 <- bX2 - X2
      t2 <- bX2
      out["pnfi"] <- (DF / bDF) * t1/t2
      if (robust) {
        t1 <- bX2.sc - X2.sc
        t2 <- bX2.sc
        out["pnfi.scaled"] <- (DF / bDF) * t1/t2
      }
    }
    if ("ifi" %in% indices) {
      t1 <- bX2 - X2
      t2 <- bX2 - DF
      out["ifi"] <- if (t2 < 0) 1 else t1/t2
      if (robust) {
        t1 <- bX2.sc - X2.sc
        t2 <- bX2.sc - DF.sc
        if (is.na(t2)) {
          out["ifi.scaled"] <- NA
        } else if (t2 < 0) {
          out["ifi.scaled"] <- 1
        } else out["ifi.scaled"] <- t1/t2
      }
    }

    N <- lavListInspect(object, "ntotal")
    Ns <- lavListInspect(object, "nobs") # N per group
    nG <- lavListInspect(object, "ngroups")
    nVars <- length(lavNames(object))
    if (!(lavoptions$likelihood == "normal" |
          lavoptions$estimator %in% c("ML","PML","FML"))) {
      N <- N - nG
      Ns <- Ns - 1
    }

    if ("mfi" %in% indices) {
      out["mfi"] <- exp(-0.5 * (X2 - DF) / N)
    }

    if ("rmsea" %in% indices) {
      N.RMSEA <- max(N, X2*4) # FIXME: good strategy??

      if (is.na(X2) || is.na(DF)) {
        out["rmsea"] <- as.numeric(NA)
      } else if (DF > 0) {
        getLambda <- function(lambda, chi, df, p) pchisq(chi, df, ncp=lambda) - p

        out["rmsea"] <- sqrt( max(0, (X2/N)/DF - 1/N) ) * sqrt(nG)
        ## lower confidence limit
        if (getLambda(0, X2, DF, .95) < 0.0) out["rmsea.ci.lower"] <- 0 else {
          lambda.l <- try(uniroot(f = getLambda, chi = X2, df = DF, p = .95,
                                  lower = 0, upper = X2)$root, silent = TRUE)
          if (inherits(lambda.l, "try-error")) lambda.l <- NA
          out["rmsea.ci.lower"] <- sqrt( lambda.l/(N*DF) ) * sqrt(nG)
        }
        ## upper confidence limit
        if (getLambda(N.RMSEA, X2, DF, .05) > 0 || getLambda(0, X2, DF, .05) < 0) {
          out["rmsea.ci.upper"] <- 0
        } else {
          lambda.u <- try(uniroot(f = getLambda, chi = X2, df = DF, p = .05,
                                  lower = 0, upper = N.RMSEA)$root, silent = TRUE)
          if (inherits(lambda.u, "try-error")) lambda.u <- NA
          out["rmsea.ci.upper"] <- sqrt( lambda.u/(N*DF) ) * sqrt(nG)
        }
        ## p value
        out["rmsea.pvalue"] <- pchisq(X2, DF, ncp = N*DF*0.05^2/nG,
                                      lower.tail = FALSE)

        ## Scaled versions (naive and robust)
        if (robust & !scaleshift) {
          ## naive
          out["rmsea.scaled"] <- sqrt( max(0, (X2/N)/d - 1/N) ) * sqrt(nG)
          ## lower confidence limit
          if (DF < 1 || d < 1 || getLambda(0, X2, d, .95) < 0.0) {
            out["rmsea.ci.lower.scaled"] <- 0
          } else {
            lambda.l <- try(uniroot(f = getLambda, chi = X2, df = d, p = .95,
                                    lower = 0, upper = X2)$root, silent = TRUE)
            if (inherits(lambda.l, "try-error")) lambda.l <- NA
            out["rmsea.ci.lower.scaled"] <- sqrt( lambda.l/(N*DF) ) * sqrt(nG)
          }
          ## upper confidence limit
          if (DF < 1|| d < 1 || getLambda(0, X2, d, .95) < 0.0 || getLambda(N.RMSEA, X2, d, .05) > 0.0) {
            out["rmsea.ci.upper.scaled"] <- 0
          } else {
            lambda.u <- try(uniroot(f = getLambda, chi = X2, df = d, p = .05,
                                    lower = 0, upper = N.RMSEA)$root, silent = TRUE)
            if (inherits(lambda.u, "try-error")) lambda.u <- NA
            out["rmsea.ci.upper.scaled"] <- sqrt( lambda.u/(N*DF) ) * sqrt(nG)
          }
          ## p value
          out["rmsea.pvalue.scaled"] <- pchisq(X2, d, ncp = N*d*0.05^2/nG,
                                               lower.tail = FALSE)

          if (!pool.robust & test.names[1] %in%
              c("satorra.bentler","yuan.bentler","yuan.bentler.mplus")) {
            ## robust
            out["rmsea.robust"] <- sqrt( max(0, (X2/N)/DF - ch/N ) ) * sqrt(nG)
            ## lower confidence limit
            if (DF.sc < 1 | getLambda(0, X2.sc, DF.sc, .95) < 0.0) {
              out["rmsea.ci.lower.robust"] <- 0
            } else {
              lambda.l <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .95,
                                      lower = 0, upper = X2.sc)$root, silent = TRUE)
              if (inherits(lambda.l, "try-error")) lambda.l <- NA
              out["rmsea.ci.lower.robust"] <- sqrt( (ch*lambda.l)/(N*DF.sc) ) * sqrt(nG)
            }
            ## upper confidence limit
            if (DF.sc < 1 | getLambda(N.RMSEA, X2.sc, DF.sc, .05) > 0.0) {
              out["rmsea.ci.upper.robust"] <- 0
            } else {
              lambda.u <- try(uniroot(f = getLambda, chi = X2.sc, df = DF.sc, p = .05,
                                      lower = 0, upper = N.RMSEA)$root, silent = TRUE)
              if (inherits(lambda.u, "try-error")) lambda.u <- NA
              out["rmsea.ci.upper.robust"] <- sqrt( (ch*lambda.u)/(N*DF.sc) ) * sqrt(nG)
            }
            ## p value
            ########## To be discovered?
          }
        } else if (robust & scaleshift) {
          ## naive only
          out["rmsea.scaled"] <- sqrt( max(0, (X2.sc/N)/DF - 1/N) ) * sqrt(nG)
          ## lower confidence limit
          if (DF < 1 | getLambda(0, X2.sc, DF, .95) < 0.0) {
            out["rmsea.ci.lower.scaled"] <- 0
          } else {
            lambda.l <- try(uniroot(f = getLambda, chi = X2.sc, df = DF, p = .95,
                                    lower = 0, upper = X2.sc)$root, silent = TRUE)
            if (inherits(lambda.l, "try-error")) lambda.l <- NA
            out["rmsea.ci.lower.scaled"] <- sqrt( lambda.l/(N*DF) ) * sqrt(nG)
          }
          ## upper confidence limit
          if (DF < 1 | getLambda(N.RMSEA, X2.sc, DF, .05) > 0.0) {
            out["rmsea.ci.upper.scaled"] <- 0
          } else {
            lambda.u <- try(uniroot(f = getLambda, chi = X2.sc, df = DF, p = .05,
                                    lower = 0, upper = N.RMSEA)$root, silent = TRUE)
            if (inherits(lambda.u, "try-error")) lambda.u <- NA
            out["rmsea.ci.upper.scaled"] <- sqrt( lambda.u/(N*DF) ) * sqrt(nG)
          }
          ## p value
          out["rmsea.pvalue.scaled"] <- pchisq(X2.sc, DF, ncp = N*DF*0.05^2/nG,
                                               lower.tail = FALSE)
        }
      }
    }

    if ("gammahat" %in% indices) {
      out["gammaHat"] <- nVars / (nVars + 2*((X2 - DF) / N))
      out["adjGammaHat"] <- 1 - (((nG * nVars * (nVars + 1)) / 2) / DF) * (1 - out["gammaHat"])
      if (robust) {
        out["gammaHat.scaled"] <- nVars / (nVars + 2*((X2.sc - DF.sc) / N))
        out["adjGammaHat.scaled"] <- 1 - (((nG * nVars * (nVars + 1)) / 2) / DF.sc) * (1 - out["gammaHat.scaled"])
      }
    }

    ## END CHI-SQUARED-BASED FIT INDICES
  } else out <- numeric(0)


  ## RESIDUALS-BASED FIT INDICES

  if (any(grepl(pattern = "rmr", x = indices))) {
    if (lavListInspect(object, "nlevels") > 1L) {
      out["srmr"] <- NA # to preserve the order in lavaan output
      out["srmr_within"] <- getSRMR(object, type = "cor", include.mean = FALSE,
                                    level = "within", omit.imps = omit.imps)
      out["srmr_between"] <- getSRMR(object, type = "cor", include.mean = FALSE,
                                     level = lavListInspect(object, "cluster"),
                                     omit.imps = omit.imps)
      out["srmr"] <- out["srmr_within"] + out["srmr_between"]
    } else {
      out["rmr"] <- getSRMR(object, type = "raw", include.mean = TRUE,
                            omit.imps = omit.imps)
      out["rmr_nomean"] <- getSRMR(object, type = "raw", include.mean = FALSE,
                                   omit.imps = omit.imps)
      out["srmr_bentler"] <- out["srmr"] <- getSRMR(object, type = "cor.bentler",
                                                    include.mean = TRUE,
                                                    omit.imps = omit.imps)
      out["srmr_bentler_nomean"] <- getSRMR(object, type = "cor.bentler",
                                            include.mean = FALSE,
                                            omit.imps = omit.imps)
      out["crmr"] <- getSRMR(object, type = "cor.bollen", include.mean = TRUE,
                             omit.imps = omit.imps)
      out["crmr_nomean"] <- getSRMR(object, type = "cor.bollen",
                                    include.mean = FALSE, omit.imps = omit.imps)
      out["srmr_mplus"] <- getSRMR(object, type = "mplus", include.mean = TRUE,
                                   omit.imps = omit.imps)
      out["srmr_mplus_nomean"] <- getSRMR(object, type = "mplus",
                                          include.mean = FALSE,
                                          omit.imps = omit.imps)
    }
    ## END RESIDUALS-BASED FIT INDICES
  }


  ## return requested measures (loosely matched)
  if ("all" %in% fit.measures) {
    fits <- out
  } else {
    fits <- out[grepl(pattern = paste(fit.measures, collapse = "|"),
                      x = names(out), ignore.case = TRUE)]
    fits <- fits[which(!is.na(names(fits)))]
  }
  class(fits) <- c("lavaan.vector","numeric")
  if (output %in% c("text","pretty")) {
    class(fits) <- c("lavaan.fitMeasures", class(fits))
  } else if (output == "matrix") {
    fits <- matrix(fits, dimnames = list(names(fits), ""))
    class(fits) <- c("lavaan.matrix","matrix")
  }
  fits
}


