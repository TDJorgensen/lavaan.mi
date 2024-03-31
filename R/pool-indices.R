### Terrence D. Jorgensen
### Last updated: 31 March 2024
### pool fit indices
### define fitMeasures() method for lavaan.mi




##' @importFrom lavaan lavNames lavListInspect
##' @importFrom stats pchisq uniroot
fitMeasures_mi <- function(object, fit.measures = "all", baseline.model = NULL,
                           fm.args = list(standard.test     = "default",
                                          scaled.test       = "default",
                                          rmsea.ci.level    = 0.90,
                                          rmsea.close.h0    = 0.05,
                                          rmsea.notclose.h0 = 0.08,
                                          cat.check.pd      = TRUE),
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
    test <- dots$test
    if (is.null(test)) {
      test <- "d4"
    } else test <- tolower(test[1])
    if (tolower(test) %in% c("cm","chan.meng","new.lrt","d4")) test <- "D4"
    if (tolower(test) %in% c("mr","meng.rubin","old.lrt","mplus","d3")) test <- "D3"
    if (tolower(test) %in% c("lmrr","li.et.al","pooled.wald","d2")) test <- "D2"
    if (test %in% c("D3","D4") && !grepl(pattern = "ML", x = lavoptions$estimator)) {
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
    # if (scaleshift) {
    #   if (test %in% c("D3","D4")) {
    #     message("If test = 'scaled.shifted' (estimator = 'WLSMV' or 'MLMV'), ",
    #             "model evaluation is only available by (re)setting ",
    #             "test = 'D2'.\nControl more options by passing arguments to ",
    #             "lavTestLRT() via the '...' argument.\n")
    #     test <- 'D2'
    #   }
    # }

    if (pool.robust && test %in% c("D3","D4")) {
      message('pool.robust = TRUE is only applicable when test = "D2". ',
              'Changed test to "D2".')
      test <- "D2"
    }

    dots$test <- test


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
      } else if (test == "D2") {
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


