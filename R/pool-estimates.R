### Terrence D. Jorgensen
### Last updated: 24 April 2024
### pool (un)standardized parameters
### analogs of parameterEstimates() and standardizedSolution()


#TODO: add FMI and RIV descriptions to help page:
# messRIV <- paste("The RIV will exceed 1 whenever between-imputation",
#                  "variance exceeds within-imputation variance",
#                  "(when FMI(1) > 50%).\n\n")
#      RIV == FMI / (1 - FMI)


##' Pooled Parameter Estimates
##'
##' This function pools parameter estimates from a lavaan model fitted to
##' multiple imputed data sets.
##'
##'
##' @aliases parameterEstimates.mi parameterestimates.mi
##' @importFrom lavaan lavNames lavListInspect parTable
parameterEstimates.mi <- function(object,
                                  # select columns
                                  se = TRUE, zstat = se, pvalue = zstat,
                                  ci = TRUE, level = .95, fmi = FALSE,
                                  standardized = FALSE, cov.std = TRUE,
                                  # add rows
                                  rsquare = FALSE,
                                  # control
                                  scale.W = !asymptotic,
                                  omit.imps = c("no.conv","no.se"),
                                  asymptotic = FALSE,
                                  # remove rows
                                  remove.system.eq = TRUE,
                                  remove.eq = TRUE,
                                  remove.ineq = TRUE,
                                  remove.def = FALSE,
                                  remove.nonfree = FALSE,
                                  remove.step1 = TRUE,
                                  remove.unused = FALSE,
                                  # output
                                  output = "data.frame",
                                  header = FALSE) {
  # check output= argument (from lavaan::parameterEstimates)
  output <- tolower(output)
  if (output %in% c("data.frame", "table")) {
    output <- "data.frame"
  } else if (output %in% c("text", "pretty")) {
    output <- "text"
  } else {
    stop(
      "lavaan ERROR: output must be ", sQuote("data.frame"),
      " or ", sQuote("text")
    )
  }

  useImps <- imps2use(object = object, omit.imps = omit.imps)
  m <- length(useImps)

  lavoptions <- lavListInspect(object, "options")

  ## extract parameter table with attributes for printing
  PE <- parTable(object)
  free <- PE$free > 0L | PE$op == ":="
  STDs <- !(PE$op %in% c("==","<",">")) # which rows can be standardized

  PE$est <- coef_lavaan_mi(object, type = "all", omit.imps = omit.imps)

  if (lavoptions$se == "none") {
    warning('pooled variances and tests unavailable when se="none" is requested')
    se <- FALSE
  }
  if (!se) {
    zstat  <- FALSE
    pvalue <- FALSE
    fmi    <- FALSE
  }
  if (!zstat) pvalue <- FALSE

  if (se) {
    VCOV <- vcov_lavaan_mi(object, scale.W = scale.W, omit.imps = omit.imps)
    PE$se <- lavaan::lav_model_vcov_se(object@Model, VCOV = VCOV,
                                       lavpartable = object@ParTable)
    W <- rowMeans(sapply(object@ParTableList[useImps], "[[", i = "se")^2)
    B <- apply(sapply(object@ParTableList[useImps], "[[", i = "est"), 1, var)
    Bm <- B + B/m
    Tot <- W + Bm
    infDF <- FALSE # Need a message about t(Inf)?
    if (asymptotic) {
      ## can't do finite-sample correction because Wald z tests have no df
      ## (see Enders, 2010, p. 231, eq. 8.13 & 8.14)
      if (zstat)  PE$z[free] <- PE$est[free] / PE$se[free]
      if (pvalue) PE$pvalue  <- pnorm(-abs(PE$z))*2
      if (ci)     crit <- qnorm(1 - (1 - level) / 2)
    } else {
      ## calculate df for t test
      if (zstat | ci)  df4t <- (m - 1) * (1 + W[free] / Bm[free])^2
      if (zstat) {
        PE$t[free] <- PE$est[free] / PE$se[free]
        PE$df[free] <- df4t
      }
      if (pvalue) PE$pvalue <- pt(-abs(PE$t), df = PE$df)*2
      if (ci)     crit <- qt(1 - (1 - level) / 2, df = df4t)

      ## Set obscenely large DF to infinity for prettier printing?
      if (zstat && output == "text" && any(PE$df > 999.5)) {
        PE$df <- ifelse(PE$df > 999, Inf, PE$df)
        infDF <- TRUE # to warn users who would panic
      }
    }
    if (ci) {
      PE$ci.lower[free] <- PE$est[free] - crit * PE$se[free]
      PE$ci.upper[free] <- PE$est[free] + crit * PE$se[free]
      PE$ci.lower[!free] <- PE$ci.upper[!free] <- PE$est[!free]
    }
  }

  if (is.logical(standardized)) {
    if (standardized) {
      standardized <- c("std.lv","std.all")
      if (length(lavNames(object, "ov.x")) && lavoptions$fixed.x) {
        standardized <- c(standardized, "std.nox")
      }
    } else standardized <- NULL
  } else {
    # !is.logical(standardized)
    standardized <- tolower(as.character(standardized))
    if ("std.nox" %in% standardized) {
      # sanity checks
      if (length(lavNames(object, "ov.x")) == 0) {
        message("`std.nox' unavailable without fixed exogenous predictors")
        standardized <- setdiff(standardized, "std.nox")
      }
      if (!object@Options$fixed.x) {
        message("`std.nox' unavailable when fixed.x=FALSE")
        standardized <- setdiff(standardized, "std.nox")
      }
    }
  }

  ## When calling standardizedSolution.mi(), cov.std= and other arguments
  ## remain unevaluated, leading to errors that the objects aren't found.
  ## Below I implement each of solutions 1-3 posted here:
  ## https://stackoverflow.com/questions/68473954/understanding-scoping-of-nested-functions

  # Add each requested type
  if ("std.lv" %in% standardized) {
    ## Solution 2: construct a call, which works because eval() sets the
    ## envir= and enclos= to the parent.frame()
    stdCall <- list(quote(standardizedSolution.mi), object = object,
                    omit.imps = omit.imps,
                    type = "std.lv", cov.std = cov.std,
                    se = FALSE, ci = FALSE,
                    remove.eq   = remove.eq,
                    remove.ineq = remove.ineq,
                    remove.def  = remove.def)
    STD.lv <- eval(as.call(stdCall))
    names(STD.lv)[names(STD.lv) == "est.std"] <- "std.lv"
    PE <- merge(PE, STD.lv, all.x = TRUE, sort = FALSE)
  }
  if ("std.all" %in% standardized) {
    ## Solution 1: set standardizedSolution.mi()'s environment to this one
    environment(standardizedSolution.mi) <- environment()
    STD.all <- standardizedSolution.mi(object, omit.imps = omit.imps,
                                       type = "std.all", cov.std = cov.std,
                                       se = FALSE, ci = FALSE,
                                       remove.eq   = remove.eq,
                                       remove.ineq = remove.ineq,
                                       remove.def  = remove.def)
    names(STD.all)[names(STD.all) == "est.std"] <- "std.all"
    PE <- merge(PE, STD.all, all.x = TRUE, sort = FALSE)
  }
  if ("std.nox" %in% standardized) {
    ## Solution 3: define a nested function (equivalent to standardizedSolution.mi)
    std.mi <- function(object, omit.imps = c("no.conv","no.se"), ...) {
      ## make fake lavaan object from pooled results
      FIT <- mi2lavaan(object, std = TRUE, omit.imps = omit.imps)

      ## pass it to lavaan::standardizedSolution()
      MC <- match.call()
      MC[[1]] <- quote(lavaan::standardizedSolution)
      MC$object <- FIT
      MC$omit.imps <- NULL
      eval(MC)
    }
    STD.nox <- std.mi(object, omit.imps = omit.imps,
                      type = "std.nox", cov.std = cov.std,
                      se = FALSE, ci = FALSE,
                      remove.eq   = remove.eq,
                      remove.ineq = remove.ineq,
                      remove.def  = remove.def)
    names(STD.nox)[names(STD.nox) == "est.std"] <- "std.nox"
    PE <- merge(PE, STD.nox, all.x = TRUE, sort = FALSE)
  }

  ## requested R-squared?
  endoNames <- c(lavNames(object, "ov.nox"), lavNames(object, "lv.nox"))
  if (rsquare & length(endoNames)) {
    isEndo <- sapply(PE$lhs, function(x) x %in% endoNames)
    rsqPE <- PE[PE$lhs == PE$rhs & PE$op == "~~" & isEndo, ]
    rsqPE$op <- "r2"
    for (i in which(!sapply(colnames(PE),
                            function(x) x %in% c("lhs","op","rhs","block",
                                                 "level","group","est","exo")))) {
      rsqPE[ , i] <- NA
    }
    ## maybe redundant, but might have extra columns
    environment(standardizedSolution.mi) <- environment()
    STD <- standardizedSolution.mi(object, omit.imps = omit.imps,
                                   type = "std.all", cov.std = cov.std,
                                   se = FALSE, ci = FALSE,
                                   remove.eq   = remove.eq,
                                   remove.ineq = remove.ineq,
                                   remove.def  = remove.def)
    isEndoSTD <- sapply(STD$lhs, function(x) x %in% endoNames)
    std.all <- STD$est.std[STD$lhs == STD$rhs & STD$op == "~~" & isEndoSTD]
    rsqPE$est <- ifelse(std.all < 0, NA, 1 - std.all) # negative variances
    if (output == "text") rsqPE$label <- ""
    PE <- rbind(PE, rsqPE)
  }

  ## fraction of missing information (and relative increase in variance)
  if (fmi) {
    PE$fmi[free] <- Bm[free] / Tot[free]
    PE$riv[free] <- Bm[free] / W[free] # (Enders, 2010, p. 226, eq. 8.10)
  }



  ## FROM HERE ON, heavily borrow / mimic lavaan::parameterEstimates
  ## (see file lav_object_methods.R)
  ## i.e., Find "tmp.list" and Replace with "PE"


  # if single level, remove level column
  if (object@Data@nlevels == 1L) PE$level <- NULL

  # if single group, remove group column
  if (object@Data@ngroups == 1L) PE$group <- NULL

  # if single everything, remove block column
  if (object@Data@nlevels == 1L &&
      object@Data@ngroups == 1L) {
    PE$block <- NULL
  }

  # if no user-defined labels, remove label column
  if (sum(nchar(object@ParTable$label)) == 0L) {
    PE$label <- NULL
  }

  # remove non-free parameters? (but keep ==, >, < and :=)
  if (remove.nonfree) {
    nonfree.idx <- which(PE$free == 0L & !PE$op %in% c("==", ">", "<", ":="))
    if (length(nonfree.idx) > 0L) {
      PE <- PE[-nonfree.idx, ]
    }
  }

  # remove 'unused' parameters
  # these are parameters that are automatically added (user == 0),
  # but with their final (est) values fixed to their default values
  # (typically 1 or 0).
  # currently only intercepts and scaling-factors (for now)
  # should we also remove fixed-to-1 variances? (parameterization = theta)?
  if (remove.unused) {
    # intercepts
    int.idx <- which(PE$op == "~1" &
                     PE$user == 0L &
                     PE$free == 0L &
                     PE$est == 0)
    if (length(int.idx) > 0L) {
      PE <- PE[-int.idx, ]
    }

    # scaling factors
    scaling.idx <- which(PE$op == "~*~" &
                         PE$user == 0L &
                         PE$free == 0L &
                         PE$est == 1)
    if (length(scaling.idx) > 0L) {
      PE <- PE[-scaling.idx, ]
    }
  }


  # remove 'free' column
  PE$free <- NULL

  # remove == rows?
  if (remove.eq) {
    eq.idx <- which(PE$op == "==" & PE$user == 1L)
    if (length(eq.idx) > 0L) {
      PE <- PE[-eq.idx, ]
    }
  }
  if (remove.system.eq) {
    eq.idx <- which(PE$op == "==" & PE$user != 1L)
    if (length(eq.idx) > 0L) {
      PE <- PE[-eq.idx, ]
    }
  }
  # remove <> rows?
  if (remove.ineq) {
    ineq.idx <- which(PE$op %in% c("<", ">"))
    if (length(ineq.idx) > 0L) {
      PE <- PE[-ineq.idx, ]
    }
  }
  # remove := rows?
  if (remove.def) {
    def.idx <- which(PE$op == ":=")
    if (length(def.idx) > 0L) {
      PE <- PE[-def.idx, ]
    }
  }

  # remove step 1 rows?
  ## only relevant for sam(), but there is no sam.mi()
  # if (remove.step1 && !is.null(PE$step)) {
  #   step1.idx <- which(PE$step == 1L)
  #   if (length(step1.idx) > 0L) {
  #     PE <- PE[-step1.idx, ]
  #   }
  #   # remove step column
  #   PE$step <- NULL
  # }

  # remove attribute for data order
  attr(PE, "ovda") <- NULL

  # remove PE$user
  PE$user <- NULL


  ## fancy or not?
  if (output == "text") {
    class(PE) <- c("lavaan.parameterEstimates","lavaan.data.frame","data.frame")

    if (header) {
      attr(PE, "categorical") <- lavoptions$categorical
      attr(PE, "parameterization") <- lavoptions$parameterization
      attr(PE, "information") <- lavoptions$information[1]
      attr(PE, "information.meat") <- lavoptions$information.meat
      attr(PE, "se") <- lavoptions$se
      attr(PE, "group.label") <- lavListInspect(object, "group.label")
      attr(PE, "level.label") <- c("within", lavListInspect(object, "cluster"))
      attr(PE, "bootstrap") <- lavoptions$bootstrap
      attr(PE, "bootstrap.successful") <- 0L #FIXME: assumes none. Implement Wei & Fan's mixing method?
      attr(PE, "missing") <- lavoptions$missing
      attr(PE, "observed.information") <- lavoptions$observed.information[1]
      attr(PE, "h1.information") <- lavoptions$h1.information[1]
      attr(PE, "h1.information.meat") <- lavoptions$h1.information.meat
      attr(PE, "header") <- header
      # FIXME: lavaan may add more!!

      ## add details about pooling options
      attr(PE, "pooled") <- TRUE
      attr(PE, "scale.W") <- scale.W
      attr(PE, "asymptotic") <- asymptotic
      attr(PE, "infDF") <- infDF
    }

  } else {
    PE$exo <- NULL
    PE$lower <- PE$upper <- NULL
    class(PE) <- c("lavaan.data.frame","data.frame")

    messPool <- paste0("Rubin's (1987) rules were used to pool point",
                       if (se) " and SE",
                       " estimates across ", m, " imputed data sets",
              if (zstat | ci & !asymptotic) ", and to calculate degrees of",
              if (zstat | ci & !asymptotic) " freedom for each parameter's",
              if (zstat & !asymptotic) " t test",
              if (zstat & ci & !asymptotic) " and",
              if (ci & !asymptotic) " CI",
              ".\n")
    attr(PE, "header") <- messPool
  }

  rownames(PE) <- NULL
  PE
}



standardizedSolution.mi <- function(object,
                                    omit.imps = c("no.conv","no.se"), ...) {
  ## make fake lavaan object from pooled results
  FIT <- mi2lavaan(object, std = TRUE, omit.imps = omit.imps)

  ## pass it to lavaan::standardizedSolution()
  MC <- match.call()
  MC[[1]] <- quote(lavaan::standardizedSolution)
  MC$object <- FIT
  MC$omit.imps <- NULL
  eval(MC)
}
## define with(out) camelCase, for consistency with lavaan
# standardizedsolution.mi <- standardizedSolution.mi


