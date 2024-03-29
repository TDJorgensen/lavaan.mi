### Terrence D. Jorgensen
### Last updated: 29 March 2024
### pool (un)standardized parameters
### analogs of parameterEstimates() and standardizedSolution()


#TODO: add FMI and RIV descriptions to help page:
# messRIV <- paste("The RIV will exceed 1 whenever between-imputation",
#                  "variance exceeds within-imputation variance",
#                  "(when FMI(1) > 50%).\n\n")
#      RIV == FMI / (1 - FMI)

## define with(out) camelCase, for consistency with lavaan
parameterestimates.mi <-
parameterEstimates.mi <- function(object,
                                  # select columns
                                  se = TRUE, zstat = se, pvalue = se,
                                  ci = TRUE, level = .95, fmi = FALSE,
                                  standardized = FALSE,
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
  PT <- parTable(object)
  myCols <- c("lhs","op","rhs","exo")
  if (lavListInspect(object, "ngroups") > 1L) myCols <- c(myCols,"block","group")
  if (lavListInspect(object, "nlevels") > 1L) myCols <- c(myCols,"block","level")
  PE <- PT[ , unique(myCols)]
  free <- PT$free > 0L | PT$op == ":="
  STDs <- !(PT$op %in% c("==","<",">")) # which rows can be standardized

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
      if (zstat)  PE$z[free] <- PE$est[free] / PE$se[free]
      if (pvalue) PE$pvalue  <- pnorm(-abs(PE$z))*2
      if (ci)     crit <- qnorm(1 - (1 - level) / 2)
    } else {
      if (zstat)  PE$t[free] <- PE$est[free] / PE$se[free]
      ## calculate df for t test
      ## can't do finite-sample correction because Wald z tests have no df
      ## (see Enders, 2010, p. 231, eq. 8.13 & 8.14)
      if (zstat)  PE$df[free] <- (m - 1) * (1 + W[free] / Bm[free])^2
      if (pvalue) PE$pvalue <- pt(-abs(PE$t), df = PE$df)*2
      if (ci)     crit <- qt(1 - (1 - level) / 2, df = PE$df)

      ## Set obscenely large DF to infinity for prettier printing?
      if (output == "text" && any(PE$df > 999.5)) {
        PE$df <- ifelse(PE$df > 999, Inf, PE$df)
        infDF <- TRUE # to warn users who would panic
      }
    }
    if (ci) {
      PE$ci.lower <- PE$est - crit * PE$se
      PE$ci.upper <- PE$est + crit * PE$se
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

  if (length(standardized) || rsquare) {
    ## pooled estimates for standardizedSolution()
    est <- coef_lavaan_mi(object, omit.imps = omit.imps)
    ## updates @Model@GLIST for standardizedSolution(..., GLIST=)
    object@Model <- lavaan::lav_model_set_parameters(object@Model, x = est)
  }

  if ("std.lv" %in% standardized) {
    PE$std.lv[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                    type = "std.lv",
                                                    GLIST = object@Model@GLIST,
                                                    est = PE$est)$est.std
  }
  if ("std.all" %in% standardized) {
    PE$std.all[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                     type = "std.all",
                                                     GLIST = object@Model@GLIST,
                                                     est = PE$est)$est.std
  }
  if ("std.nox" %in% standardized) {
    PE$std.nox[STDs] <- lavaan::standardizedSolution(object, se = FALSE,
                                                     type = "std.nox",
                                                     GLIST = object@Model@GLIST,
                                                     est = PE$est)$est.std
  }

  if (fmi) {
    PE$fmi[free] <- Bm[free] / Tot[free]
    PE$riv[free] <- Bm[free] / W[free] # (Enders, 2010, p. 226, eq. 8.10)
  }
  ## fancy or not?
  if (output == "text") {
    PE$label <- PT$label
    #FIXME: no longer needed?  PE$exo <- 0L
    class(PE) <- c("lavaan.parameterEstimates","lavaan.data.frame","data.frame")
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
  } else {
    PE$exo <- NULL
    class(PE) <- c("lavaan.data.frame","data.frame")

    messPool <- paste0("Rubin's (1987) rules were used to pool point",
                       if (se) " and SE",
                       " estimates across ", m, " imputed data sets",
                       if (se & !asymptotic) ", and to calculate degrees of",
                       if (se & !asymptotic) " freedom for each parameter's t",
                       if (se & !asymptotic) " test and CI.",
                       "\n")
    attr(PE, "header") <- messPool
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
    STD <- lavaan::standardizedSolution(object, se = FALSE, type = "std.all",
                                        GLIST = object@Model@GLIST, est = PE$est)
    isEndoSTD <- sapply(STD$lhs, function(x) x %in% endoNames)
    std.all <- STD$est.std[STD$lhs == STD$rhs & STD$op == "~~" & isEndoSTD]
    rsqPE$est <- ifelse(std.all < 0, NA, 1 - std.all) # negative variances
    if (output == "text") rsqPE$label <- ""
    PE <- rbind(PE, rsqPE)
  }

  if (output == "data.frame") PE <- PE[!(PE$op %in% c("==","<",">")), ]
  rownames(PE) <- NULL
  PE
}
