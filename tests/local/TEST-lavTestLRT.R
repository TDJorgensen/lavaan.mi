### Terrence D. Jorgensen
### Last updated: 10 June 2024
### verify that lavTestLRT.mi() can obtain any requested test from the call
#TODO: turn some of this into a set of /tests

data(HS20imps) # import a list of 20 imputed data sets

## specify CFA model from lavaan's ?cfa help page
HS.model <- '
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9
'

## lavaan.mi objects
fit1s  <- cfa.mi(HS.model, data = HS20imps) # standard stat only
fit1b  <- cfa.mi(HS.model, data = HS20imps, test = "browne.residual.adf")
fit1y  <- cfa.mi(HS.model, data = HS20imps, # mean-adjusted stat
                 test = "yuan.bentler")
fit1p  <- cfa.mi(HS.model, data = HS20imps, # mean-adjusted stat
                 test = "yuan.bentler.mplus")
fit1m  <- cfa.mi(HS.model, data = HS20imps, # mean-adjusted stat
                 estimator = "mlm")
fit1mv <- cfa.mi(HS.model, data = HS20imps, # mean- and variance-adjusted stat
                 estimator = "mlmvs")
fit1sh <- cfa.mi(HS.model, data = HS20imps, # scaled and shifted stat
                 estimator = "mlmv")

## orthogonality
fit0s  <- cfa.mi(HS.model, data = HS20imps, orthogonal = TRUE) # standard stat only
fit0b  <- cfa.mi(HS.model, data = HS20imps, orthogonal = TRUE,
                 test = "browne.residual.adf")
fit0y  <- cfa.mi(HS.model, data = HS20imps, # mean-adjusted stat
                 test = "yuan.bentler", orthogonal = TRUE)
fit0p  <- cfa.mi(HS.model, data = HS20imps, # mean-adjusted stat
                 test = "yuan.bentler.mplus", orthogonal = TRUE)
fit0m  <- cfa.mi(HS.model, data = HS20imps, # mean-adjusted stat
                 estimator = "mlm", orthogonal = TRUE)
fit0mv <- cfa.mi(HS.model, data = HS20imps, # mean- and variance-adjusted stat
                 estimator = "mlmvs", orthogonal = TRUE)
fit0sh <- cfa.mi(HS.model, data = HS20imps, # scaled and shifted stat
                 estimator = "mlmv", orthogonal = TRUE)


## check whether ?lavTestLRT example works
## (multiple tests in a single fit)
t6.1 <- cfa.mi(HS.model, data = HS20imps,
               test = c("browne.residual.adf","scaled.shifted","mean.var.adjusted",
                        "satorra.bentler", "yuan.bentler", "yuan.bentler.mplus"))
t6.0 <- cfa.mi(HS.model, data = HS20imps, orthogonal = TRUE,
               test = c("browne.residual.adf","scaled.shifted","mean.var.adjusted",
                        "satorra.bentler", "yuan.bentler", "yuan.bentler.mplus"))


# By default (test="default", type="Chisq"), the first scaled statistic
## requested will be used. Here, that is "scaled.shifted"
lavTestLRT.mi(t6.1, t6.0)
## But even if "satorra.bentler" were requested first, method="satorra.2000"
## provides the scaled-shifted chi-squared difference test:
lavTestLRT.mi(t6.1, t6.0, method = "satorra.2000")
## == lavTestLRT.mi(fit0sh, fit1sh)

## The mean- and variance-adjusted (Satterthwaite) statistic implies
## scaled.shifted = FALSE
lavTestLRT.mi(t6.1, t6.0, method = "satorra.2000", scaled.shifted = FALSE)
## == lavTestLRT.mi(fit0mv, fit1mv)

## Because "satorra.bentler" is not the first scaled test in the list,
## we MUST request it explicitly:
lavTestLRT.mi(t6.1, t6.0, test = "satorra.bentler") # method="satorra.bentler.2001"
## == lavTestLRT.mi(fit0m, fit1m)
## The "strictly-positive test" is necessary when the above test is < 0:
lavTestLRT.mi(t6.1, t6.0, test = "satorra.bentler", method = "satorra.bentler.2010")

## Likewise, other scaled statistics can be selected:
lavTestLRT.mi(t6.1, t6.0, test = "yuan.bentler")
## == lavTestLRT.mi(fit0y, fit1y)
lavTestLRT.mi(t6.1, t6.0, test = "yuan.bentler.mplus")
## == lavTestLRT.mi(fit0p, fit1p)

## To request the difference between Browne's (1984) residual-based statistics,
## rather than statistics based on the fitted model's discrepancy function,
## use the type= argument:
lavTestLRT.mi(t6.1, t6.0, type = "browne.residual.adf")
## == lavTestLRT.mi(fit0b, fit1b, type = "browne.residual.adf")

## Despite requesting multiple robust tests, it is still possible to obtain
## the standard chi-squared difference test (i.e., without a robust correction)
lavTestLRT.mi(t6.1, t6.0, method = "standard")
## == lavTestLRT.mi(fit0s, fit1s)





