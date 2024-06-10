### Terrence D. Jorgensen
### Last updated: 10 June 2024
### unit tests for NO between-level variance (should match lavaan)

library(lavaan.mi)


## ---------------
## Prepare example
## ---------------

data(HolzingerSwineford1939) # run lavaan::lavTestLRT() example
HS1 <- HolzingerSwineford1939
## create a list of 5 copies, imitating 5 imputations that don't vary
HS5 <- list(HolzingerSwineford1939,
            HolzingerSwineford1939,
            HolzingerSwineford1939,
            HolzingerSwineford1939,
            HolzingerSwineford1939)

TESTS <- c("browne.residual.adf","scaled.shifted","mean.var.adjusted",
           "satorra.bentler", "yuan.bentler", "yuan.bentler.mplus")

## fit models to complete data
HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '
fit1 <- cfa(HS.model, data = HS1, test = TESTS)
fit0 <- cfa(HS.model, data = HS1, test = TESTS, orthogonal = TRUE)

## fit same models to all "imputations"
imp1 <- cfa.mi(HS.model, data = HS5, test = TESTS)
imp0 <- cfa.mi(HS.model, data = HS5, test = TESTS, orthogonal = TRUE)

## fit same models to pre-"pooled" moments
pre <- poolSat(lapply(HS5, "[", i = paste0("x", 1:9)), meanstructure = TRUE)
pre1 <- cfa(HS.model, data = pre)
pre0 <- cfa(HS.model, data = pre, orthogonal = TRUE)


## ----------------------
## test class and methods
## ----------------------

test_that("pooled point/SE estimates match complete data", {
  sumC <- parameterEstimates(   fit1, stand = "std.all")
  sumM <- parameterEstimates.mi(imp1, stand = "std.all", asymptotic = TRUE)

  expect_true(isTRUE(all.equal(sumC$est    , sumM$est    )))
  expect_true(isTRUE(all.equal(sumC$se     , sumM$se     )))
  expect_true(isTRUE(all.equal(sumC$std.all, sumM$std.all)))

  stdC <- standardizedSolution(   fit0)
  stdM <- standardizedSolution.mi(imp0)

  expect_true(isTRUE(all.equal(stdC$est.std, unclass(stdM$est.std))))
  expect_true(isTRUE(all.equal(stdC$se     , stdM$se     )))
})

test_that("pooled test stats match complete data", {
  expect_true(isTRUE(all.equal(sumC$z, sumM$z)))
  expect_true(isTRUE(all.equal(stdC$z, unclass(stdM$z))))

  testC <- lavTestLRT(   fit1, fit0, method = "standard")
  testM <- lavTestLRT.mi(imp1, imp0, method = "standard", asymptotic = TRUE)

  expect_true(isTRUE(all.equal(testC[,"Df"], testM[,"Df"])))
  expect_true(isTRUE(all.equal(testC[,"Chisq"], testM[,"Chisq"])))
  expect_true(isTRUE(all.equal(testC[2,"Chisq diff"], testM[2,"Chisq diff"])))

  testCb <- lavTestLRT(   fit1, fit0, type = "browne.residual.adf")
  testMb <- lavTestLRT.mi(imp1, imp0, type = "browne.residual.adf", asymptotic = TRUE)

  expect_true(isTRUE(all.equal(testCb[ ,"Df"        ], testMb[ ,"Df"        ])))
  expect_true(isTRUE(all.equal(testCb[ ,"Chisq"     ], testMb[ ,"Chisq"     ])))
  expect_true(isTRUE(all.equal(testCb[2,"Chisq diff"], testMb[2,"Chisq diff"])))
})

test_that("pre-pooled analysis matches complete data", {
  expect_true(isTRUE(all.equal(unclass(pre$sample.cov),
                               cov.wt(HS1[paste0("x",1:9)], method = "ML")$cov,
                               tolerance = .0001)))
  expect_true(isTRUE(all.equal(unclass(pre$sample.mean),
                               colMeans(HS1[paste0("x",1:9)]))))
  expect_true(isTRUE(all.equal(pre$NACOV,
                               lavInspect(pre1, "gamma"))))

  testC <- lavTestLRT(fit1, fit0, method = "standard")
  testM <- lavTestLRT(pre1, pre0, method = "standard", asymptotic = TRUE)

  expect_true(isTRUE(all.equal(testC[ ,"Df"        ], testM[ ,"Df"        ]                      )))
  expect_true(isTRUE(all.equal(testC[ ,"Chisq"     ], testM[ ,"Chisq"     ], tolerance = .0000001)))
  expect_true(isTRUE(all.equal(testC[2,"Chisq diff"], testM[2,"Chisq diff"], tolerance = .0000001)))

  testCb <- lavTestLRT(fit1, fit0, type = "browne.residual.adf")
  testMb <- lavTestLRT(pre1, pre0, type = "browne.residual.adf", asymptotic = TRUE)

  expect_true(isTRUE(all.equal(testCb[ ,"Df"        ], testMb[ ,"Df"        ]                      )))
  expect_true(isTRUE(all.equal(testCb[ ,"Chisq"     ], testMb[ ,"Chisq"     ], tolerance = .0000001)))
  expect_true(isTRUE(all.equal(testCb[2,"Chisq diff"], testMb[2,"Chisq diff"], tolerance = .0000001)))
})

