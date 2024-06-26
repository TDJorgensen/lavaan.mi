### Terrence D. Jorgensen
### Last updated: 24 May 2024
### unit tests for multigroup lavaan.mi object(s)

library(lavaan.mi)


## ---------------
## Prepare example
## ---------------

data(HS20imps) # import a list of 20 imputed data sets

HS.model <- '
visual  =~ x1 + x2 + x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
'
## fit single-group model
fit1 <- cfa.mi(HS.model, data = HS20imps, estimator = "mlmv", group = "school")
fit0 <- cfa.mi(HS.model, data = HS20imps, estimator = "mlmv", group = "school",
               group.equal = c("loadings","intercepts"))


## ----------------------
## test class and methods
## ----------------------

test_that("objects returned by lavaan.mi are class lavaan.mi", {
  expect_s4_class(fit0, "lavaan.mi")
  expect_s4_class(fit1, "lavaan.mi")
})

test_that("summary returns object inheriting from lavaan.summary, with components", {
  sum1 <- summary(fit1, stand = "std.all", asymptotic = TRUE, fmi = TRUE, ci = TRUE)
  expect_true(inherits(sum1   , "lavaan.summary"))
  expect_true(inherits(sum1$pe, "data.frame"))
})

test_that("vectors are returned with(out) names", {
  expect_named(coef(fit1))
  expect_named(vcov(fit1, type = "ariv"), expected = NULL)
  expect_named(nobs(fit1, total = FALSE))
  expect_named(nobs(fit1), expected = NULL)
  expect_named(fitMeasures(fit1))
})

test_that("vcov returns lavaan matrices", {
  expect_identical(inherits(vcov(fit1),
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(vcov(fit1, type = "within"),
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(vcov(fit1, type = "between"),
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(fitMeasures(fit1, output = "matrix"),
                            c("lavaan.matrix","matrix"),
                            which = TRUE), expected = 1:2)
})

test_that("fitted and resid return lavaan matrices in $cov", {
  expect_identical(inherits(fitted(fit1)$Pasteur$cov,
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(fitted.values(fit1)$Pasteur$cov,
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(residuals(fit1)$Pasteur$cov,
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(resid(fit1, type = "cor.bentler")$Pasteur$cov,
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(lavResiduals.mi(fit1)$Pasteur$cov,
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
})

test_that("lavTestLRT and Wald return named vectors", {
  skip_on_cran()
  waldcon <- ' .p70. == 0 ; .p71. == 0 ; .p72. == 0 '
  expect_named(lavTestLRT.mi(fit1, asymptotic = TRUE, pool.method = "D2"))
  expect_named(lavTestLRT.mi(fit1, asymptotic = TRUE, pool.method = "D2",
                             pool.robust = TRUE))
  expect_named(lavTestLRT.mi(fit0, h1 = fit1, asymptotic = TRUE, pool.method = "D2"))
  expect_named(lavTestLRT.mi(fit0, h1 = fit1, asymptotic = TRUE, pool.method = "D2",
                             pool.robust = TRUE))
  expect_named(lavTestWald.mi(fit0, pool.method = "D1", asymptotic = TRUE, constraints = waldcon))
  expect_named(lavTestWald.mi(fit0, pool.method = "D2", asymptotic = TRUE, constraints = waldcon))
})

test_that("score tests match modification indices", {
  skip_on_cran()
  suppressWarnings(
    ## because I used a robust test statistic
    sc1 <- lavTestScore.mi(fit1, add = 'x1 ~~ x9 ; x4 ~~ x7', pool.method = "D1",
                           asymptotic = TRUE, epc = TRUE)
  )
  suppressWarnings(mi1 <- modindices.mi(fit1, op = '~~', pool.method = "D1"))

  sc2 <- lavTestScore.mi(fit1, add = 'x1 ~~ x9 ; x4 ~~ x7', pool.method = "D2",
                         asymptotic = TRUE, epc = TRUE)
  mi2 <- modindices.mi(fit1, op = '~~', pool.method = "D2")

  expect_equal(mi1$mi[(mi1$lhs == "x1" & mi1$rhs == "x9") | (mi1$lhs == "x4" & mi1$rhs == "x7")],
               sc1$uni$X2, tolerance = .0001)
  expect_equal(mi1$epc[(mi1$lhs == "x1" & mi1$rhs == "x9") | (mi1$lhs == "x4" & mi1$rhs == "x7")],
               sc1$uni$epc, tolerance = .0001)
  expect_equal(mi2$mi[(mi2$lhs == "x1" & mi2$rhs == "x9") | (mi2$lhs == "x4" & mi2$rhs == "x7")],
               sc2$uni$X2, tolerance = .0001)
  expect_equal(mi2$epc[(mi2$lhs == "x1" & mi2$rhs == "x9") | (mi2$lhs == "x4" & mi2$rhs == "x7")],
               sc2$uni$epc, tolerance = .0001)
})



