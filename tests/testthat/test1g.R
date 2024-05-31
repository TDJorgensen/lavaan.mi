### Terrence D. Jorgensen
### Last updated: 24 May 2024
### unit tests for single-group lavaan.mi object(s)

library(lavaan.mi)


## ---------------
## Prepare example
## ---------------

data(HS20imps) # import a list of 20 imputed data sets

HS.model <- '
visual  =~ x1 + a*x2 + b*x3
textual =~ x4 + x5 + x6
speed   =~ x7 + x8 + x9
ab := a*b
'
## fit single-group model
fit1 <- cfa.mi(HS.model, data = HS20imps, std.lv = TRUE, meanstructure = TRUE)
fit0 <- cfa.mi(HS.model, data = fit1    , std.lv = TRUE, meanstructure = TRUE,
               #se = "none",
               #estimator = "mlm",
               constraints = '.p2. == .p3. ; .p5. == .p6. ; .p8. == .p9.')


## ----------------------
## test class and methods
## ----------------------

test_that("objects returned by lavaan.mi are class lavaan.mi", {
  expect_s4_class(fit0, "lavaan.mi")
  expect_s4_class(fit1, "lavaan.mi")
})

test_that("summary returns object inheriting from lavaan.summary, with components", {
  sum1 <- summary(fit1, stand = "std.all", rsq = TRUE, fit.measures = TRUE)
  expect_true(inherits(sum1    , "lavaan.summary"))
  expect_true(inherits(sum1$fit, "lavaan.vector"))
  expect_true(inherits(sum1$pe , "data.frame"))
  expect_true(inherits(standardizedSolution.mi(fit1), "data.frame"))
})

test_that("vectors are returned with(out) names", {
  expect_named(coef(fit1))
  expect_named(vcov(fit1, type = "ariv"), expected = NULL)
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
  expect_identical(inherits(fitted(fit1)$cov,
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(fitted.values(fit1)$cov,
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(residuals(fit1)$cov,
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(resid(fit1, type = "cor.bentler")$cov,
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(lavResiduals.mi(fit1)$cov,
                            c("lavaan.matrix.symmetric","matrix"),
                            which = TRUE), expected = 1:2)
})

test_that("lavTestLRT returns named vectors", {
  skip_on_cran()
  expect_named(lavTestLRT.mi(fit1))                     # pooled LRT by D4 (default)
  expect_named(lavTestLRT.mi(fit1, pool.method = "D3")) # pooled LRT by D3
  expect_named(lavTestLRT.mi(fit1, pool.method = "D2")) # pooled Wald test (D2)
  expect_named(lavTestLRT.mi(fit1, asymptotic = TRUE))  # as chisq == F * df1
  expect_named(lavTestLRT.mi(fit0, h1 = fit1, asymptotic = TRUE))  # compare fit
  expect_named(lavTestLRT.mi(fit0, h1 = fit1, pool.method = "D2")) # using D2
})


test_that("score tests return lavaan data.frames", {
  skip_on_cran()
  sc0 <- lavTestScore.mi(fit0, epc = TRUE, cumulative = TRUE)
  expect_identical(inherits(sc0$test, c("lavaan.data.frame","data.frame"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(sc0$uni, c("lavaan.data.frame","data.frame"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(sc0$cumulative, c("lavaan.data.frame","data.frame"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(sc0$epc, c("lavaan.data.frame","data.frame"),
                            which = TRUE), expected = 1:2)
  expect_identical(inherits(modindices.mi(fit1),
                            c("lavaan.data.frame","data.frame"), which = TRUE),
                   expected = 1:2)
  expect_identical(inherits(lavTestScore.mi(fit0, pool.method = "D1")$test,
                            c("lavaan.data.frame","data.frame"), which = TRUE),
                   expected = 1:2)
})



