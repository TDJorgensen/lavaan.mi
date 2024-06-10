### Terrence D. Jorgensen
### Last updated: 6 May 2022
### remaining problems for lavaan.mi methods


## -------------------------
## fixed.x and conditional.x
## -------------------------

## Using help-page example
# data(datCat)
# datCat$u2 <- as.integer(datCat$u2) # mixture of ordered/continuous indicators
# datCat$u3 <- as.integer(datCat$u3)
# datCat$u5 <- as.integer(datCat$u5) # exogenous predictors must be numeric
# datCat$u6 <- as.integer(datCat$u6)
#
# ## impose missing values
# set.seed(456)
# for (i in 1:8) datCat[sample(1:nrow(datCat), size = .1*nrow(datCat)), i] <- NA
# library(Amelia)
# catimps <- amelia(datCat, m = 20, ords = paste0("u", 1:8), noms = "g")$imputations
#
# catmod <- '
# f =~ 1*u1 + 1*u2 + 1*u3 + 1*u4
# u1 + u2 ~ u5 + u6
# '
# fitx <- cfa.mi(catmod, data = catimps, fixed.x = TRUE, conditional.x = FALSE)
# fitxg <- cfa.mi(catmod, data = catimps, fixed.x = TRUE, conditional.x = FALSE,
#                 group = "g")
# fit.x <- cfa.mi(catmod, data = catimps, fixed.x = TRUE, conditional.x = TRUE)
# fit.xg <- cfa.mi(catmod, data = catimps, fixed.x = TRUE, conditional.x = TRUE,
#                  group = "g")
#
# fitted(fitx)
# resid(fitx)
# resid(fitx, type = "crmr")
# resid(fitx, type = "srmr")
# fitMeasures(fitx, fit.measures = "rmr")
#
# fitted(fitxg)
# resid(fitxg)
# resid(fitxg, type = "crmr")
# resid(fitxg, type = "srmr")
# fitMeasures(fitxg, fit.measures = "rmr")
#
# fitted(fit.x)
# resid(fit.x)
# resid(fit.x, type = "crmr")
# resid(fit.x, type = "srmr")
# fitMeasures(fit.x, fit.measures = "rmr")
#
# fitted(fit.xg)
# resid(fit.xg)
# resid(fit.xg, type = "crmr")
# resid(fit.xg, type = "srmr")
# fitMeasures(fit.xg, fit.measures = "rmr")
#
# summary(fitx, stand=TRUE) # only problem left: standardizing requires cov.x
# summary(fitxg, stand=TRUE) # only problem left: standardizing requires cov.x
# summary(fit.x, stand=TRUE) # only problem left: standardizing requires cov.x
# summary(fit.xg, stand=TRUE) # only problem left: standardizing requires cov.x


## ------------------------------------
## check lavaan.mi with multilevel data
## ------------------------------------

# data(Demo.twolevel)
# Demo.twolevel$id <-  paste0("id", Demo.twolevel$cluster) # character IDs
# ## assign clusters to arbitrary groups
# Demo.twolevel$g <- ifelse(Demo.twolevel$cluster %% 2L, "foo", "bar")
# ## randomize rows so cluster IDs are out of order
# set.seed(123)
# Demo.twolevel <- Demo.twolevel[sample(nrow(Demo.twolevel)), ]
# ## create missing data
# set.seed(456)
# Demo.twolevel$y1[sample(nrow(Demo.twolevel), size = .1*nrow(Demo.twolevel)) ] <- NA
#
# model <- '
# level: within
#   fw =~ y1 + L2w*y2 + L3w*y3
#   fw ~ x1 + x2 + x3
# level: between
#   fb =~ y1 + L2b*y2 + L3b*y3
#   fb ~ b1*w1 + b2*w2
#
# ## constraints
# L2w == L2b
# L3w == L3b
# '
#
# model2 <- ' group: foo
# level: within
#   fw =~ y1 + L2*y2 + L3*y3
#   fw ~ x1 + x2 + x3
# level: between
#   fb =~ y1 + L2*y2 + L3*y3
#   #fb ~ w1 + w2
#
# group: bar
#
# level: within
#   fw =~ y1 + L2*y2 + L3*y3
#   #fw ~ x1 + x2 + x3
# level: between
#   fb =~ y1 + L2*y2 + L3*y3
#   fb ~ w1 + w2
# '
#
# ## impute data
# library(mice)
# m <- 5
# mice.out <- mice(Demo.twolevel, m = m, diagnostics = TRUE)
# imputedData <- list()
# for (i in 1:m) {
#   imputedData[[i]] <- complete(data = mice.out, action = i, include = FALSE)
# }
#
# fit.mi  <- sem.mi(model, data = imputedData, cluster = "id")
# fit2.mi <- sem.mi(model2, data = imputedData, group = "g", cluster = "id")
#
# ## can't return (C/R)RMR for multigroup MLSEM
# fitMeasures(fit.mi , output = "pretty")
# fitMeasures(fit2.mi, output = "pretty")
# fitMeasures(fit2.mi, output = "pretty", # no problem for other indices
#             fit.measures = c("chisq","df","pvalue","cfi","tli","rmsea"))




