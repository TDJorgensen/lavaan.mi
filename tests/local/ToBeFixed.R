### Terrence D. Jorgensen
### Last updated: 10 February 2025
### remaining problems for lavaan.mi methods
### (both issues below are resolved, using the mi2lavaan() trick)


## -------------------------------
## fixed.x and conditional.x
## Using discretized example data
## -------------------------------

# data(binHS5imps)
#
# ## x9 (indicator) is incomplete, as well as x5 (predictor)
# catmod <- ' speed =~ x7 + x8 + x9
#             speed  ~ x4 + x5      '
# fitx <- cfa.mi(catmod, data = binHS5imps, fixed.x = TRUE, conditional.x = FALSE)
# fitxg <- cfa.mi(catmod, data = binHS5imps, fixed.x = TRUE, conditional.x = FALSE,
#                 group = "school")
# fit.x <- cfa.mi(catmod, data = binHS5imps, fixed.x = TRUE, conditional.x = TRUE)
# fit.xg <- cfa.mi(catmod, data = binHS5imps, fixed.x = TRUE, conditional.x = TRUE,
#                  group = "school")
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
# summary(fitx  , std="std.nox")
# summary(fitxg , std="std.nox")
# summary(fit.x , std="std.nox") # more efficient than cond.x=FALSE
# summary(fit.xg, std="std.nox")



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
# fitMeasures(fit.mi , output = "pretty")
# fitMeasures(fit2.mi, output = "pretty")
# fitMeasures(fit2.mi, output = "pretty", # no problem for other indices
#             fit.measures = c("chisq","df","pvalue","cfi","tli","rmsea"))




