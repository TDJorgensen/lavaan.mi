### Terrence D. Jorgensen
### Last updated: 22 November 2022
### pool imputed polychorics before fitting model in "single" step:
###    Normal data: https://doi.org/10.3102/1076998612458320
###    Categorical: https://doi.org/10.1080/00273171.2018.1523000
### using example from ?lavaan.mi help page


library(lavaan.mi)

## impose missing data
data(datCat, package = "semTools")
set.seed(123)
for (i in 1:8) datCat[sample(1:nrow(datCat), size = .1*nrow(datCat)), i] <- NA

## impute ordinal missing data using mice package
library(mice)
set.seed(456)
miceImps <- mice(datCat)
## save imputations in a list of data.frames
impList <- list()
for (i in 1:miceImps$m) impList[[i]] <- complete(miceImps, action = i)

## highest category of u4 empty in group 2,
## treat as continuous in 2-group SEM
impList1 <- lapply(impList, function(x) {
  x$g <- NULL
  x
})
impList2 <- lapply(impList, function(x) {
  x$u4 <- as.integer(x$u4)
  x
})

## fit model directly to multiple imputations
mod <- ' f1 =~ u1 + u2 + u3 + u4
         f2 =~ u5 + u6 + u7 + u8 '
fit <- cfa.mi(mod, data = impList)
summary(fit, fit = TRUE, test = "D2")

fitmg <- cfa.mi(mod, data = impList2, group = "g")
summary(fitmg, fit = TRUE, test = "D2")




## --------------------
## Pool saturated model
## --------------------

## use lavCor on each imputation
satList <- lapply(impList1, function(x) {
  lavCor(x, output = "fit", estimator = "DWLS", se = "robust.sem", h1 = TRUE)
})
m <- length(satList)

## or use lavaan.mi to pool
sat.1g <- lavaan.mi(parTable(satList[[1]]), data = impList)


## average coefficients
# sample.cov  <- Reduce("+", lapply(satList, lavInspect, what = "cov.ov")) / m
# sample.mean <- rowMeans(   sapply(satList, lavInspect, what = "mean.ov"))
# sample.th   <- rowMeans(   sapply(satList, lavInspect, what = "thresholds"))
# # from any imputation
# attr(sample.th, "th.idx") <- lavInspect(satList[[1]], "th.idx")
# N <- lavInspect(satList[[1]], "nobs")
moments <- fitted(sat.1g)
sample.cov  <- moments$cov
sample.mean <- moments$mean
sample.th   <- moments$th
attr(sample.th, "th.idx") <- lavInspect(sat.1g, "th.idx")
N <- lavListInspect(sat.1g, "nobs")

## average ACOVs
# varW <- Reduce("+", lapply(satList, vcov)) / m
# ## (co)variances of coefficients
# varB <- var(t(sapply(satList, coef)))
# ## pooled (co)variances
# varT <- varW + (1 + 1/m)*varB
varT <- vcov(sat.1g)

## format uncertainty/weights
free.idx <- lavInspect(satList[[1]], "free")
nn0      <- lavNames(satList[[1]])
mth.idx  <- lavInspect(satList[[1]], "th.idx")
num.idx  <- setdiff(seq_along(nn0), mth.idx)
mth.idx[mth.idx != 0] <- free.idx$tau[ , 1]
mth.idx[mth.idx == 0] <- free.idx$nu[num.idx, 1]
wls.idx <- c(mth.idx, # interleaved thresholds + (negative) means
             diag(free.idx$theta)[num.idx], # variances
             lav_matrix_vech(free.idx$theta, diagonal = FALSE)) # covariances
## apply a sign-switch for rows/columns of means with other rows/columns
neg1 <- matrix(1, nrow = nrow(varT), ncol = ncol(varT))
if (length(num.idx)) {
  neg1[  free.idx$nu[num.idx, 1] , -(free.idx$nu[num.idx, 1])] <- -1
  neg1[-(free.idx$nu[num.idx, 1]),   free.idx$nu[num.idx, 1] ] <- -1
}

NACOV <- (N-1L)*(neg1*varT)[wls.idx, wls.idx]
WLS.V <- diag(diag(solve(NACOV)))
dimnames(WLS.V) <- dimnames(NACOV)
class(NACOV) <- c("lavaan.matrix.symmetric","matrix")
class(WLS.V) <- c("lavaan.matrix.symmetric","matrix")

## compare to function's output (load below):
datPre1 <- poolSat(impList1)
all.equal(NACOV, datPre1$NACOV)
all.equal(WLS.V, datPre1$WLS.V)


## fit model to moments
fitpre <- cfa(mod,
              sample.cov  = datPre1$sample.cov,
              sample.mean = datPre1$sample.mean,
              sample.nobs = datPre1$sample.nobs, #N,
              sample.th   = datPre1$sample.th,
              WLS.V       = datPre1$WLS.V,
              NACOV       = datPre1$NACOV)
summary(fitpre, fit.measures = TRUE)

## using data= argument
fitpre <- cfa(mod, data = datPre1, estimator = "WLSMV")



## -----------------------------
## TODO: Repeat with group = "g"
## -----------------------------

## fit using lavaan.mi
satTemp <- lavCor(impList2[[1]], do.fit = FALSE, group = "g",
                  output = "fit", estimator = "DWLS",
                  se = "robust.sem", h1 = TRUE)
PT <- parTable(satTemp)
sat.mg <- lavaan.mi(PT, data = impList2,
                    #TODO: could use available row/colnames?
                    FUN = function(obj) list(WLS.V = lavInspect(obj, "WLS.V"),
                                             NACOV = lavInspect(obj, "gamma")),
                    group = "g")


## average coefficients
mg.moments <- fitted(sat.mg)
sample.cov  <- sapply(mg.moments, "[[", i = "cov", simplify = FALSE)
sample.mean <- sapply(mg.moments, "[[", i = "mean", simplify = FALSE)
sample.th   <- sapply(mg.moments, "[[", i = "th", simplify = FALSE)
attr(sample.th, "th.idx") <- lavInspect(sat.mg, "th.idx")

## average ACOVs
varT <- vcov(sat.mg) # gets partitioned using each group's "free" indices


N <- lavInspect(sat.mg, "nobs")
ngroups <- lavInspect(sat.mg, "ngroups")

NACOV <- WLS.V <- vector("list", length = ngroups)
for (g in 1:ngroups) {
  free.idx <- lavInspect(sat.mg, "free")[[g]]
  nn0      <- lavNames(sat.mg, group = g)
  mth.idx  <- lavInspect(sat.mg, "th.idx")[[g]]
  num.idx  <- setdiff(seq_along(nn0), mth.idx)
  mth.idx[mth.idx != 0] <- free.idx$tau[ , 1]
  mth.idx[mth.idx == 0] <- free.idx$nu[num.idx, 1]
  wls.idx <- c(mth.idx, # interleaved thresholds + (negative) means
               diag(free.idx$theta)[num.idx], # variances
               lav_matrix_vech(free.idx$theta, diagonal = FALSE)) # covariances
  ## apply a sign-switch for rows/columns of means with other rows/columns
  neg1 <- matrix(1, nrow = nrow(varT), ncol = nrow(varT))
  if (length(num.idx)) {
    neg1[  free.idx$nu[num.idx, 1] , -(free.idx$nu[num.idx, 1])] <- -1
    neg1[-(free.idx$nu[num.idx, 1]),   free.idx$nu[num.idx, 1] ] <- -1
  }
  NACOV[[g]] <- (N[g] - 1L) * (neg1*varT)[wls.idx, wls.idx]
  ## use group-1's row/colnames (no .g2 suffix)
  if (g > 1L) dimnames(NACOV[[g]]) <- dimnames(NACOV[[1]])
  WLS.V[[g]] <- diag(diag(solve(NACOV[[g]])))
  dimnames(WLS.V[[g]]) <- dimnames(NACOV[[g]])

  class(NACOV[[g]]) <- c("lavaan.matrix.symmetric","matrix")
  class(WLS.V[[g]]) <- c("lavaan.matrix.symmetric","matrix")
}


## compare to function's output (load below):
datPre2 <- poolSat(impList2, group = "g")
all.equal(NACOV, datPre2$NACOV)
all.equal(WLS.V, datPre2$WLS.V)


## fit model to moments
mgpre <- cfa(mod,
             sample.cov  = datPre2$sample.cov,
             sample.mean = datPre2$sample.mean,
             sample.nobs = datPre2$sample.nobs, #N,
             sample.th   = datPre2$sample.th,
             WLS.V       = datPre2$WLS.V,
             NACOV       = datPre2$NACOV)
#FIXME: ask Yves why this issues a warning about starting values
summary(mgpre, fit.measures = TRUE)
lavTest(mgpre, test = "Browne.residual.adf")





