### Terrence D. Jorgensen
### Last updated: 24 April 2024
### generate imputed data set for documentation examples

# data(HolzingerSwineford1939, package = "lavaan")
# ## impose missing data for example
# HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
#                                       "ageyr","agemo","school")]
# set.seed(123)
# HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
# age <- HSMiss$ageyr + HSMiss$agemo/12
# HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
# ## impute missing data with Amelia
# library(Amelia)
# set.seed(456)
# HS.amelia <- amelia(HSMiss, m = 20, noms = "school", p2s = FALSE)
# HS20imps <- HS.amelia$imputations
# save(HS20imps, file = "data/HS20imps.rda")
#
# ## ordered-categorical data
# HSbinary <- as.data.frame( lapply(HSMiss[ , paste0("x", 1:9)],
#                                   FUN = cut, breaks = 2, labels = FALSE) )
# HSbinary$school <- HSMiss$school
#
# ## impute binary missing data using mice package
# library(mice)
# set.seed(456)
# miceImps <- mice(HSbinary)
# ## save imputations in a list of data.frames
# binHS5imps <- list()
# for (i in 1:miceImps$m) binHS5imps[[i]] <- complete(miceImps, action = i)
# save(binHS5imps, file = "data/binHS5imps.rda")


##' List of imputed Holzinger & Swineford (1939) datasets
##'
##' A version of the classic Holzinger and Swineford (1939) dataset, with
##' missing data imposed on variables `x5` and `x9`:
##'
##'   - `x5` is missing not at random (MNAR) by deleting the lowest 30% of
##'     `x5` values.
##'   - `x9` is missing at random (MAR) conditional on age, by deleting `x5`
##'     values for the youngest 30% of subjects in the data.
##'
##' The data are imputed 20 times using the syntax shown in the example.
##' The data include only age and school variables, along with 9 tests
##' (`x1` through `x9`).
##'
##' @source The {lavaan} package.
##'
##' @examples
##'  \dontrun{
##' data(HolzingerSwineford1939, package = "lavaan")
##'
##' ## impose missing data for example
##' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
##'                                       "ageyr","agemo","school")]
##' set.seed(123)
##' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
##' age <- HSMiss$ageyr + HSMiss$agemo/12
##' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
##'
##' ## impute missing data with Amelia
##' library(Amelia)
##' set.seed(456)
##' HS.amelia <- amelia(HSMiss, m = 20, noms = "school", p2s = FALSE)
##' HS20imps <- HS.amelia$imputations
##' }
##'
##' @name HS20imps
##' @docType data
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##' @references
##'   Holzinger, K., & Swineford, F. (1939).
##'   *A study in factor analysis: The stability of a bifactor solution*.
##'   Supplementary Educational Monograph, no. 48.
##'   Chicago, IL: University of Chicago Press.
##' @seealso [lavaan::HolzingerSwineford1939]
##' @keywords data
NULL




##' List of imputed Holzinger & Swineford (1939) dichotomized data
##'
##' A version of the classic Holzinger and Swineford (1939) dataset, with
##' missing data imposed on variables `x5` and `x9`:
##'
##'   - `x5` is missing not at random (MNAR) by deleting the lowest 30% of
##'     `x5` values.
##'   - `x9` is missing at random (MAR) conditional on age, by deleting `x5`
##'     values for the youngest 30% of subjects in the data.
##'
##' The data are then dichotomized using a median split, and imputed 5 times
##' using the syntax shown in the example. The data include only the 9 tests
##' (`x1` through `x9`) and school.
##'
##' @source The {lavaan} package.
##'
##' @examples
##'  \dontrun{
##' data(HolzingerSwineford1939, package = "lavaan")
##'
##' ## impose missing data for example
##' HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
##'                                       "ageyr","agemo","school")]
##' set.seed(123)
##' HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
##' age <- HSMiss$ageyr + HSMiss$agemo/12
##' HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)
##'
##' ## median split
##' HSbinary <- as.data.frame( lapply(HSMiss[ , paste0("x", 1:9)],
##'                                   FUN = cut, breaks = 2, labels = FALSE) )
##' HSbinary$school <- HSMiss$school
##'
##' ## impute binary missing data using mice package
##' library(mice)
##' set.seed(456)
##' miceImps <- mice(HSbinary)
##' ## save imputations in a list of data.frames
##' binHS5imps <- list()
##' for (i in 1:miceImps$m) binHS5imps[[i]] <- complete(miceImps, action = i)
##' }
##'
##' @name binHS5imps
##' @docType data
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##' @references
##'   Holzinger, K., & Swineford, F. (1939).
##'   *A study in factor analysis: The stability of a bifactor solution*.
##'   Supplementary Educational Monograph, no. 48.
##'   Chicago, IL: University of Chicago Press.
##' @seealso [lavaan::HolzingerSwineford1939]
##' @keywords data
NULL




