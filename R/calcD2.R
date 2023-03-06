### Terrence D. Jorgensen
### Last updated: 5 May 2022
### function to calculate D2 pooled statistic

##' Calculate the "D2" statistic
##'
##' This is a utility function used to calculate the "D2" statistic for pooling
##' test statistics across multiple imputations. This function is called by
##' several functions used for \code{\linkS4class{lavaan.mi}} objects, such as
##' \code{\link{lavTestLRT.mi}}, \code{\link{lavTestWald.mi}}, and
##' \code{\link{lavTestScore.mi}}. But this function can be used for any general
##' scenario because it only requires a vector of \eqn{\chi^2} statistics (one
##' from each imputation) and the degrees of freedom for the test statistic.
##' See Li, Meng, Raghunathan, & Rubin (1991) and Enders (2010, chapter 8) for
##' details about how it is calculated.
##'
##' @importFrom stats var pf pchisq
##'
##' @param w \code{numeric} vector of Wald \eqn{\chi^2} statistics. Can also
##'   be Wald \emph{z} statistics, which will be internally squared to make
##'   \eqn{\chi^2} statistics with one \emph{df} (must set \code{DF = 0L}).
##' @param DF degrees of freedom (\emph{df}) of the \eqn{\chi^2} statistics.
##'   If \code{DF = 0L} (default), \code{w} is assumed to contain \emph{z}
##'   statistics, which will be internally squared.
##' @param asymptotic \code{logical}. If \code{FALSE} (default), the pooled test
##'   will be returned as an \emph{F}-distributed statistic with numerator
##'   (\code{df1}) and denominator (\code{df2}) degrees of freedom.
##'   If \code{TRUE}, the pooled \emph{F} statistic will be multiplied by its
##'   \code{df1} on the assumption that its \code{df2} is sufficiently large
##'   enough that the statistic will be asymptotically \eqn{\chi^2} distributed
##'   with \code{df1}.
##'
##' @return A \code{numeric} vector containing the test statistic, \emph{df},
##'   its \emph{p} value, and 2 missing-data diagnostics: the relative invrease
##'   in variance (RIV, or average for multiparameter tests: ARIV) and the
##'   fraction missing information (FMI = ARIV / (1 + ARIV)).
##'
##' @seealso \code{\link{lavTestLRT.mi}}, \code{\link{lavTestWald.mi}},
##'   \code{\link{lavTestScore.mi}}
##'
##' @author Terrence D. Jorgensen (University of Amsterdam;
##'   \email{TJorgensen314@@gmail.com})
##'
##' @references
##'   Enders, C. K. (2010). \emph{Applied missing data analysis}. New
##'   York, NY: Guilford.
##'
##'   Li, K.-H., Meng, X.-L., Raghunathan, T. E., & Rubin, D. B. (1991).
##'   Significance levels from repeated \emph{p}-values with multiply-imputed
##'   data. \emph{Statistica Sinica, 1}(1), 65--92. Retrieved from
##'   \url{https://www.jstor.org/stable/24303994}
##'
##' @examples
##' ## generate a vector of chi-squared values, just for example
##' DF <- 3 # degrees of freedom
##' M <- 20 # number of imputations
##' CHI <- rchisq(M, DF)
##'
##' ## pool the "results"
##' calculate.D2(CHI, DF) # by default, an F statistic is returned
##' calculate.D2(CHI, DF, asymptotic = TRUE) # asymptotically chi-squared
##'
##' ## generate standard-normal values, for an example of Wald z tests
##' Z <- rnorm(M)
##' calculate.D2(Z) # default DF = 0 will square Z to make chisq(DF = 1)
##' ## F test is equivalent to a t test with the denominator DF
##'
##'
##' @export
calculate.D2 <- function(w, DF = 0L, asymptotic = FALSE) {
  if (length(w) == 0L) return(NA)
  w <- as.numeric(w)
  DF <- as.numeric(DF)

  nImps <- sum(!is.na(w))
  if (nImps == 0) return(NA)

  if (DF <= 0L) {
    ## assume they are Z scores
    w <- w^2
    DF <- 1L
  }

  ## pool test statistics
  if (length(w) > 1L) {
    w_bar <- mean(w, na.rm = TRUE)
    ariv <- (1 + 1/nImps) * var(sqrt(w), na.rm = TRUE)
    test.stat <- (w_bar/DF - ((nImps + 1) * ariv / (nImps - 1))) / (1 + ariv)
  } else {
    warning('There was only 1 non-missing value to pool, leading to zero ',
            'variance, so D2 cannot be calculated.')
    test.stat <- ariv <- NA
  }
  if (test.stat < 0) test.stat <- 0
  if (asymptotic) {
    out <- c("chisq" = test.stat * DF, df = DF,
             pvalue = pchisq(test.stat * DF, df = DF, lower.tail = FALSE),
             ariv = ariv, fmi = ariv / (1 + ariv))
  } else {
    v3 <- DF^(-3 / nImps) * (nImps - 1) * (1 + (1 / ariv))^2
    out <- c("F" = test.stat, df1 = DF, df2 = v3,
             pvalue = pf(test.stat, df1 = DF, df2 = v3, lower.tail = FALSE),
             ariv = ariv, fmi = ariv / (1 + ariv))
  }
  class(out) <- c("lavaan.vector","numeric")
  out
}

