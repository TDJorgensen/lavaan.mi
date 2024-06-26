% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exampleData.R
\docType{data}
\name{binHS5imps}
\alias{binHS5imps}
\title{List of imputed Holzinger & Swineford (1939) dichotomized data}
\source{
The {lavaan} package.
}
\description{
A version of the classic Holzinger and Swineford (1939) dataset, with
missing data imposed on variables \code{x5} and \code{x9}:
}
\details{
\itemize{
\item \code{x5} is missing not at random (MNAR) by deleting the lowest 30\% of
\code{x5} values.
\item \code{x9} is missing at random (MAR) conditional on age, by deleting \code{x5}
values for the youngest 30\% of subjects in the data.
}

The data are then dichotomized using a median split, and imputed 5 times
using the syntax shown in the example. The data include only the 9 tests
(\code{x1} through \code{x9}) and school.
}
\examples{
 \dontrun{
data(HolzingerSwineford1939, package = "lavaan")

## impose missing data for example
HSMiss <- HolzingerSwineford1939[ , c(paste("x", 1:9, sep = ""),
                                      "ageyr","agemo","school")]
set.seed(123)
HSMiss$x5 <- ifelse(HSMiss$x5 <= quantile(HSMiss$x5, .3), NA, HSMiss$x5)
age <- HSMiss$ageyr + HSMiss$agemo/12
HSMiss$x9 <- ifelse(age <= quantile(age, .3), NA, HSMiss$x9)

## median split
HSbinary <- as.data.frame( lapply(HSMiss[ , paste0("x", 1:9)],
                                  FUN = cut, breaks = 2, labels = FALSE) )
HSbinary$school <- HSMiss$school

## impute binary missing data using mice package
library(mice)
set.seed(456)
miceImps <- mice(HSbinary)
## save imputations in a list of data.frames
binHS5imps <- list()
for (i in 1:miceImps$m) binHS5imps[[i]] <- complete(miceImps, action = i)
}

}
\references{
Holzinger, K., & Swineford, F. (1939).
\emph{A study in factor analysis: The stability of a bifactor solution}.
Supplementary Educational Monograph, no. 48.
Chicago, IL: University of Chicago Press.
}
\seealso{
\link[lavaan:HolzingerSwineford1939]{lavaan::HolzingerSwineford1939}
}
\author{
Terrence D. Jorgensen (University of Amsterdam;
\email{TJorgensen314@gmail.com})
}
\keyword{data}
