\name{saAlloc-package}
\alias{saAlloc-package}
\alias{saAlloc}
\docType{package}
\title{
Optimal allocation of sampling units to form strata using simulated annealing 
}
\description{
Description: This program performs optimal allocation of sampling units to a set of strata via simulated annealing.  
Constraints on maximum relative distance between stratum sizes are supported. 
}
\details{
\tabular{ll}{
Package: \tab saAlloc\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-07-28\cr
License: \tab Public Domain\cr
}
 
saMinCV: currently exchange only, uses an exchange algorithm to form strata from a set of PSU's and simultaneously updated the allocation at each iteration. 

saEqualClus: allocate a set of observations to target strata based on minimizing total sample size based on a set of target variances (for fixed size strata only).

}
\author{
Jonathan Lisic

Maintainer: <jonathan.lisic@nass.usda.gov>
}
\keyword{ package }
\keyword{ survey }
\examples{
# data set with 100 observations and 4 characteristics split between two strata
x <- matrix(rnorm( 100*4 ),ncol=4)  

# create an initial set of strata assignments 
label <- c(rep(0,50),rep(1,50)) 

#target variance
targetVar <- c(10000, 20000,30000, 40000) 

# run equal clus
a <- saEqualClus(
  x,
  label,
  cooling=20,
  targetVar=targetVar
  )

# run minCV for a sample size of 20
b <- saMinCV(
  x,
  label,
  targetCV=c( 0.08, 0.10, 1.00, 0.05),
  sampleSize=20
)


}
