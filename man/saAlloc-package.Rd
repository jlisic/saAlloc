\name{saAlloc-package}
\alias{saAlloc-package}
\docType{package}
\title{
Optimal Allocation of Sampling Units to Form Strata Using Simulated Annealing 
}
\description{
A suite of tools to perform joint optimal allocation and stratification of sampling
units using simulated annealing.  Soft and hard constraints placed on the coefficient 
of variances (CV) or variances of administrative variables are supported.
}
\details{
\tabular{ll}{
Package: \tab saAlloc\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2016-10-31\cr
License: \tab Unlimited\cr
}
 
saMinCV: uses an exchange algorithm with simulated annealing to form strata from a 
  set of sampling units and simultaneously updates the allocation at each iteration. 

saEqualClus: allocates a set of sampling units to strata based on minimizing
  total sample size based on a set of target variances (for appromately equally sized
  strata only).
  
saSampleAlloc: uses simulated annealing to perform optimal allocation.

}
\author{
Jonathan Lisic

Maintainer: <jonathan.lisic@nass.usda.gov>
}
\keyword{ survey }
\examples{
# data set with 100 observations and 4 characteristics split between two strata
x <- matrix(rnorm( 100*4 ),ncol=4)  

# create an initial set of strata assignments 
label <- c(rep(0,50),rep(1,50)) 

#target variance
targetVar <- c(10000, 20000,30000, 40000) 

# run equal cluster
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
\references{
Lisic, J.J. , J. Sang, Z. Zhu, and S. Zimmer. (2015) \emph{Optimal stratification and allocation for the June agricultural survey}.  Presented at the Federal Committee on Statistical Methodology Research Conference, Washington DC.  
}
