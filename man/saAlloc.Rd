\name{saAlloc}
\alias{saAlloc}
\alias{plot.saAlloc}
\alias{summary.saAlloc}
\title{
saAlloc Object 
}
\description{
  The saAlloc object contains a set of results obtained by either forming optimal 
  strata and allocations for fixed CV constraints and fixed total sample size \code{minCV}
  or without fixed sample size \code{saEqualClus}. 
}
\author{
  Jonathan Lisic
}
\seealso{
  \code{\link{saMinCV}}
  \code{\link{saEqualClus}}
}
\examples{

# data set with 100 observations and four 
# characteristics split between two strata
x <- matrix(rnorm( 100*4 ),ncol=4)  

# create an initial set of strata assignments 
label <- c(rep(0,50),rep(1,50)) 

# run minCV for a sample size of 20
b <- saMinCV(
  x,
  label,
  targetCV=c(0.08, 0.10, 0.05, 0.08),
  sampleSize=20
)

# both of these functions return the same 
# plot

# plot through s3 object
plot(b)

# plot explicitly
summary(b)

}

