\name{saAlloc}
\alias{saAlloc}
\alias{plot.saAlloc}
\alias{summary.saAlloc}
\title{
saAlloc Object 
}
\description{
  The saAlloc object contains a set of results obtained by either by
  forming optimal strata and allocations for fixed CV constraints and 
  fixed total sample size (minCV) or without fixed sample size (saEqualClus). 
}
\usage{

  saAlloc(
     y, 
     label, 
     targetCV,
     targetVar, 
     sampleSize, 
     weightMatrix, 
     method='saMinCV',
     iterations=1000, 
     sampleIterations=100, 
     sampleUpdateIterations=100,
     cooling=0, 
     segments=rep(1,length(label)), 
     PSUAcres=rep(1,length(label)), 
     targetVarWithin=rep(0,ncol(y)),
     preserveSatisfied=TRUE,
     tolSize=1
   )

  \method{plot}{saAlloc}(x, together = TRUE,...)
  \method{summary}{saAlloc}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  A matrix where each row is a sampling unit and each column is a 
  characteristic of the sampling unit that is to be reduced in optimal sampling
  size.  

  Multiple observations of the same PSU can be included by including additional
  rows for each observation.
 
  e.g. for 3 observations per PSU:

  PSU  Observation of PSU   .... data
  1    1
  1    2
  1    3
  2    1
  2    2
  2    3
  .    . 
  .    .
  .    .

  Note: that only x should be of this form for multiple observations, labels 
  etc...  should be of the same size as if there was a single observation 
  per PSU.
  }
  \item{label}{
  (Required) Initial stratum assignment.
  }
  \item{targetCV}{
  (Required) Target CV, must be equal to the columns of y 
  }
  \item{targetVar}{
  (Required) Target Variance, must be equal to the columns of y 
  }
  \item{sampleSize}{
  (Required) sample size, either a scalar that will be equally split between
  strata, or a sample size for each strata.
  }
  \item{weightMatrix}{
  matrix of size nrow(x) by total strata, the higher the weight for a PSU for 
  a particular strata relative to the PSU's other strata weight the more likely
  it is to be moved to that strata. 
  }
  \item{method}{
    One of two methods 'saEqualClus' or 'saMininCV'.
  }
  \item{iterations}{
  Number of iterations to perform of the Simulated Annealing procedure.
  }
  \item{sampleIterations}{
  Number of iterations for each sample Allocation update step
  (0 = no changes).
  }
  \item{sampleUpdateIterations}{
  Number of iterations to wait until performing the sample update step 
  (0 = every iteration).
  }
  \item{cooling}{
  Cooling Schedule (0 = exchange only)
  }
  \item{segments}{
  If prior clustering is performed to form sampling units, this allows 
  }
  \item{PSUAcres}{
  PSU size (originally done in acres), 
  }
  \item{targetVarWithin}{
  If prior clustering is performed and assumed constant across the population 
  then it can be handled here.
  }
  \item{preserveSatisfied}{
   TRUE implies that any constraint that is met for an accepted state is kept,
   FALSE allows for met CVs to be violated by accepted states.
  } 
  \item{x}{
  an object of class 'saAlloc', typically from a call to 'saAlloc', 'saMinCV' 
  or 'saEqualClus'.
  }
  \item{object}{
  an object of class 'saAlloc', typically from a call to 'saAlloc', 'saMinCV' 
  or 'saEqualClus'.
  }
  \item{together}{
  if TRUE then include all traces on a single plot, otherwise
  use a separate plot for each trace. 
  }
  \item{...}{ other parameters to be passed through the function
   function. 
  }
  \item{tolSize}{
    Acceptable difference between PSUacres.
  }
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

