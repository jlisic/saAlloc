\name{saMinCV}
\alias{saMinCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal Strata Formation and Allocation Using Simulated Annealing
}
\description{
This function exchanges sampling units between strata to minimize the target 
coefficient of variation (CV), every few steps (as determined by 
\code{sampleUpdateIterations}) a number of iterations are performed to improve 
the target CV by changing the sample allocation.
}
\usage{
saMinCV(
     x, 
     label, 
     targetCV, 
     sampleSize, 
     weightMatrix, 
     domainMatrixList,
     iterations=1000, 
     sampleSizeIterations=3, 
     recalcIterations=10000,
     locationAdjustment,
     scaleAdjustment,
     p = 2,                   # l2 norm of the penalty function
     penalty = -1,            # negative penalties are ignored
     cooling = 0,
     preserveSatisfied=TRUE
   )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  (Required) A matrix where each row is a sampling unit and each column is a 
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
  (Required) Target CV, must be equal to the columns of x
  }
  \item{sampleSize}{
  (Required) sample size, either a scalar that will be equally split between
  strata, or a sample size for each strata.
  }
  \item{weightMatrix}{
  matrix of size \code{nrow(x)} by total strata, the higher the weight for a PSU for 
  a particular strata relative to the PSUs other strata weight the more likely
  it is to be moved to that strata. 
  }
  \item{domainMatrixList}{
  a list of matricies of length H, where each matrix is of size K by J. 
  }
  \item{iterations}{
  Number of iterations to perform of the Simulated Annealing procedure.
  }
  \item{sampleSizeIterations}{
  Number of iterations for each sample Allocation update step
  (0 = no changes).
  }
  \item{recalcIterations}{
  Number of iterations to wait until recalculating the variance 
  (0 = every iteration).
  }
  \item{locationAdjustment}{
  A constant value added to the population variable for each administrative variable, 
  this is of the same size and layout of x.  
  }
  \item{scaleAdjustment}{
  A constant value multiplied by the population variable (before the location adjustment) 
  for each administrative variable, this is of the same size and layout of x.  
  }
  \item{penalty}{
  The penalty vector used for each of the J penalties.  Penalties with negative
  elements are ignored.  A scalar value may be used, and will be repeated for all CVs. 
  }
  \item{p}{
  The exponent for the penalty function, fractional exponents are available.
  }
  \item{cooling}{
  Cooling Schedule (0 = exchange only)
  }
  \item{preserveSatisfied}{
  TRUE implies that any constraint that is met for an accepted state is kept,
  FALSE allows for met CVs to be violated by accepted states.
 } 
 }
\value{
  \item{accept}{ Matrix with three columns and a row for each iteration, 
    "change":   change in objective function, 
    "U":        uniform random variable, 
    "accepted": 1 - accepted, 
                0 - not accepted 
    "selected": index of PSU, 
    "from":     current stratum of selected PSU , 
    "to":       candidate stratum of selected PSU, 
    "n_h":      sample size (h is the index or label associated with strata),
    column names or numbers if not presented
  }
  \item{cost}{ Objective Function }
  \item{label}{ Final Label }
  \item{sampleSize}{ Final Sample Size }
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

summary(b)

}
