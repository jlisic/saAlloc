\name{saSampleAlloc}
\alias{saSampleAlloc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal Allocation Using Simulated Annealing
}
\description{
This function exchanges sample allocations between strata to minimize the CV of
a set of known administrative variables.  This is accomplished by through simulated
annealing, where proposed changes in allocations are accepted or rejected based on 
a penalized objective function.  The penalized objective function takes the form of 
the sum of an L2 norm and a penalty function.  The first L2 norm is the vector of
differences between the CVs for the current allocation and the target CVs.  The 
penalty function is the dot product between a vector of penalty values and the 
vector of CVs for the current allocation.
}
\usage{
saSampleAlloc(
  total,
  S2,
  targetCV,
  sampleSize,
  strataSize, 
  locationAdjustment,
  scaleAdjustment,
  iterations=100,
  p = 2,                   # l2 norm of the penalty function
  penalty = -1,            # negative penalties are ignored
  cooling = 0,
  preserveSatisfied=TRUE,
  fpc = TRUE
   )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{total}{
  (Required) A matrix where each row is a replicate population total and 
  each column is a characteristic of the sampling population used for
  allocation.    

  When replicates are used, the average CV of each characteristic is used 
  in the objective function.
  }
  \item{S2}{
  (Required) A matrix where each row identifies a stratum (in order) and 
  each column is a characteristic of the sampling population used for 
  allocation.  The values of this matrix is the variance estimate for each 
  sampling unit within each stratum (S^2).

  When replicates are used for total, a variance must also be provided for
  each replicate.  This is done by appending each replicate matrix of 
  variances to the S2 matrix in the same order it appears in total.
  }
  \item{targetCV}{
  (Required) A vector of target CVs.  This vector must be of equal length 
  to the number of columns of \code{total}. 
  }
  \item{sampleSize}{
  (Required) sample size, either an integer valued scalar that will be 
  equally split between strata, or a sample size for each strata.  When 
  the sample size cannot be split excess sample allocations will be given 
  to strata in increasing order.
  }
  \item{strataSize}{
  (Required) strata size, an integer valued array or population sizes by
  stratum. 
  }
  \item{locationAdjustment}{
  A constant value added to the population variable for each administrative 
  variable, this is of the same size and layout of \code{S2}.  
  }
  \item{scaleAdjustment}{
  A constant value multiplied by the population variable (before the location 
  adjustment) for each administrative variable, this is of the same size and 
  layout of \code{S2}.  
  }
  \item{iterations}{
  Number of iterations of the allocation. 
  }
  \item{p}{
  The exponent for the penalty function, fractional exponents are available.
  }
  \item{penalty}{
  The penalty vector used for each of the CV targets penalties.  Penalties 
  with negative elements are ignored.  A scalar value may be used, and will 
  be repeated for all CVs. 
  }
  \item{cooling}{
  Cooling Schedule (0 = exchange only).
  }
  \item{preserveSatisfied}{
  \code{TRUE} implies that any constraint that is met for an accepted state is kept,
  \code{FALSE} allows for met CVs to be violated by accepted states.
 } 
  \item{fpc}{
  \code{TRUE} use a finite population correction factor (fpc),
  \code{FALSE} do not use a fpc. 
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

total <- colSums(x)

S2 <- apply( x, 2, var)

# run minCV for a sample size of 20
b <- saSampleAlloc(
  total,
  S2,
  strataSize=rep(50,2),
  targetCV=c(0.08, 0.10, 0.05, 0.08),
  sampleSize=20
)

summary(b)

}
