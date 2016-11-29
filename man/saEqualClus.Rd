\name{saEqualClus}
\alias{saEqualClus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal Equal Sized Strata Formation and Allocation Using Simulated Annealing 
}
\description{
This function minimizes the sample size for equally sized strata by exchanging
sampling units between strata.  Optimization is performed under variance 
constraints for known administrative variables.
}
\usage{
saEqualClus(
      x, 
      label, 
      iterations=1000, 
      cooling=0, 
      segments=rep(1,nrow(x)), 
      controlVariable=rep(1,nrow(x)), 
      targetVar, 
      varWithin=rep(0,ncol(x)), 
      tolSize=0.05
    )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  (Required) A matrix where each row is a sampling unit and each column is a 
  characteristic of the sampling unit that is to be reduced in optimal 
  sampling size.
  }
  \item{label}{
  (Required) Initial stratum assignment.
  }
  \item{iterations}{
  Number of iterations to perform of the simulated annealing procedure.
  }
  \item{cooling}{
  Cooling schedule (0 = exchange only).
  }
  \item{segments}{
  Number of secondary sampling units. 
  }
  \item{controlVariable}{
  Control variable that enforces equality constraints.  If a control variable
  is specified it must exist for each member of the population. 
  }
  \item{targetVar}{
  (Required) A vector of target variances for each administrative variable.  
  This vector must be of equal length to the number of columns of \code{x}. 
}
  \item{varWithin}{
  This parameter is used to account for within PSU variance.  It is a constant
  value added to the target variance and PSU variance at each iteration.
  }
  \item{tolSize}{
  Acceptable difference proportion between strata sums of the control variable.
  }
}
\value{
  \item{accept}{ Matrix with three columns and a row for each iteration, 

      "change": change in objective function, 

      "U": uniform random variable, 

      "accepted": 
      
      1 - accepted, 

      0 - not accepted 
  }
  \item{cost}{ Objective Function }
  \item{label}{ Final Label }
}

\author{
Jonathan Lisic <jonathan.lisic@nass.usda.gov> 
}
\examples{
# data set with 100 observations and four characteristics 
# split between two strata
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


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survey }
