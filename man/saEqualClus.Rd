\name{saEqualClus}
\alias{saEqualClus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal Equal Sized Strata Formation and Allocation via Simulated Annealing 
}
\description{
This function minimizes the sample size for a fixed set of target variances by moving sampling units between strata.
}
\usage{
saEqualClus(
      x, 
      label, 
      iterations=1000, 
      cooling=0, 
      segments=rep(1,nrow(x)), 
      PSUAcres=rep(1,nrow(x)), 
      targetVar, 
      targetVarWithin=rep(0,ncol(x)), 
      tolSize=0.05
    )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  A matrix where each row is a sampling unit and each column is a characteristic of the sampling unit that is to be reduced in optimal sampling size.
  }
  \item{label}{
  Initial stratum assignment.
  }
  \item{iterations}{
  Number of iterations to perform of the Simulated Annealing procedure.
  }
  \item{cooling}{
  Cooling Schedule (0 = exchange only).
  }
  \item{segments}{
  If prior clustering is performed to form sampling units, this allows 
  }
  \item{PSUAcres}{
  PSU size (originally done in acres), 
  }
  \item{targetVar}{
%%     ~~Describe \code{targetVar} here~~
}
  \item{targetVarWithin}{
  If prior clustering is performed and assumed constant across the population then it can be handled here.
  }
  \item{tolSize}{
    Acceptable difference between PSUacres.
  }
}
\value{
  \item{accept}{ Matrix with three columns and a row for each iteration, "change": change in objective function, "U": uniform random variable, "accepted": 1 - accepted, 0 - not accepted }
  \item{cost}{ Objective Function }
  \item{label}{ Final Label }
}

\author{
Jonathan Lisic <jonathan.lisic@nass.usda.gov> 
}
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


}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
