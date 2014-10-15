\name{saMinCV}
\alias{saMinCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimum CV Strata Formation and Allocation via Simulated Annealing
}
\description{
This function exchanges sampling units between strata to minimzie the target CV, for each accepted exchange a number of iterations are performed to improve the target CV by changing the sample allocation.
}
\usage{
saMinCV (x, label, targetCV, sampleSize, weightMatrix, iterations=1000, sampleIterations=100, cooling=0, segments=rep(1,nrow(x)), PSUAcres=rep(1,nrow(x)), targetVarWithin=rep(0,ncol(x)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  (Required) A matrix where each row is a sampling unit and each column is a characteristic of the sampling unit that is to be reduced in optimal sampling size.
  }
  \item{label}{
  (Required) Initial stratum assignment.
  }
  \item{targetCV}{
  (Required) Target CV, must be equal to the columns of x
  }
  \item{sampleSize}{
  (Required) sample size, either a scalar that will be equally split between strata, or a sample size for each strata.
  }
  \item{weightMatrix}{
  matrix of size nrow(x) by total strata, the higher the weight for a PSU for a particular strata relative to the PSU's other strata weight the more likely it is to be moved to that strata. 
  }
  \item{iterations}{
  Number of iterations to perform of the Simulated Annealing procedure.
  }
  \item{sampleIterations}{
  Number of iterations to perform on minimizing the sample Allocation (0 = no changes).
  }
  \item{cooling}{
  Cooling Schedule (0 = exchange only, as of this release only exchange is supported).
  }
  \item{segments}{
  If prior clustering is performed to form sampling units, this allows 
  }
  \item{PSUAcres}{
  PSU size (originally done in acres), 
  }
  \item{targetVarWithin}{
  If prior clustering is performed and assumed constant across the population then it can be handled here.
  }
}
\value{
  \item{accept}{ Matrix with three columns and a row for each iteration, "change": change in objective function, "U": uniform random variable, "accepted": 1 - accepted, 0 - not accepted }
  \item{cost}{ Objective Function }
  \item{label}{ Final Label }
  \item{sampleSize}{ Final Sample Size }
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

# data set with 100 observations and 4 characteristics split between two strata
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


# compare results
print("before CV")
print( sqrt(colSums(50*aggregate(x,by=list(label),var))[-1]) / colSums(x) )

print("Min CV")
print( sqrt(colSums(50*aggregate(x,by=list(b$label),var))[-1]) / colSums(x) )

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
