\name{saMinCV}
\alias{saMinCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Minimum CV Strata Formation via Simulated Annealing
}
\description{
This function exchanges sampling units between strata to minimzie the target CV.
}
\usage{
saMinCV (x, label, iterations=1000, cooling=0, segments=rep(1,nrow(x)), PSUAcres=rep(1,nrow(x)), targetCV, sampleSize, targetVarWithin=rep(0,ncol(x)), tolSize=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  (Required) A matrix where each row is a sampling unit and each column is a characteristic of the sampling unit that is to be reduced in optimal sampling size.
  }
  \item{label}{
  (Required) Initial stratum assignment.
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
  \item{targetCV}{
  (Required) Target CV, must be equal to the columns of x
  }
  \item{sampleSize}{
  (Required) sample size, either a scalar that will be equally split between strata, or a sample size for each strata.
  }
  \item{targetVarWithin}{
  If prior clustering is performed and assumed constant across the population then it can be handled here.
  }
  \item{tolSize}{
  Acceptable difference between PSUacres. 
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{accept}{ Matrix with three columns and a row for each iteration, "change": change in objective function, "U": uniform random variable, "accepted": 1 - accepted, 0 - not accepted }
  \item{cost}{ Objective Function }
  \item{label}{ Final Label }
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Jonathan Lisic <jonathan.lisic@nass.usda.gov>
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
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
