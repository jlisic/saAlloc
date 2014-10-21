\name{plot.saAlloc}
\alias{plot.saAlloc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Diagnostics for an saAlloc Object 
}
\description{
  Produced a single plot for the trace of the objective function, or
  multiple plots on a single plotting device for each marginal objective
  function.  Multiple plots are only available for saMinCV generated
  objects.
}
\usage{
  plot.saAlloc(x, together = TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  an object of class 'saAlloc', typically from a call to 'saMinCV' 
  or 'saEqualClus'.
}
  \item{together}{
  if TRUE then include all traces on a single plot, otherwise
  use a separate plot for each trace. 
  }
  \item{...}{ other parameters to be passed through to the plotting
   function. 
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
plot.saAlloc(b)

}
