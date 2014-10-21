\name{summary.saAlloc}
\alias{summary.saAlloc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizing an saAlloc Object 

}
\description{
 'summary' method for class 'saAlloc'
}
\usage{
summary.saAlloc(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  an object of class 'saAlloc', typically from a call to 'saMinCV' 
  or 'saEqualClus'.
  }
}
\details{
  There is no difference between the summaries of either method.
}
\value{
  This summary is simply a side effect, any values desired
  from the summary are contained in the original saAlloc object. 
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
# result

# print summary through s3 object
summary(b)

# print summary explicitly
summary.saAlloc(b)

}
