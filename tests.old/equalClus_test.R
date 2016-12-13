# example for saEqualClus
library(saAlloc)

# data set with 100 observations and four characteristics
# split between two strata
x <- matrix(rnorm( 100*4 ),ncol=4)

# create an initial set of strata assignments
label <- c(rep(0,50),rep(1,50))

#target variance
targetVar <- c(0.75, 1, 2, 3)

# run equal cluster
a <- saEqualClus(
  x,
  label,
  iter= 100,
  cooling=20,
  targetVar=targetVar
)

print(
sqrt(colSums(aggregate(x,by=list(label),var)[,-1]))
/
colSums(x)
)

print(
sqrt(colSums(aggregate(x,by=list(a$label),var)[,-1]))/
colSums(x)
)

print( head(a$accept))
