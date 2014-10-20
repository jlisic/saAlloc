library(saAlloc)

################# functions ######################

# function to calculate cv
cv <- function( x, strata = rep(1,nrow(x)), sampleSize=1) {
   
  # get standard deviation
  vars <- aggregate( x,by=list(strata), var)[,-1]
 
  # get number of strata
  H <- length(unique(strata))                                                                                                                                               
    # if only one sample is provided, samples are distributed uniformly
  if( length(sampleSize) != H) {
       
    if( length(sampleSize) != 1 ) stop("Invalid number of sample sizes")
  
    sampleSize = rep(sampleSize/H,H)
  }
   
  # get total
  totals <- colSums( x )
  Nh <- aggregate( x,by=list(strata), length)[,-1]
  return(sqrt(colSums(vars * Nh^2/sampleSize)) / totals )
}


################# generate data ######################

set.seed(400)

nHigh   <- 40
nMed    <- 20 
nLow <- 60 

# high density cultivation, few farms
highCult <- cbind( 
                rnorm(nHigh, 500, 80),
                rlnorm(nHigh, 1,.3)
              )

# medium density cultivation
mediumCult <- cbind( 
                rnorm(nMed,200,100),
                rlnorm(nMed,1.5,.3)
              )

# low density cultivation, many farms and non-farms 
lowCult <- cbind( 
             rnorm( nLow,10,40),
             rlnorm( nLow,2)
           )
               

# handle possible negative values
highCult[highCult[,1] < 0,1] <- 0 
mediumCult[mediumCult[,1] < 0,1] <- 0 
lowCult[lowCult[,1] < 0,1] <- 0 

highCult[highCult[,2] > 5,2] <- 10 
mediumCult[mediumCult[,2] > 10,2] <- 10 
lowCult[lowCult[,2] > 40 ,2] <- 40 

# create population
x <- rbind( highCult, mediumCult, lowCult)


################# set parameters ######################

N <- nrow(x)
H <- 3 
sampleSize <- rep( 4, H) 
iterations <- 50000
sampleIterations <- 100
sampleUpdateIterations <- 100
cool = 3
cvTargetX = 0.1
cvTargetY = 0.1
 
# initial k-means clutering
kMeansCluster <- kmeans(x, H)$cluster 

# sa strata
saStrata <- saMinCV(
  x,
  kMeansCluster,
  iterations=iterations,
  sampleIterations=sampleIterations,
  sampleUpdateIterations=sampleUpdateIterations,
  cooling=cool,
  targetCV=c(cvTargetX, cvTargetY),
  sampleSize=sampleSize
)

# calculate cv
cv.minCV <- cv(
  x,
  strata=currentCluster,
  sampleSize=clustersOpt()$sampleSize 
)






