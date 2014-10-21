library(saAlloc)



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
iterations <- 500000
sampleIterations <- 100 
sampleUpdateIterations <- 300
cool = 10 
cvTargetX = 0.07
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


# plot and summary
plot(saStrata)
summary(saStrata)

