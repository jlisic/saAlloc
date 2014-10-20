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
iterations <- 50000
sampleIterations <- 0
sampleUpdateIterations <- 300
cool = 0 
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


print( cbind( cv.minCV, cv.kMeansCV, c(cvTargetX, cvTargetY)))



# plot for s3 saAlloc object
# together == T, all vars on one plot
# together == F, separate plots for each commodity
plot.saAlloc <- function( x, together=TRUE ) {
  
  # get variables 
  accept <- x$accept
  variables <- x$variables

  ### pars for multiplots
  if(!together ) {
    plotRow <-  ceiling(sqrt(length(variables)))
    plotCol <-  ceiling(length(variables) / plotRow)
    par(mfrow= c(plotRow,plotCol) )
  }


  #### plot all variables on the same plot

  # plot first variable 
  plot(
    1:iterations,
    accept[,variables[1]],
    type='l', 
    ylim=c(0,max(accept[,variables])) 
  )

  # add title
  if(together) {
    title("Trace of CVs")
  } else {
    title(sprintf("Trace of CV for %s",variables[1]))
  }

  # plot other variables 
  if(length(variables) > 1 ) {

    if( together ) {
      for( i in 2:length(variables) ) {
        points(
            1:iterations,
            accept[,variables[i] ],
            type='l',
            col=i
          )
      }
    } else {
      for( i in 2:length(variables) ) {
        plot(
          1:iterations,
          accept[,variables[i] ],
          type='l',
          col=i
        )
        title(sprintf("Trace of CV for %s",variables[i]))
      }
    }


  }

  # add legend
  if( together ) {
    legend( 
      "topright",
      variables,                     # add names
      lty=rep(1,length(variables)),  # add lines
      col = 1:length(variables)      # add color
    )
  } 
}



# summary for s3 saAlloc object
summary.saAlloc <- function(x) {

  print("CVs")
  y <- cbind( x$CVStart, x$CV, x$targetCV ) 
  colnames(y) <- c("Initial", "Final", "Target" ) 
  print(y)
 
  print("Sample Size") 
  y <- cbind( x$sampleSizeStart, x$sampleSize ) 
  colnames(y) <- c("Initial", "Final" ) 
  print(y)


  print("Strata Size:")
  y <- cbind( x$strataSizeStart, x$strataSize ) 
  colnames(y) <- c("Initial", "Final" ) 
  print(y)


  print("Auxiliary Size Constraint:")
  y <- cbind( x$acresStart, x$acres ) 
  colnames(y) <- c("Initial", "Final" ) 
  print(y)


  print( x$runTime)

}





# plot and summary
plot(saStrata)
summary(saStrata)

