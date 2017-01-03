#
## function to calculate cv
#.cv <- function( x, strata = rep(1,nrow(x)), sampleSize=1, average=FALSE) {
#
#
#  # get k 
#  if(average) {
#    k <- nrow(x) / length(strata)
#    if(as.integer(k) != k) stop(sprintf("average distance cannot be calculated, non-integer k=%f", k))
#  } else {
#    k <- 1
#  }
#
#  # get variance
#  if(k > 1) {
#    vars <- 
#      stats::aggregate( 
#        x,
#        by= list( rep(strata,each=k), rep(1:k,length(strata))),
#        stats::var
#      )
#    
#    vars <- stats::aggregate( 
#              vars[,-c(1,2),drop=FALSE], 
#              by=list(vars$Group.1), 
#              mean)[,-1]
#    
#
#  } else {
#    # get variance
#    vars <- stats::aggregate( x,by=list(strata), stats::var)[,-1]
#  }
# 
#  # get number of strata
#  H <- length(unique(strata))                                                                                                                                               
#    # if only one sample is provided, samples are distributed uniformly
#  if( length(sampleSize) != H) {
#       
#    if( length(sampleSize) != 1 ) stop("Invalid number of sample sizes")
#  
#    sampleSize = rep(sampleSize/H,H)
#  }
#   
#  # get total
#  totals <- colSums( x ) /k
#  Nh <- stats::aggregate( strata,by=list(strata), length)[,-1]
#
#
#  return(sqrt(colSums(vars * Nh^2/sampleSize)) / totals )
#}



.cv2 <- function( x, strata = rep(1,nrow(x)), sampleSize=1, average=FALSE, locationAdjustment, scaleAdjustment, fpc=TRUE) {


  # get R 
  if(average) {
    R <- nrow(x) / length(strata)
    if(as.integer(R) != R) stop(sprintf("average distance cannot be calculated, non-integer R=%f", R))
  } else {
    R <- 1
  }

  vars <- 
    stats::aggregate( 
      x,
      by= list( rep(strata,each=R), rep(1:R,length(strata))),
      stats::var
    )[,c(-1,-2)]

 vars.init <<- vars 
 R <<- R

  # perform scale adjustment on S^2 
  if( !missing(scaleAdjustment) ) {
    scaleAdjustment <- 
    stats::aggregate( 
      scaleAdjustment,
      by= list( rep(strata,each=R), rep(1:R,length(strata))),
      mean 
    )[,c(-1,-2)]
    vars = vars * scaleAdjustment 
  }
  
  # perform location adjustment on S^2 
  if( !missing(locationAdjustment) ) {
      list( rep(strata,each=R), rep(1:R,length(strata)))

    locationAdjustment <- 
    stats::aggregate( 
      locationAdjustment,
      by= list( rep(strata,each=R), rep(1:R,length(strata))),
      mean 
    )[,c(-1,-2)]
    vars = vars + locationAdjustment 
  }

  # get number of strata
  H <- length(unique(strata))                                                                                                                                               
  # if only one sample is provided, samples are distributed uniformly
  if( length(sampleSize) != H) {
    if( length(sampleSize) != 1 ) stop("Invalid number of sample sizes")
    sampleSize = rep(sampleSize/H,H)
  }
 
  # get total
  totals <- stats::aggregate(x, by=list(rep(1:R,length(strata))),sum)[,-1] 
  Nh <- stats::aggregate( strata,by=list(strata), length)[,-1]

  vars <<- vars
  totals <<- totals
  Nh <<- totals


  if(fpc) {
    varsAdj <- vars*(Nh^2/sampleSize) * (1 - sampleSize/Nh)
  } else {
    varsAdj <- vars*(Nh^2/sampleSize)
  }

  result <- sqrt( stats::aggregate( varsAdj , by=list(rep(1:R,each=H)), sum)[,-1])  

  varsAdj <<- varsAdj 

  if( !is.null(dim(result)) ) {
    return( colMeans( result / totals ) )
  } else {
    return( mean( result / totals ) )
  }
}




