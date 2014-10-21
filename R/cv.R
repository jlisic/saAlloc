
# function to calculate cv
.cv <- function( x, strata = rep(1,nrow(x)), sampleSize=1, average=FALSE) {


  # get k 
  if(average) {
    k <- nrow(x) / label(x)
    if(!is.integer(k)) stop("average distance cannot be calculated, non-constant k")
  } else {
    k <- 1
  }

  # get variance
  if(k > 1) {
    vars <- 
      aggregate( 
        x,
        by= list( rep(strata,each=k), rep(1:k,length(strata))),
        var
      )

    vars <- aggregate( vars$x, by=list(vars$Group.1), mean)$x
    

  } else {
    # get variance
    vars <- aggregate( x,by=list(strata), var)$x
  }
 
  # get number of strata
  H <- length(unique(strata))                                                                                                                                               
    # if only one sample is provided, samples are distributed uniformly
  if( length(sampleSize) != H) {
       
    if( length(sampleSize) != 1 ) stop("Invalid number of sample sizes")
  
    sampleSize = rep(sampleSize/H,H)
  }
   
  # get total
  totals <- colSums( x ) /k
  Nh <- aggregate( x,by=list(strata), length)[,-1]
  return(sqrt(colSums(vars * Nh^2/sampleSize)) / totals )
}

