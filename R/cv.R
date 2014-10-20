
# function to calculate cv
.cv <- function( x, strata = rep(1,nrow(x)), sampleSize=1) {
   
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

