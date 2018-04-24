# function to perform simulated annealing based sample size allocation
saSampleAlloc <-
function(
  total,
  S2,
  targetCV,
  sampleSize,
  strataSize, # size of the stratum
  locationAdjustment,
  scaleAdjustment,
  iterations=100,
  p = 2,                   # l2 norm of the penalty function
  penalty = -1,            # negative penalties are ignored
  cooling = 0,
  preserveSatisfied=TRUE,
  fpc = TRUE
  ) {
 
  tolSize <- 1 # not used

  # get the number of variables
  J <- length(targetCV)

  # get the number of strata
  H <- length(strataSize)

  # get the number of reps
  R <- NROW(S2) / H 

  if( round(H) != H) stop("Number of rows of total are not an integer") 
  if( NCOL(total) != NCOL(S2) ) stop("Number of columns between total and S2 do not match") 

  #################### LOCATION ADJUSTMENT ######################################
 
  # location adjustment is N x J  
  if( missing( locationAdjustment ) ) {
    locationAdjustment = -1
  } else {
    if( !identical(dim(locationAdjustment),dim(S2)) ) stop("Dimensions between locationAdjustment and S2 do not match") 
    locationAdjustment = c(locationAdjustment) 
  }
 

  #################### SCALE ADJUSTMENT #########################################
 
  # location adjustment is H x J  
  if( missing( scaleAdjustment ) ) {
    scaleAdjustment = -1
  } else {
    if( !identical(dim(scaleAdjustment),dim(S2)) ) stop("Dimensions between scaleAdjustment and S2 do not match") 
    scaleAdjustment = c(scaleAdjustment) 
  }
  
  
  #################### DOMAIN MATRIX ######################################
  # not used
    domainMatrix <- rep(c(diag(J)),H)

  #################### PENALTY ######################################

  if( sum( penalty < 0) >= 1 ) {
    penalty <- rep(0,J)
  } else {

    if( length(penalty) == 1 )  {
      penalty <- rep(penalty,J)
    } else {
      if( length(penalty) != J ) {
        stop("penalty is not of length 1 or J")
      } 
    }

  }
  
  
  #################### SAMPLE ITERATIONS ######################################
  
  if(iterations < 0 ) stop("iterations is not positive") 
  if( round(iterations) != iterations ) stop("iterations is not an integer ") 
  

  
  
  #################### INITIAL SAMPLE SIZE ######################################
  
  # handle sample size 
  if( length(sampleSize) > 1 ) {

    # stop if the sampleSize length does not equal the number of strata 
    if( H != length(sampleSize) ) {
      stop( "H != length(sampleSize)" )  
    # warning if there are no names 
    } 

  } else if(length(sampleSize) == 1) {

    sampleSize.n <- sampleSize # make a copy of sampleSize

    sampleSize <- rep(floor(sampleSize/H) , H)

    if( sum(sampleSize) < sampleSize.n ) {
      sampleSizeAddTo <- sample(1:H,size=sampleSize.n - sum(sampleSize))
      sampleSize[ sampleSizeAddTo ] <- sampleSize[ sampleSizeAddTo ] + 1  
    }
    sampleSize <- sapply(sampleSize, function(x) max(2,x) )

  } else {
      stop( "0 == length(sampleSize)" ) 
  }
  # final sample size check 
  if( sum(sampleSize < 2) > 0 ) stop("Each element of sampleSize must be of at least 2")


  ### allocate space to return objective functions
  a <- rep(0,iterations)

  #################### RUN C FUNCTION ######################################
  
Cprog <- proc.time()
  
r.result <- .C("R_sampleAlloc",
  as.double(total),                 # 1 the column major listing of totals J*H*R 
  as.double(S2),                    # 2 the column major listing of variances J*H*R 
  as.integer(J),                    # 5 number of variables 
  as.integer(H),                    # 6 number of strata 
  as.integer(R),                    # 7 number of observatinos of each PSU (e.g. years) 
  as.integer(iterations), # 8 number of iterations     
  as.integer(domainMatrix),         # 10 not really used right now
  as.double(targetCV),              # 11 target CV
  as.double(locationAdjustment),    # 12 location adjustment J*H*R
  as.double(scaleAdjustment),       # 13 scale adjustment J*H*R
  as.double(p),                     # 14 exponent
  as.double(penalty),               # 15 lambda (penalty value, length == J)
  as.double(sampleSize),            # 16 sample size by stratum
  as.integer(strataSize),            # 17 stratum size
  as.double(a),                     # 18 what we return
  as.double(cooling),               # 19 cooling
  as.integer(fpc)                   # 20 fpc (1 = yes, 0 = no)
)


runTime <- proc.time() - Cprog
 

  #################### RETURN DATA ######################################


#  ## sample size 
#  newSampleSize <- unlist(r.result[16]) 
#  names(newSampleSize) <- sprintf("n_%d",unique.label)
#  names(sampleSize) <- sprintf("n_%d",unique.label)
#
#  ## final and initial CVs
#  CVStart <- .cv2( x, rlabel, sampleSize, average=TRUE,fpc=fpc) 
#  CV      <- .cv2( x, newRlabel, newSampleSize, average=TRUE,fpc=fpc) 
#
#  ## strata Size
#  strataSizeStart <- stats::aggregate(rlabel, by=list(rlabel), length)
#  rownames(strataSizeStart) <- strataSizeStart$Group.1
#  strataSizeStart <- strataSizeStart[,2,drop=FALSE] 
#
#  strataSize      <- stats::aggregate(newRlabel, by=list(newRlabel), length)
#  rownames(strataSize) <- strataSize$Group.1
#  strataSize <- strataSize[,2,drop=FALSE] 
#
#
#
#  myList <- list(
#    "objective"=r.result[[17]],
#    "label"=newRlabel,
#    "sampleSize"=newSampleSize,
#    "sampleSizeStart"=sampleSize,
#    "CV"=CV,
#    "CVStart"=CVStart,
#    "targetCV"=targetCV,
#    "strataSize"=strataSize,
#    "strataSizeStart"=strataSizeStart,
#    "runTime"=runTime,
#    "variables"=x.colnames,
#    "fpc"=TRUE
#  )
#
#
#  class(myList) <- "saAlloc"
#
#  return(myList)

  return(r.result)
}






























