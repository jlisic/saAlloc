# function to perform simulated annealing based sample size allocation
saSampleAlloc <-
function(
  x,
  label,
  targetCV,
  sampleSize,
  weightMatrix,            # missing handled
#  domainMatrixList,
  locationAdjustment,
  scaleAdjustment,
  sampleSizeIterations=100,
  p = 2,                   # l2 norm of the penalty function
  penalty = -1,            # negative penalties are ignored
  cooling = 0,
  preserveSatisfied=TRUE
  ) {

  tolSize <- 1 # not used


  
  #################### DATA SET ######################################

  # check if X is a matrix
  if( !is.matrix(x) ) {
    stop("x must be a matrix") 
  }
 
  # get rows etc... 
  N <- length(label)            # number of observations
  K <- ncol(x)                  # number of distinct PSUs 
  R <- length( c(x) ) / ( N * K)  # number of observations of each PSU
  
  
  # transform label  
  unique.label <- sort(unique(label))
  rlabel <- sapply(label, function(x){ which(x == unique.label) } ) 
  rlabel <- rlabel - 1


  H <- length(unique.label)    # number of strata

  #################### PROBABILITY ######################################



  # get prob
  if( missing(weightMatrix) ) {
    weightMatrix <- matrix( 1/N, nrow=N, ncol=H)
  }

  # check prob
  if(( nrow(weightMatrix) != N) | (ncol(weightMatrix) != H )) {
    stop( sprintf("probMatrixability matrix (weightMatrix) does not have correct dimensions\n Needed (nrow = %d, ncol = %d), provided (nrow = %d, ncol = %d)\n",
                 N, K, nrow(weightMatrix), ncol(weightMatrix) ) ) 
  }

  # get max prob
  prob <- 1 - apply(cbind(rlabel,weightMatrix), 1, function(x) x[x[1]+2] ) 

  totalProb <- sum(prob)
  
  # create row major matrix for input to C program
  weightMatrix <- c(t(weightMatrix))

  # get colnames for x 
  if( is.null(colnames(x)) ) {
    x.colnames <- sprintf("%d",1:K)
  } else {
   x.colnames <- colnames(x)
  }

  #################### LOCATION ADJUSTMENT ######################################
 
  # location adjustment is N x K  
  if( missing( locationAdjustment ) ) {
    locationAdjustment = -1
  } else {
    locationAdjustment = c(locationAdjustment) 
  }
 

  #################### SCALE ADJUSTMENT #########################################
 
  # location adjustment is H x K  
  if( missing( scaleAdjustment ) ) {
    scaleAdjustment = -1
  } else {
    scaleAdjustment = c(scaleAdjustment) 
  }
  
  
  #################### DOMAIN MATRIX ######################################
  # it is assumed that the domain matrix is a list of matricies of dim K x J of length H
  # where K is the number of coi, 
  # J is the number of domains 
  # and H is the number of strata
  # if it is missing this list will be substituted with a list of identity matricies.

  #if( missing(domainMatrixList) ) {
    J <- K
    domainMatrixList.names <- x.colnames
    domainMatrix <- rep(c(diag(K)),H)
#  } else {
#    domainMatrixList.names <- names(domainMatrixList) 
#    J <- ncol(domainMatrixList[[1]])
#    if( J < K) stop("Error J < K")
#    if( J == K) {
#      domainMatrix <- rep(c(diag(K)),H)
#    } else {
#      domainMatrixTmp <- c()
#      for( h in 1:H ) domainMatrixTmp <- c(domainMatrixTmp, c( t(domainMatrixList[[h]]) ))
#      domainMatrix <- domainMatrixTmp 
#    }
#
#  }

  #################### PENALTY ######################################

  if( sum( penalty < 0) > 1 ) {
    penalty <- rep(0,J)
  } else {

    if( length(penalty) == 1 )  {
      penalty <- rep(0,J)
    } else {
      if( length(penalty) != J ) {
        stop("penalty is not of length 1 or J")
      } 
    }

  }
  
  
  #################### SAMPLE ITERATIONS ######################################
  
  if(sampleSizeIterations < 0 ) stop("sampleSizeIterations is not positive") 
  if( round(sampleSizeIterations) != sampleSizeIterations ) stop("sampleSizeIterations is not an integer ") 
   
  
  #################### TARGET CV ######################################
  # handle sample size 

#  # stop if the CV length does not equal the number of domains 
#  if( d != length(targetCV) ) {
#    stop( "ncol(x) != length(targetCV)" )  
#  # warning if there are no names 
#  } else if( length(names(targetCV)) == 0 ) {
#    warning( "targetCV has no names, assuming that the targetCV are in the same order as x." ) 
#  # error if there are names but don't match what we find in unique.label
#  } else if( !identical( sort(names(targetCV)), sort(colnames(x)))  ) {
#    stop( "targetCV names do match column names of x" ) 
#  } else if( anyDuplicated( names(targetCV) ) != 0  ) {
#    stop( "targetCV names do exist but are not unique" ) 
#  # at this point everthing seems ok, so reorder the sample size
#  } else {
#    targetCV <- targetCV[ colnames(x) ] 
#  }
  
  
  #################### INITIAL SAMPLE SIZE ######################################
  
  # handle sample size 
  if( length(sampleSize) > 1 ) {

    # stop if the sampleSize length does not equal the number of strata 
    if( H != length(sampleSize) ) {
      stop( "H != length(sampleSize)" )  
    # warning if there are no names 
    } else if( length(names(sampleSize)) == 0 ) {
      #warning( "sampleSize has no names, assuming that the sample sizes are in asscending order with respect to the strata in label." ) 
    # error if there are names but don't match what we find in unique.label
    } else if( !identical( sort(as.numeric(names(sampleSize))), unique.label)  ) {
      stop( "sampleSize names do match elements in label" ) 
    # at this point everthing seems ok, so reorder the sample size
    } else {
      sampleSize <- sampleSize[ as.character( unique.label ) ] 
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
  a <- rep(0,sampleSizeIterations)

  #################### RUN C FUNCTION ######################################
  
Cprog <- proc.time()
  
r.result <- .C("R_sampleAlloc",
  as.double(c(x)),                  # 1 
  as.integer(N),                    # 2 
  as.integer(J),                    # 2
  as.integer(K),                    # 4 
  as.integer(H),                    # 5
  as.integer(R),                    # 6
  as.integer(sampleSizeIterations), # 7     
  as.integer(rlabel),               # 8 
  as.integer(domainMatrix),         # 9
  as.double(weightMatrix),          # 10
  as.double(targetCV),              # 11
  as.double(locationAdjustment),    # 12
  as.double(scaleAdjustment),       # 13
  as.double(p),                     # 14
  as.double(penalty),               # 15
  as.double(sampleSize),            # 16
  as.double(a),                     # 17
  as.double(cooling)                # 18
)

runTime <- proc.time() - Cprog
 

  #################### RETURN DATA ######################################


  ## label
  newRlabel <- sapply(unlist(r.result[8]), function(x) unique.label[x+1] ) 
   
  ## sample size 
  newSampleSize <- unlist(r.result[16]) 
  names(newSampleSize) <- sprintf("n_%d",unique.label)
  names(sampleSize) <- sprintf("n_%d",unique.label)

  ## final and initial CVs
  CVStart <- .cv( x, rlabel, sampleSize, average=TRUE) 
  CV      <- .cv( x, newRlabel, newSampleSize, average=TRUE) 

  ## strata Size
  strataSizeStart <- stats::aggregate(rlabel, by=list(rlabel), length)
  rownames(strataSizeStart) <- strataSizeStart$Group.1
  strataSizeStart <- strataSizeStart[,2,drop=FALSE] 

  strataSize      <- stats::aggregate(newRlabel, by=list(newRlabel), length)
  rownames(strataSize) <- strataSize$Group.1
  strataSize <- strataSize[,2,drop=FALSE] 



  myList <- list(
    "objective"=r.result[[17]],
    "label"=newRlabel,
    "sampleSize"=newSampleSize,
    "sampleSizeStart"=sampleSize,
    "CV"=CV,
    "CVStart"=CVStart,
    "targetCV"=targetCV,
    "strataSize"=strataSize,
    "strataSizeStart"=strataSizeStart,
    "runTime"=runTime,
    "variables"=x.colnames
  )


  class(myList) <- "saAlloc"

  return(myList)
}






























