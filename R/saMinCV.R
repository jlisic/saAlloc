saMinCV <- function(
  x,
  label,
  targetCV,
  sampleSize,
  weightMatrix,            # missing handled
  domainMatrixList,
  iterations=1000,
  sampleSizeIterations=3,
  recalcIterations=10000,
  locationAdjustment,
  scaleAdjustment,
  p = 2,                   # l2 norm of the penalty function
  penalty = -1,            # negative penalties are ignored
  cooling = 0,
  preserveSatisfied=TRUE
  ) {

  tolSize <- 1 # not used


  # check if X is a matrix
  if( !is.matrix(x) ) {
    stop("x must be a matrix") 
  }
  
  Cprog <- proc.time()
  
  # transform label  
  unique.label <- sort(unique(label))
  rlabel <- sapply(label, function(x){ which(x == unique.label) } ) 
  rlabel <- rlabel - 1


  # get rows etc... 
  N <- length(label)            # number of observations
  K <- ncol(x)                  # number of distinct characteristics  (independent subsets of the covariance matrix)
  R <- length(c(x)) / (N * K)  # dimension of each of the distinct characteristics 
                                # (dependent members of the subsets of the covariance matrix)
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
  prob <- apply(cbind(rlabel,weightMatrix), 1, function(x) max(x[-1* c(1, x[1]+2)]) ) 

  totalProb <- sum(prob)
  
  # create row major matrix for input to C program
  weightMatrix <- c(t(weightMatrix))
  
  #################### LOCATION ADJUSTMENT ######################################
 
  # location adjustment is N x K  
  if( missing( locationAdjustment ) ) {
    C_locationAdjustment = -1
  } else {
    C_locationAdjustment = c(locationAdjustment) 
  }
 

  #################### SCALE ADJUSTMENT #########################################
 
  # location adjustment is H x K  
  if( missing( scaleAdjustment ) ) {
    C_scaleAdjustment = -1
  } else {
    C_scaleAdjustment = c(scaleAdjustment) 
  }
  
  
  #################### DOMAIN MATRIX ######################################
  # it is assumed that the domain matrix is a list of matricies of dim K x J of length H
  # where K is the number of coi, 
  # J is the number of domains 
  # and H is the number of strata
  # if it is missing this list will be substituted with a list of identity matricies.

  if( missing(domainMatrixList) ) {
    J <- K
    domainMatrix <- rep(c(diag(K)),H)
  } else {
    J <- ncol(domainMatrixList[[1]])
    if( J < K) stop("Error J < K")
    if( J == K) {
      domainMatrix <- rep(c(diag(K)),H)
    } else {
      domainMatrixTmp <- c()
      for( h in 1:H ) domainMatrixTmp <- c(domainMatrixTmp, c( t(domainMatrixList[[h]]) ))
      domainMatrix <- domainMatrixTmp 
    }

  }

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
  
  #################### INITIAL SAMPLE SIZE ######################################
  
  # handle sample size 
  if( length(sampleSize) > 1 ) {

    # stop if the sampleSize length does not equal the number of strata 
    if( H != length(sampleSize) ) {
      stop( "H != length(sampleSize)" )  
    # warning if there are no names 
    } else if( length(names(sampleSize)) == 0 ) {
      warning( "sampleSize has no names, assuming that the sample sizes are in asscending order with respect to the strata in label." ) 
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
  
  
  #################### ITERATIONS ######################################
  
  if(iterations < 0 ) stop("Iterations is not positive") 
  if( round(iterations) != iterations ) stop("Iterations is not an integer ") 
 

  #################### SAMPLE ITERATIONS ######################################
  
  if(sampleSizeIterations < 0 ) stop("sampleSizeIterations is not positive") 
  if( round(sampleSizeIterations) != sampleSizeIterations ) stop("sampleSizeIterations is not an integer ") 
 

  #################### SAMPLE SIZE ######################################
  
  if(recalcIterations < 0 ) stop("recalcIterations is not positive") 
  if( round(recalcIterations) != recalcIterations ) stop("recalcIterations is not an integer ") 
   
  
  #################### TARGET CV ######################################
  # handle sample size 

  # stop if the sampleSize lenght does not equal the number of distinct elements
  if( J != length(targetCV) ) {
    stop( "ncol(x) != length(targetCV)" )  
  # warning if there are no names 
  } else if( length(names(targetCV)) == 0 ) {
    warning( "targetCV has no names, assuming that the targetCV are in the same order as x." ) 
  # error if there are names but don't match what we find in unique.label
  } else if( !identical( sort(names(targetCV)), sort(names(domainMatrix)))  ) {
    stop( "targetCV names do match column names of x" ) 
  } else if( anyDuplicated( names(targetCV) ) != 0  ) {
    stop( "targetCV names do exist but are not unique" ) 
  # at this point everthing seems ok, so reorder the sample size
  } else {
    targetCV <- targetCV[ names(domainMatrix) ] 
  }
  
  costChangeSize <- 6 + H + J  
  
  # need to re-label
  acceptRate <- rep(0,costChangeSize*(iterations+1)) 

#  print(paste( "N", N))
#  print(paste( "l(x)", length(c(x))))
#  print(paste( "K", K))
#  print(paste( "R", R))
#  print(paste( "H", H))
#  print(paste( "penalty: ", penalty))


  #################### RUN C FUNCTION ######################################
  r.result <- .C("R_minCV",
    as.double(c(x)),                  # 1 
    as.integer(N),                    # 2 
    as.integer(J),                    # 2
    as.integer(K),                    # 4 
    as.integer(H),                    # 5
    as.integer(R),                    # 6
    as.integer(sampleSizeIterations), # 7     
    as.integer(rlabel),               # 8 
    as.integer(domainMatrix),         # 9
    as.double(prob),                  #10
    as.double(weightMatrix),          #11
    as.double(targetCV),              #12
    as.double(C_locationAdjustment),    #13
    as.double(C_scaleAdjustment),       #14
    as.double(p),                     #15
    as.double(penalty),               #16
    as.double(sampleSize),            #17
    as.double(acceptRate),            #18
    as.double(cooling),                  #19
    as.integer(recalcIterations),     #20  /* how often to update the counter */
    as.integer(costChangeSize),       #21
    as.integer(iterations),           #22     
    as.integer(preserveSatisfied)     #23
  )

  runTime <- proc.time() - Cprog
 

  #################### RETURN DATA ######################################

  # get colnames for x 
  if( is.null(colnames(x)) ) {
    x.colnames <- sprintf("%d",1:K)
  } else {
    x.colnames <- colnames(x)
  }

  a <<- unlist(r.result[18])

  ## change and other accept data
  a <- matrix(unlist(r.result[18]),ncol=costChangeSize ,byrow=T)
  colnames(a) <- c( 'change', 'U', 'T','selected','from','to', 
                   sprintf("n_%d",unique.label), 
                   x.colnames
                   )

  
  a.names <- c(colnames(a),'accepted')
  a <- cbind(a, as.numeric(a[,'U'] <= a[,'T']))
  colnames(a) <- a.names


  ## label
  newRlabel <- sapply(unlist(r.result[8]), function(x) unique.label[x+1] ) 
   
  ## sample size 
  newSampleSize <- unlist(r.result[17]) 
  names(sampleSize) <- sprintf("n_%d",unique.label)
  names(newSampleSize) <- sprintf("n_%d",unique.label)

  ## final and initial CVs
  if( missing(locationAdjustment) & missing(scaleAdjustment) ) { 
    print("no adjustment")
    CVStart  <- .cv2( x, rlabel, newSampleSize, average=TRUE)
    CV  <- .cv2( x, newRlabel, newSampleSize, average=TRUE )
  } else if( missing(locationAdjustment) )  { 
    print("just scale adjustment")
    CVStart  <- .cv2( x, rlabel, newSampleSize, average=TRUE, scaleAdjustment=scaleAdjustment)
    CV  <- .cv2( x, newRlabel, newSampleSize, average=TRUE, scaleAdjustment=scaleAdjustment)
  } else if( missing(scaleAdjustment) )  { 
    print("just location adjustment")
    CVStart  <- .cv2( x, rlabel, newSampleSize, average=TRUE, locationAdjustment=locationAdjustment)
    CV  <- .cv2( x, newRlabel, newSampleSize, average=TRUE, locationAdjustment=locationAdjustment)
  } else {
    print("both adjustments")
    CVStart  <- .cv2( x, rlabel, newSampleSize, average=TRUE, locationAdjustment=locationAdjustment, scaleAdjustment=scaleAdjustment)
    CV  <- .cv2( x, newRlabel, newSampleSize, average=TRUE, locationAdjustment=locationAdjustment, scaleAdjustment=scaleAdjustment)
  }


  ## strata Size
  strataSizeStart <- aggregate(rlabel, by=list(rlabel), length)
  rownames(strataSizeStart) <- strataSizeStart$Group.1
  strataSizeStart <- strataSizeStart[,2,drop=FALSE] 

  strataSize      <- aggregate(newRlabel, by=list(newRlabel), length)
  rownames(strataSize) <- strataSize$Group.1
  strataSize <- strataSize[,2,drop=FALSE] 


  myList <- list(
    "accept"=a, 
#    "cost"=unlist(r.result[7]), 
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
