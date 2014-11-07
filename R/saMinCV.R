saMinCV <-
function(
  x,
  label,
  targetCV,
  sampleSize,
  weightMatrix,            # missing handled
  iterations=1000,
  sampleIterations=100,
  sampleUpdateIterations=100,
  cooling=0,
  segments=rep(1,length(label)),
  PSUAcres=rep(1,length(label)),
  targetVarWithin=rep(0,ncol(x)),
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
  d <- ncol(x)                  # number of distinct characteristics  (independent subsets of the covariance matrix)
  k <- length(c(x)) / ( N * d)  # dimension of each of the distinct characteristics 
                                # (dependent members of the subsets of the covariance matrix)
  H <- length(unique.label)    # number of strata

  # need to re-label
  acceptRate <- rep(0,(5 + H + d)*iterations) 
  cost <- rep(1,d) 
  
  #################### PROBABILITY ######################################
 
  # get prob
  if( missing(weightMatrix) ) {
    weightMatrix <- matrix( 1/N, nrow=N, ncol=H)
  }

  # check prob
  if(( nrow(weightMatrix) != N) | (ncol(weightMatrix) != H )) {
    stop( sprintf("probMatrixability matrix (weightMatrix) does not have correct dimensions\n Needed (nrow = %d, ncol = %d), provided (nrow = %d, ncol = %d)\n",
                 N, d, nrow(weightMatrix), ncol(weightMatrix) ) ) 
  }

  # get max prob
  prob <- 1 - apply(cbind(rlabel,weightMatrix), 1, function(x) x[x[1]+2] ) 

  totalProb <- sum(prob)
  
  # create row major matrix for input to C program
  weightMatrix <- c(t(weightMatrix))


  #################### TOTALS ######################################

  # get total 
  total <- colSums(x)/k
  
  #################### SAMPLE SIZE ######################################
  
  if(sampleUpdateIterations < 0 ) stop("sampleUpdateIter is not positive") 
  if( round(sampleUpdateIterations) != sampleUpdateIterations ) stop("sampleUpdateIter is not an integer ") 
   
  
  #################### TARGET CV ######################################
  # handle sample size 

  # stop if the sampleSize lenght does not equal the number of distinct elements
  if( d != length(targetCV) ) {
    stop( "ncol(x) != length(targetCV)" )  
  # warning if there are no names 
  } else if( length(names(targetCV)) == 0 ) {
    warning( "targetCV has no names, assuming that the targetCV are in the same order as x." ) 
  # error if there are names but don't match what we find in unique.label
  } else if( !identical( sort(names(targetCV)), sort(colnames(x)))  ) {
    stop( "targetCV names do match column names of x" ) 
  } else if( anyDuplicated( names(targetCV) ) != 0  ) {
    stop( "targetCV names do exist but are not unique" ) 
  # at this point everthing seems ok, so reorder the sample size
  } else {
    targetCV <- targetCV[ colnames(x) ] 
  }
  
  
  
  #################### OPTIMAL SAMPLE SIZE ######################################
  
  # handle sample size 
  if( length(sampleSize) > 1 ) {

    # stop if the sampleSize lenght does not equal the number of distinct elements
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


  # group data together for input
  adminDbl <- c( 
    PSUAcres, 
    targetCV, 
    targetVarWithin, 
    total, 
    sampleSize, 
    prob, 
    weightMatrix, 
    totalProb,
    cooling, 
    tolSize 
  )
  adminDblLength <- length( adminDbl )
  adminInt <- c(segments, sampleIterations, as.integer(preserveSatisfied)) 
  adminIntLength <- length( adminInt )

  dup <- c() 

  costChangeSize <- 5 + H + d 


  #################### RUN C FUNCTION ######################################
  
r.result <- .C("R_minCV",
  as.double(c(x)),            #checked       1
  as.integer(k),              #checked       2
  as.integer(d),              #checked       3
  as.integer(N),              #checked       4
  as.integer(iterations),     #checked       5
  as.integer(rlabel),          #checked       6
  as.double(cost),            #checked Q     7
  as.double(adminDbl),        #checked       8
  as.integer(adminInt),       #checked       9
  as.integer(adminIntLength), #checked      10
  as.integer(adminDblLength), #checked      11
  as.integer(dup),            #checked      12
  as.double(acceptRate),      #checked      13
  as.double(sampleSize),      #checked      14
  as.integer(sampleUpdateIterations),
  as.integer(costChangeSize)
)

  runTime <- proc.time() - Cprog
 

  #################### RETURN DATA ######################################

  # get colnames for x 
  if( is.null(colnames(x)) ) {
    x.colnames <- sprintf("%d",1:d)
  } else {
    x.colnames <- colnames(x)
  }

  ## change and other accept data
  a <- matrix(unlist(r.result[13]),ncol=5 + H + d,byrow=T)
  colnames(a) <- c( 'change', 'U', 'accepted','from','to', 
                   sprintf("n_%d",unique.label), 
                   x.colnames
                   )
  ## label
  newRlabel <- sapply(unlist(r.result[6]), function(x) unique.label[x+1] ) 
   
  ## sample size 
  newSampleSize <- unlist(r.result[14]) 
  names(newSampleSize) <- sprintf("n_%d",unique.label)
  names(sampleSize) <- sprintf("n_%d",unique.label)

  ## final and initial CVs
  CVStart <- .cv( x, rlabel, sampleSize, average=TRUE) 
  CV      <- .cv( x, newRlabel, newSampleSize, average=TRUE) 

  ## strata Size
  strataSizeStart <- aggregate(rlabel, by=list(rlabel), length)
  rownames(strataSizeStart) <- strataSizeStart$Group.1
  strataSizeStart <- strataSizeStart[,2,drop=FALSE] 

  strataSize      <- aggregate(newRlabel, by=list(newRlabel), length)
  rownames(strataSize) <- strataSize$Group.1
  strataSize <- strataSize[,2,drop=FALSE] 


  ## auxiliary size constraint (acres)
  acresStart <- aggregate(PSUAcres*segments, by=list(rlabel), sum)
  rownames(acresStart) <- acresStart$Group.1
  acresStart <- acresStart[,2,drop=FALSE] 

  acres      <- aggregate(PSUAcres*segments, by=list(newRlabel), sum)
  rownames(acres) <- acres$Group.1
  acres <- acres[,2,drop=FALSE] 



  myList <- list(
    "accept"=a, 
    "cost"=unlist(r.result[7]), 
    "label"=newRlabel, 
    "sampleSize"=newSampleSize,
    "sampleSizeStart"=sampleSize,
    "CV"=CV,
    "CVStart"=CVStart,
    "targetCV"=targetCV,
    "strataSize"=strataSize,
    "strataSizeStart"=strataSizeStart,
    "acres"=acres,
    "acresStart"=acresStart,
    "runTime"=runTime,
    "variables"=x.colnames
  )


  class(myList) <- "saAlloc"

  return(myList)
}
