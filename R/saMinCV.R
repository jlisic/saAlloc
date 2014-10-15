saMinCV <-
function(
  x,
  label,
  targetCV,
  sampleSize,
  weightMatrix,            # missing handled
  iterations=1000,
  sampleIterations=100,
  cooling=0,
  segments=rep(1,nrow(x)),
  PSUAcres=rep(1,nrow(x)),
  targetVarWithin=rep(0,ncol(x))
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
  acceptRate <- rep(0,3*iterations) 
  cost <- rep(1,d+1) 
  
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
  total <- colSums(x)
  
  
  #################### OPTIMAL SAMPLE SIZE ######################################
  
  # handle sample size 
  if( length(sampleSize) > 1 ) {
    if( H != length(sampleSize) ) {
      stop( "H != length(sampleSize)" ) 
    }
  } else if(length(sampleSize) == 1) {

    sampleSize.n <- sampleSize # make a copy of sampleSize

    sampleSize <- rep(floor(sampleSize/H) , H)
    print(sampleSize)

    if( sum(sampleSize) < sampleSize.n ) {
      sampleSizeAddTo <- sample(1:H,size=sampleSize.n - sum(sampleSize))
      sampleSize[ sampleSizeAddTo ] <- sampleSize[ sampleSizeAddTo ] + 1  
    }
    sampleSize <- sapply(sampleSize, function(x) max(2,x) )

  } else {
      stop( "0 == length(sampleSize)" ) 
  }
    

  print("sampleSize")
  print(sampleSize)

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
  adminInt <- c(segments, sampleIterations) 
  adminIntLength <- length( adminInt )

  dup <- c() 


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
  as.double(sampleSize)       #checked      14
)

  print("C running time")
  print(proc.time() - Cprog) 
 

  #################### RETURN DATA ######################################
  
  a <- matrix(unlist(r.result[13]),ncol=3,byrow=T)
  colnames(a) <- c( 'change', 'U', 'accepted')

  rlabel <- sapply(unlist(r.result[6]), function(x) unique.label[x+1] ) 

  myList <- list("accept"=a, "cost"=unlist(r.result[7]), "label"=rlabel, "sampleSize"=unlist(r.result[14]) )
    
  return(myList)
}
