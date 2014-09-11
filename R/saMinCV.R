saMinCV <-
function(
  x,
  label,
  targetCV,
  sampleSize,
  prob,                  # missing handled
  neighbors,             # missing handled
  nNeighbors,            # missing handled
  iterations=1000,
  cooling=0,
  segments=rep(1,nrow(x)),
  PSUAcres=rep(1,nrow(x)),
  targetVarWithin=rep(0,ncol(x)),
  tolSize=1
  ) {

  # dependency on Fast Nearest Neighbors
  requires(FNN)

  # check if X is a matrix
  if( !is.matrix(x) ) {
    stop("x must be a matrix") 
  }
  
  Cprog <- proc.time()
 
  # get rows etc... 
  rlabel <- label
  N <- length(label)            # number of observations
  d <- ncol(x)                  # number of distinct characteristics  (independent subsets of the covariance matrix)
  k <- length(c(x)) / ( N * d)  # dimension of each of the distinct characteristics (dependent members of the subsets of the covariance matrix)

  H <- length(unique(label))    # number of strata

  # need to re-label

  acceptRate <- rep(0,3*iterations) 
  cost <- rep(1,d+1) 
  
  
  #################### PROBABILITY ######################################
 
  # get prob
  if( is.missing(prob) ) {
    prob <- matrix( 1/N, nrow=N, ncol=d)
  }

  # check prob
  if(( nrow(prob) != N) | (ncol(prob) != d )) {
    stop( sprintf("probability matrix (prob) does not have correct dimensions\n Needed (nrow = %d, ncol = %d), provided (nrow = %d, ncol = %d)\n",
                 N, d, nrow(prob), ncol(prob) ) ) 
  }

  # get max prob
  maxProb <- apply(prob, 1, max) 

  totalProbability <- sum(maxProb)

  # create row major matrix for input to C program
  prob <- c(t(prob))


  #################### NEIGHBORS ######################################
 
  # get neighbors if they are missing 
  if( is.missing(neighbors) ) {
 
    # if the number of neighbors is missign obtain it here 
    if( is.missing(nNeighbors) ) nNeighbors = min( 10, N ) 

    # approximate the nearest neighbors
    neighbors <- knnx.index(x, x, k=nNeighbors)  
    
    # get neighbors and correct for C indexing
    neighbors <- c( 
      sapply( 
        1:nrow(neighborsIndex), 
        function(x) { 
          z <- neighborsIndex[x,] 
          z[ z != x][1:(length(z) - 1)] 
        } ) 
      ) - 1  # minus one here to adjust for 0 indexes in C
  } else {
    
    # if the number of neighbors is missign obtain it here 
    if( is.missing(nNeighbors) ) nNeighbors = ncol(neighbors) 

    # if neighbors are not missing make sure they are 'correct' 
    if(( nrow(neighbors) != N) | (ncol(neighbors) != nNeighbors )) {
      stop( sprintf("neighbors matrix (neighbors) does not have correct dimensions\n Needed (nrow = %d, ncol = %d), provided (nrow = %d, ncol = %d)\n",
                 N, nNeighbors, nrow(neighbors), ncol(neighbors) ) ) 
    }

    # check if neighbors are too big
    if( max(neighbors) > N ) stop( sprintf("max neighbor %d exceeds N = %d", max(neighbors), N) ) 

    # check if neighbors are too small
    if( min(neighbors) < 1 ) stop( sprintf("min neighbor %d is less than N = %d", max(neighbors), N) ) 

    # check neighbors for integers
    if( sum(!is.integer(neighbors) ) ) stop( "Some neighbors are not integers")

    # correct for C indexing, and transpose to make the matrix row major
    neighbors <- c(t(neighbors)) - 1
  }


  #################### TOTALS ######################################

  # get total 
  total <- colSums(x)
  
  # handle sample size 
  if( length(sampleSize) > 1 ) {
    if( H != length(sampleSize) ) {
      stop( "H != length(sampleSize)" ) 
    }
  } else {
    sampleSize <- rep(sampleSize/H,H )  
  }

  # group data together for input
  adminDbl <- c( PSUAcres, targetCV, targetVarWithin, total, sampleSize, cooling, tolSize, maxProb, prob, totalProb)
  adminDblLength <- length( adminDbl )
  adminInt <- c(segments, nNeighbors, neighbors) 
  adminIntLength <- length( adminInt )

  dup <- c() 

  print(total)
  print(k)
  print(d)
  print(N)
  print("AdminIntLength:")
  print(adminIntLength)
  print("AdminDoubleLength:")
  print(adminDblLength)

  print(cost)
 

  #################### RUN C FUNCTION ######################################
  
r.result <- .C("R_minCV",
  as.double(c(x)),      #checked
  as.integer(k),           #checked
  as.integer(d),           #checked 
  as.integer(N),           #checked
  as.integer(iterations),  #checked
  as.integer(label),       #checked
  as.double(cost),         #checked Q
  as.double(adminDbl),     #checked
  as.integer(adminInt),     #checked
  as.integer(adminIntLength), #checked
  as.integer(adminDblLength), #checked
  as.integer(dup),        #checked
  as.double(acceptRate)   #checked 
)

  print("C running time")
  print(proc.time() - Cprog) 
 

  #################### RETURN DATA ######################################
  
  a <- matrix(unlist(r.result[13]),ncol=3,byrow=T)
  colnames(a) <- c( 'change', 'U', 'accepted')

  myList <- list("accept"=a, "cost"=unlist(r.result[7]), "label"=unlist(r.result[6]))
    
  return(myList)
}
