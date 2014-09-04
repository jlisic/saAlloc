saMinCV <-
function(
  x,
  label,
  iterations=1000,
  cooling=0,
  segments=rep(1,nrow(x)),
  PSUAcres=rep(1,nrow(x)),
  targetCV,
  sampleSize,
  targetVarWithin=rep(0,ncol(x)),
  tolSize=1
  ) {

  
  Cprog <- proc.time()
  
  rlabel <- label
  N <- length(label)
  d <- ncol(x)
  k <- length(c(x)) / ( N * d)

  H <- length(unique(label))
 
  acceptRate <- rep(0,3*iterations) 
  cost <- rep(1,d+1) 
  
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
  adminDbl <- c( PSUAcres, targetCV, targetVarWithin, total, sampleSize, cooling, tolSize)
  adminDblLength <- length( adminDbl )
  adminInt <- c(segments) 
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
  
  a <- matrix(unlist(r.result[13]),ncol=3,byrow=T)
  colnames(a) <- c( 'change', 'U', 'accepted')

  myList <- list("accept"=a, "cost"=unlist(r.result[7]), "label"=unlist(r.result[6]))
    
  return(myList)
}
