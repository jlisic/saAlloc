saEqualClus <-
function(
  x,
  label,
  iterations=1000,
  cooling=0,
  segments=rep(1,nrow(x)),
  PSUAcres=rep(1,nrow(x)),
  targetVar,
  targetVarWithin=rep(0,ncol(x)),
  tolSize=0.05
  ) {

  
  Cprog <- proc.time()
  
  rlabel <- label
  n <- length(label)
  d <- ncol(x)
  k <- length(c(x)) / ( n * d)
 
  acceptRate <- rep(0,3*iterations) 
  cost <- rep(1,d) 

  adminDbl <- c( PSUAcres, targetVar, targetVarWithin, cooling, tolSize)
  adminDblLength <- length( adminDbl )
  adminInt <- c(segments) 
  adminIntLength <- length( adminInt )


  # CHECK FOR VIOLATION OF TOLSIZE
  checkAcres <- aggregate(PSUAcres,by=list(label),sum)$x
  if( 1 - min(checkAcres) / max(checkAcres) > tolSize ) { 
    print("Initial Acreages Assumptions violated")
    print("Returning NULL")
    return(NULL)
  }
  
  dup <- c() 

  print(k)
  print(d)
  print(n)
  print(adminIntLength)
  print(adminDblLength)

  print(cost)
  
r.result <- .C("R_substrata2",
  as.double(c(x)),      #checked
  as.integer(k),           #checked
  as.integer(d),           #checked 
  as.integer(n),           #checked
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

  myList <- list("accept"=a,
                 "cost"=r.result[7], 
                 "label"=r.result[6])
    
  return(myList)
}
