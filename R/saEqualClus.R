saEqualClus <-
function(
  x,
  label,
  iterations=1000,
  cooling=0,
  segments=rep(1,nrow(x)),
  controlVariable=rep(1,nrow(x)),
  targetVar,
  varWithin=rep(0,ncol(x)),
  tolSize=0.05
  ) {

  Cprog <- proc.time()
  
  unique.label <- unique(label)
  rlabel <- sapply(label, function(x){ which(x == unique.label) - 1} ) 

  n <- length(label)
  d <- ncol(x)
  k <- length(c(x)) / ( n * d) # number of repeat obs
  
  # get colnames for x
  if( is.null(colnames(x)) ) {
     x.colnames <- sprintf("%d",1:d)
  } else {
     x.colnames <- colnames(x)
  }
 
  cost <- rep(1,d) 

  adminDbl <- c( controlVariable, targetVar, varWithin, cooling, tolSize)
  adminDblLength <- length( adminDbl )
  adminInt <- c(segments) 
  adminIntLength <- length( adminInt )

  # CHECK FOR VIOLATION OF TOLSIZE
  checkAcres <- stats::aggregate(controlVariable,by=list(rlabel),sum)$x
  if( 1 - min(checkAcres) / max(checkAcres) > tolSize ) { 
    print("Initial acreages assumptions violated,")
    print("returning NULL")
    return(NULL)
  }
  
  dup <- c() 
  costChangeSize <- 6 + d 
  acceptRate <- rep(0,costChangeSize*iterations) 
  
  r.result <- .C("R_substrata2",
    as.double(c(x)),         # 1 checked
    as.integer(k),           # 2 checked
    as.integer(d),           # 3 checked 
    as.integer(n),           # 4 checked
    as.integer(iterations),  # 5 checked
    as.integer(rlabel),      # 6 checked
    as.double(cost),         # 7 checked Q
    as.double(adminDbl),     # 8 checked
    as.integer(adminInt),    # 9 checked
    as.integer(adminIntLength), # 10 checked
    as.integer(adminDblLength), # 11 checked
    as.integer(dup),            # 12 checked
    as.double(acceptRate),      # 13 checked 
    as.integer(iterations),     # 14 needed but not used
    as.integer(costChangeSize)
  )

  #print("C running time")
  runTime <- (proc.time() - Cprog)[3]
  
  newLabel <- sapply(unlist(r.result[6]), function(x) unique.label[x+1] ) 
                 
  a <- matrix(unlist(r.result[13]),ncol=costChangeSize,byrow=T) 
  colnames(a) <- c( 'change', 'U', 'T', 'selected', 'from', 'to',x.colnames)
  
  a.names <- c(colnames(a),'accepted')
  a <- cbind(a, as.numeric(a[,'U'] <= a[,'T']))
  colnames(a) <- a.names

  # calculate strata size
  strataSizeStart <- stats::aggregate(rlabel, by=list(rlabel), length)$x
  strataSize <- stats::aggregate(rlabel, by=list(newLabel), length)$x
  
  # calculate CV
  CVStart  <- .cv2( x, rlabel, strataSize, average=TRUE)
  CV  <- .cv2( x, newLabel, strataSizeStart, average=TRUE)
  CVTarget <- k*sqrt((targetVar + varWithin)*strataSizeStart)/ colSums(x)


  myList <- list(accept=a,
                 cost=r.result[[7]], 
                 label=newLabel,
                 sampleSize=NULL,
                 sampleSizeStart=NULL,
                 CV =CV,
                 CVStart=CVStart,
                 targetCV=CVTarget , 
                 strataSize=strataSize,
                 strataSizeStart=strataSizeStart,
                 runTime = runTime,
                 variables=x.colnames 
                 )
  class(myList) <- "saAlloc" 
    
  return(myList)
}
