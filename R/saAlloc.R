
saAlloc <- function(
   y, 
   label, 
   targetCV, 
   targetVar,
   sampleSize, 
   weightMatrix, 
   method='saMinCV',
   iterations=1000, 
   sampleIterations=100, 
   sampleUpdateIterations=100,
   cooling=0, 
   segments=rep(1,length(label)), 
   PSUAcres=rep(1,length(label)), 
   targetVarWithin=rep(0,ncol(y)),
   preserveSatisfied=TRUE,
   tolSize=1
) {
    
  if(method=='saMinCV') {
    x <- saMinCV(
      y, 
      label, 
      targetCV, 
      sampleSize, 
      weightMatrix, 
      iterations, 
      sampleIterations, 
      sampleUpdateIterations,
      cooling, 
      segments, 
      PSUAcres, 
      targetVarWithin,
      preserveSatisfied
    )

    return(x)
  }
 
   

  if(method=='saMinCV') { 
    x <- saEqualClus(
      y,
      label,
      iterations,
      cooling,
      segments,
      PSUAcres,
      targetVar,
      targetVarWithin,
      tolSize
    )
    
    return(x)
  }

  stop("Invalid Method")

}

