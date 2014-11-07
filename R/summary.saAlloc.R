
# summary for s3 saAlloc object
summary.saAlloc <- function(object,...) {

  cat("\nSummary of ") 
  cat(as.list(sys.call())[[2]])
  cat( ":\n\n")

  # print title
  cat("CVs:\n")
  # print data
  y <- cbind( object$CVStart, object$CV, object$targetCV ) 
  colnames(y) <-   c("Initial", "Final", "Target" ) 
  print( data.frame(y, stringsAsFactors=FALSE) ) 
 
  # print title
  cat("\nSample Size:\n") 
  # print data 
  y <- cbind( object$sampleSizeStart, object$sampleSize ) 
  colnames(y) <-   c("Initial", "Final") 
  print( data.frame(y, stringsAsFactors=FALSE) ) 


  cat("\nStrata Size:\n")
  y <- cbind( object$strataSizeStart, object$strataSize ) 
  colnames(y) <- c("Initial", "Final" ) 
  print( data.frame(y, stringsAsFactors=FALSE) ) 


  cat("\nAuxiliary Size Constraint:\n")
  y <- cbind( object$acresStart, object$acres ) 
  colnames(y) <- c("Initial", "Final" ) 
  print( data.frame(y, stringsAsFactors=FALSE) ) 

  cat("\nRun Time:\n")
  print(object$runTime)

}

