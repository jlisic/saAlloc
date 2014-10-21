
# summary for s3 saAlloc object
summary.saAlloc <- function(x) {

  cat("\nSummary of ") 
  cat(as.list(sys.call())[[2]])
  cat( ":\n\n")

  # print title
  cat("CVs:\n")
  # print data
  y <- cbind( x$CVStart, x$CV, x$targetCV ) 
  colnames(y) <-   c("Initial", "Final", "Target" ) 
  print( data.frame(y, stringsAsFactors=FALSE) ) 
 
  # print title
  cat("\nSample Size:\n") 
  # print data 
  y <- cbind( x$sampleSizeStart, x$sampleSize ) 
  colnames(y) <-   c("Initial", "Final") 
  print( data.frame(y, stringsAsFactors=FALSE) ) 


  cat("\nStrata Size:\n")
  y <- cbind( x$strataSizeStart, x$strataSize ) 
  colnames(y) <- c("Initial", "Final" ) 
  print( data.frame(y, stringsAsFactors=FALSE) ) 


  cat("\nAuxiliary Size Constraint:\n")
  y <- cbind( x$acresStart, x$acres ) 
  colnames(y) <- c("Initial", "Final" ) 
  print( data.frame(y, stringsAsFactors=FALSE) ) 

  cat("\nRun Time:\n")
  print(x$runTime)

}

