
# summary for s3 saAlloc object
summary.saAlloc <- function(object,...) {

  syscall <- as.list(sys.call())[[2]]

  y <- cbind( object$CVStart, object$CV, object$targetCV ) 
  colnames(y) <-   c("Initial", "Final", "Target" ) 
  CVs <- data.frame(y, stringsAsFactors=FALSE) 

  y <- cbind( object$sampleSizeStart, object$sampleSize ) 
  colnames(y) <-   c("Initial", "Final") 
  sampleSize <- data.frame(y, stringsAsFactors=FALSE)  

  y <- cbind( object$strataSizeStart, object$strataSize ) 
  colnames(y) <- c("Initial", "Final" ) 
  strataSize <-  data.frame(y, stringsAsFactors=FALSE)  


  if( !is.null(object$acres) ) { 
    y <- cbind( object$acresStart, object$acres ) 
    colnames(y) <- c("Initial", "Final" ) 
    acres <-  data.frame(y, stringsAsFactors=FALSE)  
  } else {
    aces <- NULL
  }

  runTime <- object$runTime

  ans <- list( syscall=syscall, CVs=CVs, sampleSize=sampleSize, strataSize=strataSize, runTime=runTime ) 
  class(ans) <- "summary.saAlloc"
  ans

}

# print for s3 summary object
print.summary.saAlloc <- function(x,...) {
  
  cat("\nSummary of ") 
  cat(x$syscall)
  cat( ":\n\n")

  # print title
  cat("CVs:\n")
  # print data
  print( x$CVs ) 

  # print title
  cat("\nSample Size:\n") 
  # print data 
  print( x$sampleSize )

  cat("\nStrata Size:\n")
  print( x$strataSize )


  if( !is.null(x$acres) ) { 
    cat("\nAuxiliary Size Constraint:\n")
    print( x$acres ) 
  }

  cat("\nRun Time:\n")
  print(x$runTime)

}


