
# summary for s3 saAlloc object
summary.saAlloc <- function(object,...) {

  syscall <- as.list(sys.call())[[2]]

  y <- cbind( object$CVStart, object$CV, object$targetCV ) 
  colnames(y) <-   c("Initial", "Final", "Target" ) 
  CVs <- data.frame(y, stringsAsFactors=FALSE) 

  if( !is.null(object$sampleSizeStart) ) {
    y <- cbind( object$sampleSizeStart, object$sampleSize ) 
    colnames(y) <-   c("Initial", "Final") 
    sampleSize <- data.frame(y, stringsAsFactors=FALSE)  
  } else {
    sampleSize <- NULL
  }

  y <- cbind( object$strataSizeStart, object$strataSize ) 
  colnames(y) <- c("Initial", "Final" ) 
  strataSize <-  data.frame(y, stringsAsFactors=FALSE)  


  if( !is.null(object$control) ) { 
    y <- cbind( object$controlStart, object$control ) 
    colnames(y) <- c("Initial", "Final" ) 
    control <-  data.frame(y, stringsAsFactors=FALSE)  
  } else {
    control <- NULL
  }

  runTime <- object$runTime

  ans <- list( 
              syscall=syscall, 
              CVs=CVs, 
              sampleSize=sampleSize, 
              strataSize=strataSize, 
              runTime=runTime,
              control=control ) 
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

  if( !is.null(x$sampleSize) ) { 
    cat("\nSample Size:\n") 
    print( x$sampleSize )
  }

  cat("\nStrata Size:\n")
  print( x$strataSize )

  cat("\nRun Time:\n")
  print(x$runTime)

}


