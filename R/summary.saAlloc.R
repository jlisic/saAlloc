
# summary for s3 saAlloc object
summary.saAlloc <- function(x) {

  # function to cat matrix by row
  catMatrix <- function(x) {
    apply(x,1,function(z) {
          for( i in 1:length(z) ) {

            # determine how to print
            if(is.integer(z)) {
              cat( sprintf("%d\t",z[i]) ) 
            } else if(is.double(z) ) {
              cat( sprintf("%.5f\t",z[i]) ) 
            } else cat( sprintf("hi%s\t",z[i]) ) 

          }
          cat('\n')
    })
  }

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

