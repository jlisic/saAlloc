
# plot for s3 saAlloc object
# together == T, all vars on one plot
# together == F, separate plots for each commodity
plot.saAlloc <- function( x, together=TRUE ) {
  
  # get variables 
  accept <- x$accept
  variables <- x$variables

  ### pars for multiplots
  if(!together ) {
    plotRow <-  ceiling(sqrt(length(variables)))
    plotCol <-  ceiling(length(variables) / plotRow)
    par(mfrow= c(plotRow,plotCol) )
  }


  #### plot all variables on the same plot

  # plot first variable 
  plot(
    1:iterations,
    accept[,variables[1]],
    type='l', 
    ylim=c(0,max(accept[,variables])) ,
    xlab='Iterations',
    ylab='CV'
  )

  # add title
  if(together) {
    title("Trace of CVs")
  } else {
    title(sprintf("Trace of CV for %s",variables[1]))
  }

  # plot other variables 
  if(length(variables) > 1 ) {

    if( together ) {
      for( i in 2:length(variables) ) {
        points(
            1:iterations,
            accept[,variables[i] ],
            type='l',
            col=i
          )
      }
    } else {
      for( i in 2:length(variables) ) {
        plot(
          1:iterations,
          accept[,variables[i] ],
          type='l',
          col=i,
          xlab='Iterations',
          ylab='CV'
        )
        title(sprintf("Trace of CV for %s",variables[i]))
      }
    }


  }

  # add legend
  if( together ) {
    legend( 
      "topright",
      variables,                     # add names
      lty=rep(1,length(variables)),  # add lines
      col = 1:length(variables)      # add color
    )
  } 
}


