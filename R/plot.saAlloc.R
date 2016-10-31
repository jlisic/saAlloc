

# plot for s3 saAlloc object
# together == T, all vars on one plot
# together == F, separate plots for each commodity
plot.saAlloc <- function( x, together=TRUE, ... ) {
  
  # get variables 
  accept <- x$accept
  variables <- x$variables
  iterations <- nrow(accept)

  ### pars for multiplots
  if(!together ) {
    plotRow <-  ceiling(sqrt(length(variables)))
    plotCol <-  ceiling(length(variables) / plotRow)
    graphics::par(mfrow= c(plotRow,plotCol) )
  }


  #### plot all variables on the same plot
  
  accepted <- accept[,'accepted']
  
  
  if( sum(accepted) == 0 ) {
    print("No change to plot")
    return
  }

  a <- accept[,'accepted']
  changed <- a

for( i in 1:iterations ) {
  if( a[i] == 1 ) lastAccept = i
  changed[i] = lastAccept
}


# print("getting changes")
#  # get changes 
#  for( i in 1:iterations )  {
#    if( i == 1) {  
#      changed <- 1
#      lastAccept <- changed 
#    } else {
#      a <- accept[i,]
#     if( a['U'] <= a['T'] ) {
#        lastAccept <- i
#        changed <- c(changed, lastAccept) 
#      } else {
#        # make a copy of values so we don't change the accepted sequence
#        changed <- c(changed, lastAccept) 
#      }
#    }
#  }
# print("done getting changes")

  # plot first variable 
  graphics::plot(
    0:(iterations-1),
    accept[changed,variables[1]],
    type='l', 
    ylim=range(accept[changed,variables]) ,
    xlab='Iterations',
    ylab='CV'
  )
  

  # add title
  if(together) {
    graphics::title("Trace of CVs")
  } else {
    graphics::title(sprintf("Trace of CV for %s",variables[1]))
  }

  # plot other variables 
  if(length(variables) > 1 ) {

    if( together ) {
      for( i in 2:length(variables) ) {
        graphics::points(
            0:(iterations-1),
            accept[changed ,variables[i] ],
            type='l',
            col=i
          )
      }
    } else {
      for( i in 2:length(variables) ) {
        graphics::plot(
          0:(iterations-1),
          accept[changed,variables[i] ],
          type='l',
          col=i,
          xlab='Iterations',
          ylab='CV'
        )
        graphics::title(sprintf("Trace of CV for %s",variables[i]))
      }
    }
  }

  # add legend
  if( together ) {
    graphics::legend( 
      "topright",
      variables,                     # add names
      lty=rep(1,length(variables)),  # add lines
      col = 1:length(variables)      # add color
    )
  } 
}


