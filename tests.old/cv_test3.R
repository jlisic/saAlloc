#source('~/src/saAlloc/R/saSampleAlloc.R')
#source('~/src/saAlloc/R/saMinCV.R')
#source('~/src/saAlloc/R/cv.R')
library(saAlloc)
library(MASS)
library(ggplot2)


plotDir <- "~/src/saAllocDoc/plots"

set.seed(400)

Nh <- c(1,1,8)*10
K  <- 2 
R  <- 1 

H <- length(Nh)

# create x
x <-c()
for(k in 1:K){
x <- cbind(x, c(unlist(apply( cbind(Nh,1:H), 1, function(x) rnorm(x[1]*R,mean=100 + x[2]^2 ,sd=x[2])))))
}

# create strata
strata <- c(unlist(apply( cbind(Nh,1:H) ,1, function(x) rep(x[2],x[1]) ))) -1
                      
y <- t( matrix( t(x) , nrow=K*R) ) 

init.kmeans <- kmeans(y, H ) 
strata <- init.kmeans$cluster - 1

# create prob matrix
probMatrix <- matrix( 1/H, nrow=sum(Nh) ,ncol=H)


penalty  <- rep(100,K) 
sampleSize <- rep(4,H)
#expectedVar <- sqrt(   sum( ( (1:H)^2 + 100) * (Nh^2/sampleSize)  )  ) / (100 * sum(Nh)) 
#targetCV <- rep(expectedVar,K)
targetCV <- c( 0.0040, 0.0050 )


print( saAlloc:::.cv2( x, strata=strata, sampleSize=sampleSize, average=TRUE) )


if ( F ) {
test <- saSampleAlloc( 
  x,
  label=strata,
  targetCV=targetCV,
  sampleSize=sampleSize,
  sampleSizeIterations=10000,
  penalty=penalty            # negative penalties are ignored
  ) 

print( saAlloc:::.cv2( x, strata=strata, sampleSize=sampleSize, average=TRUE))
}



#iterations <- 41024
iterations <- 1000
iterations <- 225

# test location adjustment
#locationAdjustment <- sqrt(x)
#scaleAdjustment <- sqrt(x)
p <- 2





colnames(x) <- c('x','y')

set.seed(100)
test2 <- saMinCV( 
  x,
  label=strata,
  iterations=iterations,
  targetCV=targetCV,
  sampleSize=sampleSize,
  weightMatrix=probMatrix,            # missing handled
  sampleSizeIterations=3,
  penalty=penalty,            # negative penalties are ignored
  cooling=0.0,
  #locationAdjustment=locationAdjustment,
  #scaleAdjustment=scaleAdjustment,
  preserveSatisfied=TRUE
  ) 

print(summary(test2))
stop("hi2u")

x.var <- aggregate( x, by=list(strata), var)[,-1]
x.N <- aggregate( x, by=list(strata), length)[,-1]
print( sqrt(colSums( x.var * x.N^2 / sampleSize)) / colSums(x) )






if( T ) {
  # compare saMinCV output to R output
  print("initial")
  print( saAlloc:::.cv2( x, strata=strata, sampleSize=sampleSize, average=TRUE ))
  
  strataNew <- strata
  sampleSizeNew <- sampleSize
  iterations <- iterations 
  acceptSequence <- c() 
  lastAccept <- 1 
  
  for( i in 1:(iterations+1) )  {
 
    if( i %% 1000 == 0) print(i)

    if( i == 1) {  
      output <- c(1, saAlloc:::.cv2( x, strata=strata, sampleSize=sampleSize, average=TRUE))
    } else {
      a <- test2$accept[i,]
     if( a['U'] <= a['T'] ) {
        strataNew[a['selected']] = a['to']
        sampleSizeNew <- a[sprintf("n_%d",0:(H-1))] 
        output <- rbind(output, c(i, saAlloc:::.cv2( x, strata=strataNew, sampleSize=sampleSizeNew, average=TRUE)))
        lastAccept <- i
        acceptSequence <- c(acceptSequence, lastAccept) 
      } else {
        # make a copy of values so we don't change the accepted sequence
        strataNewTmp <- strataNew
        sampleSizeNewTmp <- sampleSizeNew
       strataNewTmp[a['selected']] = a['to']
        sampleSizeNewTmp <- a[sprintf("n_%d",0:(H-1))] 
        output <- rbind(output, c(i, saAlloc:::.cv2( x, strata=strataNewTmp, sampleSize=sampleSizeNewTmp, average=TRUE)))
        acceptSequence <- c(acceptSequence, lastAccept) 
      }
    }
  }

  output.trial <- test2$accept[,c('x','y')]

  # get objective functions
  a <- t( (t(output.trial) - targetCV) ) 
  obj <- rowSums( output.trial^p + t( t(a)^p * penalty ) * (a >= 0) )
  
  obj.change <- obj[-1] - obj[acceptSequence]
  obj.change <- c(-1, obj.change)


  output.diff.change <- abs(cbind( obj.change, test2$accept[,'change']))
  print(max( output.diff.change[,1] - output.diff.change[,2]))

  output.diff <- abs(output[,-1] - output.trial)
  print(max(output.diff))

}



# distribution of accepted moves

moveFreq <- function( y ) {
  # get transitions
  y <- y$accept[-1,]

  # get accepted transitions
  y <- y[ y[,'accepted'] == 1, ] 

  # create a table of accepted transitions
  a <- aggregate( y[,'accepted'],
           by=list(
             y[,'selected'],
             y[,'from'],
             y[,'to']
             ), sum)

  colnames(a) <- c('selected','from','to','count')
  return(a)
}


a <- moveFreq(test2)
b <- a[ a$from != a$to,]


if ( T ) {


# initial
df <- data.frame( (y[,1:2]) )
colnames(df) <- c('x','y')
df$Stratum <- as.factor(strata)

png(filename=sprintf("%s/init.png",plotDir))
print( ggplot(df, aes(x=x,y=y,color=Stratum),size=1.5) + geom_point())
dev.off()

# final
df <- data.frame( (y[,1:2]) )
colnames(df) <- c('x','y')
df$Stratum <- as.factor(test2$label)

png(filename=sprintf("%s/final.png",plotDir))
print( ggplot(df, aes(x=x,y=y,color=Stratum),size=1.5) + geom_point())
dev.off()







############## frequency of exchange plots ############################
freq <- aggregate( b$count, by=list(b$selected), sum)
freq.all <- rep(0,length(strata))

freq.all[freq$Group.1] <- freq$x

freq.density <- c(unlist(apply( cbind(freq.all + 1,1:length(strata) ) ,1, function(x) rep(x[2],x[1]) ))) 

df <- y[freq.density,1:2]
df <- data.frame(df)
colnames(df) <- c('x','y')

df.y <- cbind( y[,1:2], test2$label+1)
df.y <- data.frame(df.y)
colnames(df.y) <- c('x','y','c')
df.y$c <- as.factor(df.y$c)

png(filename=sprintf("%s/density.png",plotDir))
print(
ggplot(df, aes(x=x,y=y)) + 
  stat_density2d(aes(fill=..level..), geom='polygon') + scale_fill_gradient("Frequency",low="grey",high='white')  
  + geom_point(data=df.y, aes(x=x,y=y,color=c),size=1.5) + guides(color=FALSE)
)
dev.off()


png(filename=sprintf("%s/hist.png",plotDir))
hist(freq.all, main="Histogram of Exchanges per Sampling Unit", xlab="Exchanges", ylab="Sampling Units")
dev.off()


png(filename=sprintf("%s/trace.png",plotDir))
plot( test2 )
dev.off()
}


#print("update")
#
#strata  = c(  1,   0,   0,   0,  0,  1,  1,   1,   1,   1, 2, 2, 2, 2, 2 );
#
#
#x.mean <- aggregate( x, by=list( c(strata,strata+10)), mean)
#x.var <- aggregate( x, by=list( c(strata,strata+10)), var)
#print(x.mean)
#print(x.var)

