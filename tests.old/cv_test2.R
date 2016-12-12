#source('~/src/saAlloc/R/saSampleAlloc.R')
#source('~/src/saAlloc/R/saMinCV.R')
#source('~/src/saAlloc/R/cv.R')
library(saAlloc)

set.seed(400)

Nh <- c(10,10,8)
K  <- 5 
R  <- 4

H <- length(Nh)

# create x
x <-c()
for(k in 1:K){
x <- cbind(x, c(unlist(apply( cbind(Nh,1:H), 1, function(x) rnorm(x[1]*R,mean=100,sd=x[2])))))
}

# create strata
strata <- c(unlist(apply( cbind(Nh,1:H),1, function(x) rep(x[2],x[1]) ))) -1
#bylist = list( rep(strata,each=2), rep( c(0,1),length(strata)) )

# create prob matrix
probMatrix <- matrix( 1/H, nrow=sum(Nh) ,ncol=H)


penalty  <- rep(10,K) 
sampleSize <- rep(4,H)
expectedVar <- sqrt(   sum( ( (1:H)^2 + 100) * (Nh^2/sampleSize)  )  ) / (100 * sum(Nh)) 
targetCV <- rep(expectedVar,K)


write.csv(x ,'test.csv')




if ( F ) {
test <- saSampleAlloc( 
  x,
  label=strata,
  targetCV=targetCV,
  sampleSize=sampleSize,
  weightMatrix=probMatrix,            # missing handled
  sampleSizeIterations=300,
  penalty=penalty,            # negative penalties are ignored
  preserveSatisfied=FALSE,
  locationAdjustment=x,
  scaleAdjustment=scaleAdjustment
  ) 

print( .cv2( x, strata=strata, sampleSize=sampleSize, average=TRUE, locationAdjustment=x, scaleAdjustment=scaleAdjustment) )
}


iterations <- 13
iterations <- 200


# test location adjustment
locationAdjustment <- sqrt(x)
scaleAdjustment <- sqrt(x)
p <- 2

test2 <- saMinCV( 
  x,
  label=strata,
  iterations=iterations,
  targetCV=targetCV,
  sampleSize=sampleSize,
  weightMatrix=probMatrix,            # missing handled
  sampleSizeIterations=10,
  penalty=penalty,            # negative penalties are ignored
  locationAdjustment=locationAdjustment,
  scaleAdjustment=scaleAdjustment,
  preserveSatisfied=TRUE
  ) 



if( F ) {
# compare saMinCV output to R output
print("initial")
print( saAlloc:::.cv2( x, strata=strata, sampleSize=sampleSize, average=TRUE, locationAdjustment=x))     #, locationAdjustment=x) )
      #, scaleAdjustment=x) 

strataNew <- strata
sampleSizeNew <- sampleSize
iterations <- iterations 

for( i in 1:(iterations+1) )  {

  if( i == 1) {  
    output <- c(1, saAlloc:::.cv2( x, strata=strata, sampleSize=sampleSize, average=TRUE, locationAdjustment=locationAdjustment, scaleAdjustment=scaleAdjustment) )
  } else {
    a <- test2$accept[i,]
    if( a['U'] <= a['T'] ) {
      strataNew[a['from']] = a['to']
      sampleSizeNew <- a[sprintf("n_%d",0:2)] 
      output <- rbind(output, c(i, saAlloc:::.cv2( x, strata=strataNew, sampleSize=sampleSizeNew, average=TRUE, locationAdjustment=locationAdjustment, scaleAdjustment=scaleAdjustment) ))
    } else {
      # make a copy of values so we don't change the accepted sequence
      strataNewTmp <- strataNew
      sampleSizeNewTmp <- sampleSizeNew
      strataNewTmp[a['from']] = a['to']
      sampleSizeNewTmp <- a[sprintf("n_%d",0:2)] 
      output <- rbind(output, c(i, saAlloc:::.cv2( x, strata=strataNewTmp, sampleSize=sampleSizeNewTmp, average=TRUE, locationAdjustment=locationAdjustment, scaleAdjustment=scaleAdjustment) ))
    }
  }
}

output.trial <- test2$accept[test2$accept[,'accepted']==1,sprintf("%d",1:5)]
output.trial <- test2$accept[,sprintf("%d",1:5)]

# get objective functions
a <- t( (t(output.trial) - targetCV) ) 
obj <- rowSums( output.trial^p + t( t(a)^p * penalty ) * (a > 0) )

obj.change <- obj[-1] - obj[cumsum(test2$accept[,'accepted'])][1:(length(obj)-1)]
obj.change <- c(-1, obj.change)


output.diff.change <- cbind( obj.change, test2$accept[,'change'])
print(max( output.diff.change[,1] - output.diff.change[,2]))

output.diff <- abs(output[,-1] - output.trial)

print(max(output.diff))

}



# distribution of accepted moves






#print("update")
#
#strata  = c(  1,   0,   0,   0,  0,  1,  1,   1,   1,   1, 2, 2, 2, 2, 2 );
#
#
#x.mean <- aggregate( x, by=list( c(strata,strata+10)), mean)
#x.var <- aggregate( x, by=list( c(strata,strata+10)), var)
#print(x.mean)
#print(x.var)

