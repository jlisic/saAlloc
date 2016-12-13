source('~/src/saAlloc/R/saSampleAlloc.R')
source('~/src/saAlloc/R/cv.R')


x <- matrix( c(

# strata 0
    1.0, 2.0, 3.0, 
    4.0, 4.0, 5.0, 

    5.0, 5.0, 1.0, 
    2.0, 3.0, 4.0, 

    2.0, 3.0, 4.0, 
    4.0, 5.0, 5.0, 

    5.0, 1.0, 2.0, 
    3.0, 4.0, 4.0, 

    5.0, 5.0, 5.0, 
    5.0, 5.0, 5.0, 

# strata 1
    0.0, 2.0, 3.0,  #a
    4.0, 5.0, 5.0,  #b
    
    5.0, 5.0, 0.0,  #a 
    2.0, 3.0, 4.0,  #b 
    
    2.0, 3.0, 4.0, 
    4.0, 6.0, 5.0, 
    
    5.0, 0.0, 2.0, 
    3.0, 4.0, 4.0, 
    
    3.0, 4.0, 4.0, #a 
    5.0, 7.0, 5.0, #b

# strata 2 
    0.0, 2.0, 3.0, 
    4.0, 5.0, 5.0, 
    
    5.0, 5.0, 0.0, 
    2.0, 3.0, 4.0, 
   
    2.0, 3.0, 4.0, 
    4.0, 6.0, 5.0, 
    
    5.0, 0.0, 2.0, 
    3.0, 4.0, 4.0, 
   
    3.0, 4.0, 4.0, 
    5.0, 7.0, 5.0 
), byrow=T, ncol=3)
  
probMatrix <- matrix( c( 
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,

    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,

    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2,
    0.5, 0.3, 0.2
  ), byrow=T,ncol=3) 

strata  = c(  0,   0,   0,   0,  0,  1,  1,   1,   1,   1, 2, 2, 2, 2, 2 );

bylist = list( rep(strata,each=2), rep( c(0,1),length(strata)) )

x.mean <- aggregate( x, by=bylist, mean)
x.var <- aggregate( x, by=bylist, var)
print(x.mean)
print(x.var)

cv <- 
sqrt(aggregate( x.var[,3:5]/3, by=list(x.var$Group.2), sum)) /
aggregate( x.mean[,3:5], by=list(x.var$Group.2), sum) 

cv <- colMeans(cv)

targetCV <- c( 1.0, 1.1, 1.2 )
penalty  <- c(10,10,10)
sampleSize <- c(3,3,3)





saSampleAlloc( 
  x,
  label=strata,
  targetCV=targetCV,
  sampleSize=sampleSize,
  weightMatrix=probMatrix,            # missing handled
  sampleSizeIterations=3,
  penalty=penalty,            # negative penalties are ignored
  preserveSatisfied=FALSE
  ) 

print(cv)
print( .cv2( x, strata=strata, sampleSize=sampleSize, average=TRUE) )






#print("update")
#
#strata  = c(  1,   0,   0,   0,  0,  1,  1,   1,   1,   1, 2, 2, 2, 2, 2 );
#
#
#x.mean <- aggregate( x, by=list( c(strata,strata+10)), mean)
#x.var <- aggregate( x, by=list( c(strata,strata+10)), var)
#print(x.mean)
#print(x.var)

