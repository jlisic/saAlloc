library(saAlloc)

set.seed(400)

# data set with 100 observations and 4 characteristics split between two strata

x <- matrix(1:16,ncol=4)

prob <- c(
           .95, .05,
           .95, .05,
           .95, .05,
           .95, .05,
           .95, .05,
           .95, .05,
           .95, .05,
           .95, .05,
           .05, .95,
           .05, .95,
           .05, .95,
           .05, .95,
           .05, .95,
           .05, .95,
           .05, .95,
           .05, .95
           )

prob <- matrix(prob,ncol=2,byrow=T)


label <- c(
           0,
           0,
           0,
           0, 
           0,
           1,
           1,
           1, 
           0,
           0,
           0,
           1,
           1,
           1,
           1,
           1
           )

x1 <- c(
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        2,
        2,
        2,
        2,
        2,
        2,
        2,
        2
        )
x2 <- c(
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        3,
        4,
        4,
        4,
        4,
        4,
        4,
        4,
        4
        )

x <- cbind(x1,x2,x2)

#target variance
targetCV <- c(.01,.01,.01)

# run minCV
b <- saMinCV(
x,
label,
iter=92,
cooling=0,
targetCV=targetCV,
sampleSize=8
, prob=prob
)

# compare results
print("before CV")
print( 2*sqrt(colSums( aggregate(x,by=list(label),var))[-1]) / colSums(x) )

print("Target")
print( sqrt(targetCV) )

print("Min CV")
print( 2*sqrt(colSums(aggregate(x,by=list(b$label),var))[-1]) / colSums(x) )


nHigh   <- 10
nMed    <- 5 
nLow <- 15

# high density cultivation, few farms
highCult <- cbind( 
                rnorm(nHigh, 500, 80),
                rlnorm(nHigh, 1,.3)
              )

# medium density cultivation
mediumCult <- cbind( 
                rnorm(nMed,200,100),
                rlnorm(nMed,1.5,.3)
              )

# low density cultivation, many farms and non-farms 
lowCult <- cbind( 
             rnorm( nLow,10,40),
             rlnorm( nLow,2)
           )
               

# handle possible negative values
highCult[highCult[,1] < 0,1] <- 0 
mediumCult[mediumCult[,1] < 0,1] <- 0 
lowCult[lowCult[,1] < 0,1] <- 0 

highCult[highCult[,2] > 5,2] <- 10 
mediumCult[mediumCult[,2] > 10,2] <- 10 
lowCult[lowCult[,2] > 40 ,2] <- 40 


x <- rbind( highCult, mediumCult, lowCult)

# targetCV
targetCV <- c(.05,1)

# label
kMean <- kmeans(scale(x),3)

label <- kMean$cluster

# find optimal sample size
library(Rsolnp)

# initial Allocation 
y <- c(5,5,5) 

# items Needed for allocation 
N     <- aggregate(label,by=list(label),length)[-1]
S2    <- aggregate(x,by=list(label),var)[-1]
Total <- aggregate(x,by=list(label),sum)[-1]
sampleSize <- sum(y)

# constraint function
constraintFunc <- function(y) return( (sum(y) - sampleSize)^2 )

# objective function
objFunc <- function(y){
  a <- colSums(S2 * N^2/y)/colSums(Total)^2  - targetCV^2
  v <- a * 1/(1+exp(-a))
  return( sum(v * v) + constraintFunc(y) )
}


# get optimal sample size
optSampleSize <- solnp( 
  y,
  fun=objFunc, 
#  eqfun=constraintFunc, 
  LB=rep(2,length(y)),
  UB=rep(50,length(y)) 
  )  


# adjust 
nOpt <- optSampleSize$pars 
print(nOpt)


# get prob
prob <- t(apply( scale(x), 1, function(x) dnorm(diag( t(t(kMean$centers) - x) %*% (t(kMean$centers) - x) )) ))

# make suer the narginal probability sums to 1
prob <- prob/rowSums(prob)

# run minCV
b <- saMinCV(
x,
label-1,
iter=30,
cooling=0,
targetCV=targetCV,
sampleSize=nOpt  
#,prob=prob
)

# compare results
print("before CV")
print( sqrt(colSums( aggregate(x,by=list(label),var)* N^2/nOpt)[-1]) / colSums(x) )

print("Target")
print( targetCV )

print("Min CV")
print( sqrt(colSums( aggregate(x,by=list(b$label),var)* N^2/nOpt)[-1]) / colSums(x) )


#par(mfrow=c(1,2))
#
#plot( x, col=label)
#plot( x, col=b$label+1)



