# this is a test for averaging S^2 (e.g. k > 1)  

library(saAlloc)

set.seed(400)

# data set with 100 observations and 4 characteristics split between two strata

x <- matrix(1:16,ncol=4)


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

x3 <- 1:16

#a <- cbind(x1 + runif(16),x2 + runif(16), x3 + runif(16))
#b <- cbind(x1 + runif(16),x2 + runif(16), x3 + runif(16))

a <- cbind(x1 ,x2 , x3 )
b <- cbind(x1 ,x2 , x3 + label )


colnames(a) <- c('Turtle','Trout','Tuna')
colnames(b) <- c('Turtle','Trout','Tuna')


# interleave a and b
x <- t(matrix( rbind( t(a), t(b)), nrow=ncol(a)))
colnames(x) <- c('Turtle','Trout','Tuna')


# target variance
targetCV <- c(.01,.02,.03)
names(targetCV) <- c('Tuna','Turtle','Trout')



# run minCV
result <- saMinCV(
  x,
  label+1,
  iterations=1,
  cooling=0,
  targetCV=targetCV,
  sampleSize=8
)

#summary(result) 
#

truVar <- aggregate(x,         
by=list(        
rep(label+1,each=2)*100      
+ rep(1:2,length(label))      
),var)                   



