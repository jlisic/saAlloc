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

x <- cbind(x1,x2,x2)

colnames(x) <- c('Turtle','Trout','Tuna')


#target variance
targetCV <- c(.01,.02,.03)
names(targetCV) <- c('Tuna','Turtle','Trout')



# run minCV
b <- saMinCV(
  x,
  label,
  iterations=20,
  cooling=0,
  targetCV=targetCV,
  sampleSize=8
)

summary(b) 

sampleSize <- c(6,2)
names(sampleSize) <- c(1,0)


# run minCV
b <- saMinCV(
  x,
  label,
  iterations=20,
  cooling=10,
  targetCV=targetCV,
  sampleSize=sampleSize
)

summary(b) 

