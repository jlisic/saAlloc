library(saAlloc)

set.seed(400)

# data set with 100 observations and 4 characteristics split between two strata

x <- matrix(1:16,ncol=4)


label <- c(
           0, # 1
           0, # 2
           0, # 3
           0, # 4
           0, # 5
           1, # 6
           1, # 7
           1, # 8
           0, # 9
           0, # 10
           0, # 11
           1, # 12
           1, # 13
           2, # 14
           2, # 15
           2  # 15
           )

x1 <- c(
        1,  # 1
        1,  # 2
        1,  # 3
        1,  # 4
        1,  # 5
        1,  # 6
        1,  # 7
        1,  # 8
        2,  # 9
        2,  # 10
        2,  # 11
        2,  # 12
        2,  # 13
        2,  # 14
        2,  # 15
        2   # 16
        )
x2 <- c(
        3, # 1
        3, # 2
        3, # 3
        3, # 4
        3, # 5
        3, # 6
        3, # 7
        3, # 8
        4, # 9
        4, # 10
        4, # 11
        4, # 12
        4, # 13
        4, # 14
        4, # 15
        4  # 16
        )

x <- cbind(x1,x2,x2)

#target variance
targetCV <- c(.01,.01,.01)

# run minCV
b <- saMinCV(
  x,
  label,
  iter=20000,
  cooling=1/10,
  targetCV=targetCV,
  sampleSize=8
)


# get accepted rows
accept <- b$accept[b$accept[,'accepted'] == 1,]
label.cv <- label

# init check data set
check <- c()

# calculate the CV's from the changes in the labels directly
for( i in 1:nrow(accept) ) { 
  label.cv[ accept[i,'from'] ] <- accept[i,'to'] 
  test.cv  <- saAlloc:::.cv( x, strata=label.cv, sampleSize=accept[i,names(b$sampleSize)] )

  check <- rbind( check , c( accept[i,9:11], test.cv, accept[i,9:11] - test.cv) )
}

print(check)



