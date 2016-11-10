library(saAlloc)

set.seed(400)
iter <- 1000

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
           2, # 13
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
targetCV <- c(.00,.00,.00)

# run minCV
b <- saMinCV(
  x,
  label,
  iter=iter,
  #cooling=0,
  targetCV=targetCV,
  sampleSize=8
)


# get accepted rows
accept <- b$accept[b$accept[,'accepted'] == 1,]
label.cv <- label

# init check data set
check <- c()
  

test.cv <- saAlloc:::.cv2( x, strata=label.cv, sampleSize=accept[1,names(b$sampleSize)] )
check <- rbind( check , c( accept[1,10:12], test.cv, accept[1,10:12] - test.cv) )

# calculate the CV's from the changes in the labels directly
if( iter > 1 ) {
  for( i in 2:nrow(accept) ) { 
    label.cv[ accept[i,'selected'] ] <- accept[i,'to'] 
    test.cv  <- saAlloc:::.cv2( x, strata=label.cv, sampleSize=accept[i,names(b$sampleSize)] )
  
    check <- rbind( check , c( accept[i,10:12], test.cv, accept[i,10:12] - test.cv) )
  }
}
print(check)



