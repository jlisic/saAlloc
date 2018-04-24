## test for anticipated variance


set.seed(100)

N <- 1000
Beta <- 1 
Sigma <- 10
Gam   <- 0.75
fpc <- TRUE
iter <- 1000


x <- runif(N, min=0, max=300)

locAdj <- x^Gam * Sigma 


run_it <- FALSE
result <- c()

for( i in 1:10) {


y <-  x * Beta + rnorm(N,sd=locAdj)

#target variance
targetCV <- c(0.1)

label <- rep(1:3,length.out=N)

library(saAlloc)

# run minCV with adjustmetn
strata_without_adjustment <- saMinCV(
  matrix(x,ncol=1),
  label,
  iter=iter,
  targetCV=targetCV,
  sampleSize=8,
  fpc=fpc
)

# run minCV with adjustmetn
strata_with_adjustment <- saMinCV(
  matrix(x,ncol=1),
  label,
  iter=iter,
  targetCV=targetCV,
  locationAdjustment= locAdj^2,
  sampleSize=8,
  fpc=fpc
)

table_with_adjustment <- data.frame(
  Nh=strata_with_adjustment$strataSize,
  nh=strata_with_adjustment$sampleSize,
  var=aggregate(y,by=list(strata_with_adjustment$label),var)$x
)

table_without_adjustment <- data.frame(
  Nh=strata_without_adjustment$strataSize,
  nh=strata_without_adjustment$sampleSize,
  var=aggregate(y,by=list(strata_without_adjustment$label),var)$x
)

check_adjustment <- data.frame(
  Nh=strata_with_adjustment$strataSize,
  nh=strata_with_adjustment$sampleSize,
  var=aggregate(locAdj^2,by=list(strata_with_adjustment$label),mean)$x
)

check_x <- data.frame(
  Nh=strata_with_adjustment$strataSize,
  nh=strata_with_adjustment$sampleSize,
  var=aggregate(x,by=list(strata_with_adjustment$label),var)$x
)

colnames(check_adjustment) <- c('Nh','nh','adj')
colnames(check_x) <- c('Nh','nh','var')

colnames(table_with_adjustment) <- c('Nh','nh','var')
colnames(table_without_adjustment) <- c('Nh','nh','var')


result <- rbind(result,
c(  
      sqrt( sum( table_with_adjustment$Nh^2/table_with_adjustment$nh * 
                (1-table_with_adjustment$nh/table_with_adjustment$Nh) * table_with_adjustment$var)) / sum(y),
      sqrt( sum( table_without_adjustment$Nh^2/table_without_adjustment$nh * 
                (1-table_without_adjustment$nh/table_without_adjustment$Nh) * table_without_adjustment$var)) / sum(y) 
)
)
      

#check_cv <- sqrt( sum( 
#                      check_adjustment$Nh^2/check_adjustment$nh * (1-check_adjustment$nh/check_adjustment$Nh) * (check_x$var + check_adjustment$adj)
#                      )) / sum(x)

}

#print( sprintf("CV of y: with adjustment %f without %f", 
#      sqrt( sum( table_with_adjustment$Nh^2/table_with_adjustment$nh * 
#                (1-table_with_adjustment$nh/table_with_adjustment$Nh) * table_with_adjustment$var)) / sum(y),
#      sqrt( sum( table_without_adjustment$Nh^2/table_without_adjustment$nh * 
#                (1-table_without_adjustment$nh/table_without_adjustment$Nh) * table_without_adjustment$var)) / sum(y) 
#      )
#)


set.seed(100)

av_check <- c()
for( i in 1:100000) {
y <-  x * Beta + rnorm(N,sd=locAdj)

label <- rep(1:3,length.out=N)

table_av_check <- data.frame(
  Nh=aggregate(label,by=list(label),length)$x,
  nh=rep(20,3),
  var=aggregate(y,by=list(label),var)$x
)

colnames(table_av_check) <- c('Nh','nh','var')
      
av_check <- c( av_check, sum( table_av_check$Nh^2/table_av_check$nh * 
                (1-table_av_check$nh/table_av_check$Nh) * table_av_check$var))   

}

print(mean(av_check))
print(sd(av_check)/sqrt(100000))




table_av_true <- data.frame(
  Nh=aggregate(label,by=list(label),length)$x,
  nh=rep(20,3),
  var=aggregate(x,by=list(label),var)$x
)
check_adjustment <- data.frame(
  Nh=aggregate(label,by=list(label),length)$x,
  nh=rep(20,3),
  var=aggregate(locAdj^2,by=list(label),mean)$x
)
colnames(check_adjustment) <- c('Nh','nh','adj')
colnames(table_av_true) <- c('Nh','nh','var')
      
av_true <- sum( table_av_true$Nh^2/table_av_true$nh * 
                (1-table_av_true$nh/table_av_true$Nh) * (check_adjustment$adj + table_av_true$var))

print(av_true)
print(
(1.96 * sd(av_check) + mean(av_check) > av_true ) & 
(-1.96 * sd(av_check) + mean(av_check) < av_true)
)




