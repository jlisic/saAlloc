library(saAlloc)
set.seed(45)

x <- matrix(rnorm( 100*4 ),ncol=4) + 10 

# create an initial set of strata assignments
label <- c(rep(0,50),rep(1,50))

x_total <- matrix(colSums(x),nrow=1)
x_var   <- as.matrix(aggregate( x,by=list(label), var)[,-1])
Nh_size <- aggregate( label,by=list(label), length)[,-1]
nh_size <- c(10,10)
                        
target_CV <- c(0.015, 0.015, 0.05, 0.08)

# run minCV for a sample size of 20
b <- saSampleAlloc(
                        x_total,
                        x_var,
                        targetCV=target_CV,
                        sampleSize=nh_size,
                        strataSize=Nh_size,
                        iterations=1,
                        penalty=100
                      )

nh_size_new <- b[[13]]

print(nh_size_new)

cv_init <- sqrt(colSums(x_var * Nh_size^2 * (1-nh_size/Nh_size) / nh_size )) / x_total 
print(cv_init)

cv_init_new <- sqrt(colSums(x_var * Nh_size^2 * (1-nh_size_new/Nh_size) / nh_size_new )) / x_total 
print(cv_init_new)

penalty_val <- ifelse( target_CV - cv_init_new < 0, 100 * (cv_init_new - target_CV),0)

print( sqrt(sum(cv_init_new^2)) + sum(penalty_val) )

