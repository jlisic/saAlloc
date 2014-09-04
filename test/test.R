# this is a test of SA
# currently this test just checks if things don't fail 
# todo: verify results
# todo: compare against verified results

# test the SA package
x <- matrix( c( 1, 3, 6,
                2, 3, 7,
                3, 4, 8,
                2, 5, 9
                ), byrow = T, nrow=4)

# create an initial set of strata assignments 
label <- c(0,0,1,1) 

cooling <- 20
seg <- c(5,6,7,9)
iterations <- 10

# useful for spatial sampling, set the size of the units
PSUArea <- c(3200, 3839, 4478, 5760)
targetVar <- c(10000, 20000, 30000, 40000) 



# run equal clus
print(" TESTING SA EQUAL CLUS (4 - Value Test)" )
a <- saEqualClus(
  x,
  label,
  iterations,
  cooling,
  seg,
  PSUArea,
  targetVar,
  rep(0,ncol(x)),
  0.05 
  )



# number of iterations
iterations <- 1000 
x <- matrix(rnorm( 100*4,mean=10 ),ncol=4)  
label <- c(rep(0,50),rep(1,50)) 
cooling <- 0
seg <- rep(5,100)
PSUAcres <- rep(1,100)
targetVar <- c(30,30,30,30)



# run equal clus
print(" TESTING SA EQUAL CLUS (Simulated Data Test)" )
a <- saEqualClus(
  x,
  label,
  iterations,
  cooling,
  seg,
  PSUArea,
  targetVar,
  rep(0,ncol(x)),
  0.05 
  )



# run minCV
print(" TESTING SA MIN CV(Simulated Data Test)" )
b <- saMinCV(
  x,
  label,
  targetCV=c( 0.08, 0.10, 1.00, 0.05),
  sampleSize=20
)






