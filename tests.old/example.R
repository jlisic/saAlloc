library(saAlloc)

x <- matrix(rnorm( 100*4 ),ncol=4)

# create an initial set of strata assignments
label <- c(rep(0,50),rep(1,50))

# run minCV for a sample size of 20
b <- saSampleAlloc(
                        x,
                        label,
                        targetCV=c(0.08, 0.10, 0.05, 0.08),
                        sampleSize=20
                      )
