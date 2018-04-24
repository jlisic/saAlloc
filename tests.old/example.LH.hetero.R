library(saAlloc)
library(SamplingStrata)
library(stratification)

set.seed(200)

pop_size <- 5000


beta1 <- c( 1, 0)   # correlated with x1
beta5 <- c( 0, 1)   # not correlated with x1

beta2 <- c(3/4,1/4)  # neg, cor. with x1 and x2 to some extent (iterated over)
beta3 <- c(1/2,1/2)  # neg, cor. with x1 and x2 to some extent (iterated over)
beta4 <- c(1/4,3/4)  # neg, cor. with x1 and x2 to some extent (iterated over)



v <- rep(1,5) # sd

# parameters
P <- 5 # number of Y's 
R <- 2 # number of X's
iterations_sa <- 250000
iterations_ga <- 5
cool = 1/10 
h <- 6 
cvs <- rep(.04,5)
names(cvs) <- sprintf("CV%d",1:P)
gamm <- 3/4

# fixed quantity
z <- cbind(matrix(rchisq(pop_size*R,3),ncol=R))*50


# deviates and model
x <- matrix(0,nrow=NROW(z),ncol=P) 
x[,1] <- z %*% beta1
x[,2] <- z %*% beta2
x[,3] <- z %*% beta3
x[,4] <- z %*% beta4
x[,5] <- z %*% beta5



# heteroscedastic model
x_adj <- matrix(0,nrow=NROW(x),ncol=P) 
x_adj[,1] <- v[1] * z[,1]^gamm
x_adj[,2] <- v[2] * sqrt(rowSums(z^2))^gamm 
x_adj[,3] <- v[3] * sqrt(rowSums(z^2))^gamm 
x_adj[,4] <- v[4] * sqrt(rowSums(z^2))^gamm 
x_adj[,5] <- v[5] * z[,2]^gamm 


##############################################################
# UNIVARIATE
##############################################################

# LH
p <- 1
LHkozak <- strata.LH(x = x[,p], CV=cvs[p], Ls = h, alloc = c(0.5, 0, 0.5), takeall = 0, algo = "Kozak",
           model='linear', 
           model.control=list(
             beta=1,
             sig2=1,
             gamma=1.5  
           )
           )

# simulated annealing (ACV)
penalty1 <-c(100,0,0,0,0) 
sa_solution_acv <- saMinCV(
  x,
  label=LHkozak$stratumID,
  targetCV=cvs,
  sampleSize=LHkozak$nh,
  iterations=iterations_sa,
  locationAdjustment=x_adj^2,
  cooling=cool,
  preserveSatisfied=FALSE,
  penalty=penalty1
)

# simulated annealing (CV)
sa_solution_cv <- saMinCV(
  x,
  label=LHkozak$stratumID,
  targetCV=cvs,
  sampleSize=LHkozak$nh,
  sampleSizeIterations=5,
  iterations=iterations_sa,
  cooling=cool,
  preserveSatisfied=FALSE,
  penalty=penalty1
)

############ Collect univariate results ############

result_lh_acv_uni <- saAlloc:::.cv2( x=x, strata=LHkozak$stratumID, sampleSize=LHkozak$nh, locationAdjustment=x_adj^2 )

result_sa_acv_uni <- saAlloc:::.cv2( x=x, strata=sa_solution_acv$label, sampleSize=sa_solution_acv$sampleSize, locationAdjustment=x_adj^2 )

result_sa_cv_uni <- saAlloc:::.cv2( x=x, strata=sa_solution_cv$label, sampleSize=sa_solution_cv$sampleSize, locationAdjustment=x_adj^2 )

result_uni <- rbind( result_lh_acv_uni,result_sa_cv_uni, result_sa_acv_uni)

rownames(result_uni)  <- c("LH (ACV)","SA CV-Optimized (ACV)","SA ACV-Optimized (ACV)")
print(result_uni)
      
      
      
##############################################################
# MULTIVARIATE 
##############################################################
h <- 5 


######## K MEANS STRATIFICATOIN ############

## k-means clutering
kMeansCluster <- kmeans(x, h)$cluster 

colnames(x) <- sprintf("X%d",1:5)

# create some helpful auxiliary variables
# to speed things up
bins <- 20
X <- x
X[,1] <- var.bin( x[,1], bins= bins) 
X[,2] <- var.bin( x[,2], bins= bins) 
X[,3] <- var.bin( x[,3], bins= bins) 
X[,4] <- var.bin( x[,4], bins= bins) 
X[,5] <- var.bin( x[,5], bins= bins) 

frame <- data.frame(
    Y1 = x[,1],
    Y2 = x[,2],
    Y3 = x[,3],
    Y4 = x[,4],
    Y5 = x[,5],
    M1 = x[,1],
    M2 = x[,2],
    M3 = x[,3],
    M4 = x[,4],
    M5 = x[,5],
#    X1 = X[,1], 
#    X2 = X[,2], 
#    X3 = X[,3], 
#    X4 = X[,4], 
#    X5 = X[,5], 
    X1 = 1:NROW(x),
    S1 = rep(0, NROW(x)),
    S2 = rep(0, NROW(x)),
    S3 = rep(0, NROW(x)),
    S4 = rep(0, NROW(x)),
    S5 = rep(0, NROW(x)),
    cens = rep(0, NROW(x)),
    cost = rep(1, NROW(x)),
    N = rep(1, NROW(x)),
    domainvalue = rep(1, NROW(x)))

rownames(frame) <- 1:NROW(x)

strata <- buildStrataDF(frame)

cv <- data.frame(DOM = "DOM1", t(cvs), domainvalue = 1)



# multivariate

ga_solution <- optimizeStrata(
  cv, 
  strata, 
  cens = NULL, 
  strcens = FALSE,
  alldomains = TRUE, 
  dom = NULL, 
  initialStrata = 5, 
  addStrataFactor = 0.0,
  minnumstr = 2, 
  iter = iterations_ga, 
  mut_chance = 0.0005,
  elitism_rate = 0.2, 
  highvalue = 1e08, 
  suggestions = NULL,
  realAllocation = TRUE, 
  writeFiles = FALSE)


# create sample size
sampleSizeMulti=sum(ceiling(ga_solution$aggr_strata$SOLUZ))

penalty2 <-c(100,100,100,100,100)

# simulated annealing ACV
sa_solution_acv_multi <- saMinCV(
  x,
  label=kMeansCluster,
  targetCV=cvs,
  sampleSize=sampleSizeMulti,
  iterations=iterations_sa,
  locationAdjustment=x_adj^2,
  cooling=cool,
  preserveSatisfied=FALSE,
  penalty=penalty2
)

# simulated annealing ACV
sa_solution_cv_multi <- saMinCV(
  x,
  label=kMeansCluster,
  targetCV=cvs,
  sampleSize=sampleSizeMulti,
  iterations=iterations_sa,
  cooling=cool,
  preserveSatisfied=FALSE,
  penalty=penalty2
)

######### Collect result for Multivariate ############
result_ga_cv_multi <- saAlloc:::.cv2( x=x, strata=ga_solution$indices, sampleSize=ceiling(ga_solution$aggr_strata$SOLUZ))
result_ga_acv_multi <- saAlloc:::.cv2( x=x, strata=ga_solution$indices, sampleSize=ceiling(ga_solution$aggr_strata$SOLUZ), locationAdjustment=x_adj^2)
result_sa_acv_multi <- saAlloc:::.cv2( x=x, strata=sa_solution_acv_multi$label, sampleSize=sa_solution_acv_multi$sampleSize, locationAdjustment=x_adj^2)
result_sa_cv_multi <- saAlloc:::.cv2( x=x, strata=sa_solution_cv_multi$label, sampleSize=sa_solution_cv_multi$sampleSize, locationAdjustment=x_adj^2)
result_multi <- rbind( result_ga_cv_multi,result_ga_acv_multi, result_sa_cv_multi, result_sa_acv_multi)

rownames(result_multi)  <- c("GA (CV)", "GA (ACV)","SA CV-Optimized (ACV)","SA ACV-Optimized (ACV)")
print(result_multi)

