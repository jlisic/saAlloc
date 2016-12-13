library(SamplingStrata)

# load test data set
data("swissmunicipalities", package = "sampling")


frame <- NULL
frame$id <- swissmunicipalities$Nom


frame$Y1 <- swissmunicipalities$Pop020
frame$Y2 <- swissmunicipalities$Pop2040
frame$Y3 <- swissmunicipalities$Pop4065
frame$Y4 <- swissmunicipalities$Pop65P

library("SamplingStrata")
set.seed(1508)
frame$X1 <- var.bin(swissmunicipalities$POPTOT, bins = 18)
frame$X2 <- var.bin(swissmunicipalities$Surfacesbois, bins = 3)
frame$X3 <- var.bin(swissmunicipalities$Surfacescult, bins = 3)
frame$X4 <- var.bin(swissmunicipalities$Alp, bins = 3)
frame$X5 <- var.bin(swissmunicipalities$Airbat, bins = 3)
frame$X6 <- var.bin(swissmunicipalities$Airind, bins = 3)

frame$domainvalue <- swissmunicipalities$REG
frame <- as.data.frame(frame)

# build atomic strata
strata <- buildStrataDF(frame)

cv <- data.frame(DOM = "DOM1", CV1 = 0.05, CV2 = 0.05, CV3 = 0.05, CV4 = 0.05, domainvalue = 1:7)


errors <- cv[1, 1:5]
allocation <- bethel(strata, errors)
length(allocation)



solution <- optimizeStrata(errors = cv, strata = strata, cens = NULL,
    strcens = FALSE, initialStrata = nrow(strata), addStrataFactor = 0.00,
    minnumstr = 2, iter = 400, pops = 20, mut_chance = 0.005,
    elitism_rate = 0.2, highvalue = 1e+08, suggestions = NULL,
    realAllocation = TRUE, writeFiles = TRUE)


samp <- read.delim("outstrata.txt")



library(saAlloc)

#solution_saAlloc <- saMinCV(
#  x[,c("X1","X2" 



