library(MatchMixeR)
library(PPCAXNORM)

load("data/simulated_dataset.RData")

######################################################################
## Analysis II. Using the second set of data (BCD, n=60). Focusing on
## the accuracy of parameter estimation.
######################################################################

######################################################################
## Using10-fold CV to confirm that PPCA-GD is also the best method
######################################################################

rr <- list()
ss10 <- cv.split(ns["BCD"], n.folds=10)
rr$MM <- CV(Xs$BCD, Ys$BCD, K=3, ss10, method="MM")
rr$OLS <- CV(Xs$BCD, Ys$BCD, K=3, ss10, method="OLS")
rr$PPCA.NOX <- CV(Xs$BCD, Ys$BCD, K=3, ss10, method="PPCA-NOX", L=Lstars[["BCD"]])
rr$PPCA.X <- CV(Xs$BCD, Ys$BCD, K=3, ss10, method="PPCA-X", L=Lstars[["BCD"]])
rr$PPCA.NOX.GD <- CV(Xs$BCD, Ys$BCD, K=3, ss10, method="PPCA-NOX-GD", L=Lstars[["BCD"]])
rr$PPCA.X.GD <- CV(Xs$BCD, Ys$BCD, K=3, ss10, method="PPCA-X-GD", L=Lstars[["BCD"]])

##
BCD.MSEs <- t(sapply(lapply(rr, "[[", 1), colMeans))
BCD.MSEs
kableExtra::kbl(BCD.MSEs, digits = 4, format = "latex")

MSE1(Ys$BCD[1:1000,], truth=Ys$BCD[1001:2000,])
MSE1(Ys$BCD[1:1000,], truth=Ys$BCD[2001:3000,])
MSE1(Ys$BCD[1001:2000,], truth=Ys$BCD[2001:3000,])


######################################################################
## Now use all data (BCD) to produce estimates, and check the accuracy
## of parameter estimations
######################################################################
Ests <- list()
t4 <- system.time(Ests$PPCA.NOX <- InitEst(X=NULL, Ys$BCD, K=3, Lstars[["BCD"]])); t4
t5 <- system.time(Ests$PPCA.X <- InitEst(Xs$BCD, Ys$BCD, K=3, Lstars[["BCD"]])); t5
## GD refinement
t6 <- system.time(Ests$PPCA.NOX.GD <- GDfun(X=NULL, Ys$BCD, K=3, Ests$PPCA.NOX)); t6
t7 <- system.time(Ests$PPCA.X.GD <- GDfun(Xs$BCD, Ys$BCD, K=3, Ests$PPCA.X)); t7

## check dd
L.oracle <- length(est1$dd)
ddest <- sapply(Ests, function(ee) ee$dd)
ddmat <- cbind(Oracle=est1$dd, rbind(ddest, matrix(NA, L.oracle-Lstars[["BCD"]], ncol(ddest))))
rownames(ddmat) <- paste0("L", 1:L.oracle)
## For our own understanding: 1. As of 05/10/2024, GD does not update
## covariance-related parameters (including dd, U, sigmaU2s, and
## sigma2s). 2. Using the covariates (X) improves the accuracy of
## estimating dd.
ddmat

## check the estimates of b1. As of 05/10/2024, b0 is not updated in
## GD nor with the inclusion of X (centered) so no need to compare
## between the three methods.
b0.mses <- sapply(Ests, function(ee) MSE1(ee$b0, b0s$BCD))
b0.mses                                 #for our own understanding

b1.mses <- sapply(Ests, function(ee) MSE1(ee$b1, b1s$BCD))
b1.mses                                 #GD is the best

## Visualization. For clarity, only use a subset of randomly selected
## 200 genes in plots.
sigma2s.mses <- sapply(Ests, function(ee) MSE1(ee$sigma2s, est1$sigma2s))
sigmaU2s.mses <- sapply(Ests, function(ee) MSE1(ee$sigmaU2s, est1$sigmaU2s))
## Using X, we have better estimates of sigma2s and sigmaU2s. Note
## that GD does not change these variables.
rbind(sigma2s.mses, sigmaU2s.mses)

## check Beta
Beta.mses <- c(PPCA.X=MSE1(Ests$PPCA.X$Beta, Beta),
               PPCA.X.GD=MSE1(Ests$PPCA.X.GD$Beta, Beta))
Beta.mses                               #GD improves Beta estimation

kableExtra::kbl(rbind(b0.mses, b1.mses, sigma2s.mses, sigmaU2s.mses, Beta.mses),
                digits = 4,
                format = "latex")
