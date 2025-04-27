library(MatchMixeR)
library(PPCAXNORM)

load("data/simulated_dataset.RData")

######################################################################
## Analysis III. Model integration
######################################################################

## Train individual PPCA models with a common, weighted average L.
## Overall (in terms of the estimated AA'), the results are better than using Lmin.
Lmean <- round(sum(ns*Lstars)/sum(ns))
##
mod0.AB <- InitEst(Xs$AB, Ys$AB, K=2, L=Lmean)
mod.AB <- GDfun(Xs$AB, Ys$AB, K=2, mod0.AB)
##
mod0.BCD <- InitEst(Xs$BCD, Ys$BCD, K=3, L=Lmean)
mod.BCD <- GDfun(Xs$BCD, Ys$BCD, K=3, mod0.BCD)
##
mod0.EF <- InitEst(Xs$EF, Ys$EF, K=2, L=Lmean)
mod.EF <- GDfun(Xs$EF, Ys$EF, K=2, mod0.EF)

######################################################################
## Integrate three models into one grand model
######################################################################

## modList is a list of trained PPCA models (with the same L)
modL <- list(AB=mod.AB, BCD=mod.BCD, EF=mod.EF)
t7 <- system.time(ModAll <- ModIntegrate(modList = modL, Xlist = Xs, Ylist = Ys, PairList = PFs))

######################################################################
## check the accuracy of parameter estimation
######################################################################
b0.all <- cbind(b0s$AB, b0s$BCD[, 2:3], b0s$EF)
b1.all <- cbind(b1s$AB, b1s$BCD[, 2:3], b1s$EF)
AA <- est1$U %*% diag(est1$dd) %*% t(est1$U) #true AA' matrix
Grps <- names(modL); L <- length(ModAll$dd)

## accuracy of AA' needs some special care
mse.int <- c(b0=MSE1(ModAll$b0, b0.all),
             b1=MSE1(ModAll$b1, b1.all),
             Beta=MSE1(ModAll$Beta, Beta),
             sigmaU2s=MSE1(ModAll$sigmaU2s, est1$sigmaU2s),
             sigma2s=MSE1(ModAll$sigma2s, est1$sigma2s),
             d2=MSE1(ModAll$dd, est1$dd[1:L]),
             AA=mean( ((ModAll$U %*% diag(ModAll$dd))%*% t(ModAll$U) -AA)^2 ) )
##
MSEtab3 <- data.frame(b0=sapply(Grps, function(gg) MSE1(modL[[gg]]$b0, b0s[[gg]])),
                      b1=sapply(Grps, function(gg) MSE1(modL[[gg]]$b1, b1s[[gg]])),
                      Beta=sapply(Grps, function(gg) MSE1(modL[[gg]]$Beta, Beta)),
                      sigmaU2s=sapply(Grps, function(gg) MSE1(modL[[gg]]$sigmaU2s, est1$sigmaU2s)),
                      sigma2s=sapply(Grps, function(gg) MSE1(modL[[gg]]$sigma2s, est1$sigma2s)),
                      d2=sapply(Grps, function(gg) MSE1(modL[[gg]]$dd[1:L], est1$dd[1:L])),
                      AA=sapply(Grps, function(gg) {
                        U <- modL[[gg]]$U; dd <- modL[[gg]]$dd
                        AA.g <- U%*%diag(dd)%*%t(U)
                        return(mean((AA.g -AA)^2))
                      }))
##
MSEtab3 <- rbind(MSEtab3, average=colMeans(MSEtab3))
MSEtab3 <- rbind(MSEtab3, Integrated=mse.int)
round(MSEtab3, 4)

kableExtra::kbl(MSEtab3, digits = 4, format = "latex")


######################################################################
## Predictions
##
## We want to demonstrate that: (a) predictions related to platform B
## should be more accurate based on the integrated model than both
## mod.AB and mod.BCD (we study the predictions A1-->B1 and B2-->C2,
## with uni-models and integrated model), and (b) the integrated model
## is capable of making predictions between platform A2-->C2 (no
## direct connection) and A3-->F3 (no direct nor indirect
## connections), which are very difficult or impossible with
## MatchMixer or individually trained PPCA-Xnorm models. In addition,
## the integrated model should be able to improve the prediction
## accuracy for E3-->F3, even though models AB and BCD are
## *independent* of model EF.
##
## As a remark, we do not need to compare the integrated PPCA-Xnorm
## model with MatchMixeR, because we have already demonstrated that
## PPCA-Xnorm is better at making predictions than MatchMixeR.
####################################################################

## 10-fold CV using all three sets of data. To make the splits more
## even, we conduct the train/test split for each dataset separately.
set.seed(22222)
nFolds <- 10
ss10.all <- lapply(ns, function(n) cv.split(n, n.folds=nFolds))

YA1 <- Ys[["AB"]][1:m,]
YA2 <- Y.all[1:m, seq(slist.l["BCD"], slist.u["BCD"])]
YA3 <- Y.all[1:m, seq(slist.l["EF"], slist.u["EF"])]
YB1 <- Ys[["AB"]][(m+1):(2*m),]
YB2 <- Ys[["BCD"]][1:m,]
YC2 <- Ys[["BCD"]][(m+1):(2*m),]
YD2 <- Ys[["BCD"]][(2*m+1):(3*m),]
YE3 <- Ys[["EF"]][1:m,]
YF3 <- Ys[["EF"]][(m+1):(2*m),]

preds3 <- c("A1->B1 (AB)", "A1->B1 (Int1)", "A1->B1 (Int2)",
            "A2->C2 (BCD)", "A2->C2 (Int1)", "A2->C2 (Int2)",
            "B2->C2 (BCD)", "B2->C2 (Int1)", "B2->C2 (Int2)",
            "A3->F3 (EF)", "A3->F3 (Int1)", "A3->F3 (Int2)",
            "E3->F3 (EF)", "E3->F3 (Int1)", "E3->F3 (Int2)")

## Before the start of the 10-fold CV, let us compute the prediction
## MSE when *no normalization* is applied. For example, A1->B1 without
## normalization means to use A1 as the prediction of B1.
mses3.nonorm <- c("A1->B1"=MSE1(YA1, YB1),
                  "A2->C2"=MSE1(YA2, YC2),
                  "B2->C2"=MSE1(YB2, YC2),
                  "A3->F3"=MSE1(YA3, YF3),
                  "E3->F3"=MSE1(YE3, YF3))
round(mses3.nonorm,digits = 4)
## the 10-fold CV loop
MSES3 <- matrix(NA, length(preds3), nFolds); rownames(MSES3) <- preds3
for (k in 1:nFolds){
  test.idx <- lapply(ss10.all, function(ff) ff[[k]]); names(test.idx) <- DS
  train.idx <- lapply(DS, function(nn) setdiff(1:ns[nn], test.idx[[nn]])); names(train.idx) <- DS
  Xs.test <- lapply(DS, function(nn) Xs[[nn]][, test.idx[[nn]], drop=FALSE]); names(Xs.test) <- DS
  Xs.train <- lapply(DS, function(nn) Xs[[nn]][, train.idx[[nn]], drop=FALSE]); names(Xs.train) <- DS
  Ys.test <- lapply(DS, function(nn) Ys[[nn]][, test.idx[[nn]], drop=FALSE]); names(Ys.test) <- DS
  Ys.train <- lapply(DS, function(nn) Ys[[nn]][, train.idx[[nn]], drop=FALSE]); names(Ys.train) <- DS
  ## train three individual models and then integrate them.
  mods <- lapply(DS, function(nn) {
    mod0 <- InitEst(Xs.train[[nn]], Ys.train[[nn]], K=Ks[nn], L=Lmean)
    mod <- GDfun(Xs.train[[nn]], Ys.train[[nn]], K=Ks[nn], mod0)
    return(mod)
  }); names(mods) <- DS
  mod.int1 <- ModIntegrate(mods, Xs.train, Ys.train, PFs)
  ## Now, pretend that all platforms are independent
  YT <- list(A=YA1[,train.idx[["AB"]]],
             B=cbind(YB1[,train.idx[["AB"]]], YB2[,train.idx[["BCD"]]]),
             C=YC2[,train.idx[["BCD"]]], D=YD2[,train.idx[["BCD"]]],
             E=YE3[,train.idx[["EF"]]], F=YF3[,train.idx[["EF"]]])
  XT <- list(A=Xs.train[["AB"]],
             B=cbind(Xs.train[["AB"]], Xs.train[["BCD"]]),
             C=Xs.train[["BCD"]], D=Xs.train[["BCD"]],
             D=Xs.train[["EF"]], F=Xs.train[["EF"]])
  mods.nopair <- lapply(names(YT), function(pf) {
    mod0 <- InitEst(XT[[pf]], YT[[pf]], K=1, L=Lmean)
    return(GDfun(XT[[pf]], YT[[pf]], K=1, mod0))
  })
  mod.int2 <- ModIntegrate(mods.nopair, XT, YT, PairList=names(YT))
  ## making predictions and compute the MSEs
  ##
  ## A1-->B1, Model AB or the integrated model
  YA1.test <- YA1[,test.idx[["AB"]]]
  YB1.test <- YB1[,test.idx[["AB"]]]
  Yhat.B1 <- Prediction(YA1.test, Xs.test[["AB"]], mods[["AB"]], 1, 2)
  Yhat.B1.int1 <- Prediction(YA1.test, Xs.test[["AB"]], mod.int1, 1, 2)
  Yhat.B1.int2 <- Prediction(YA1.test, Xs.test[["AB"]], mod.int2, 1, 2)
  MSES3["A1->B1 (AB)",k] <- MSE1(Yhat.B1, YB1.test)
  MSES3["A1->B1 (Int1)",k] <- MSE1(Yhat.B1.int1, YB1.test)
  MSES3["A1->B1 (Int2)",k] <- MSE1(Yhat.B1.int2, YB1.test)
  ## A2-->C2
  YA2.test <- YA2[,test.idx[["BCD"]]]
  YC2.test <- YC2[,test.idx[["BCD"]]]
  ## in mod.int1/int2, platfroms A and C are platforms 1 and 3
  Yhat.C2.int1.A2 <- Prediction(YA2.test, Xs.test[["BCD"]], mod.int1, 1, 3)
  Yhat.C2.int2.A2 <- Prediction(YA2.test, Xs.test[["BCD"]], mod.int2, 1, 3)
  MSES3["A2->C2 (Int1)",k] <- MSE1(Yhat.C2.int1.A2, YC2.test)
  MSES3["A2->C2 (Int2)",k] <- MSE1(Yhat.C2.int2.A2, YC2.test)
  ## B2-->C2
  YB2.test <- YB2[,test.idx[["BCD"]]]
  Yhat.C2 <- Prediction(YB2.test, Xs.test[["BCD"]], mods[["BCD"]], 1, 2)
  Yhat.C2.int1.B2 <- Prediction(YB2.test, Xs.test[["BCD"]], mod.int1, 2, 3)
  Yhat.C2.int2.B2 <- Prediction(YB2.test, Xs.test[["BCD"]], mod.int2, 2, 3)
  MSES3["B2->C2 (BCD)",k] <- MSE1(Yhat.C2, YC2.test)
  MSES3["B2->C2 (Int1)",k] <- MSE1(Yhat.C2.int1.B2, YC2.test)
  MSES3["B2->C2 (Int2)",k] <- MSE1(Yhat.C2.int2.B2, YC2.test)
  ## A3->F3
  YA3.test <- YA3[,test.idx[["EF"]]]
  YF3.test <- YF3[,test.idx[["EF"]]]
  Yhat.F3.int1.A3 <- Prediction(YA3.test, Xs.test[["EF"]], mod.int1, 1, 6)
  Yhat.F3.int2.A3 <- Prediction(YA3.test, Xs.test[["EF"]], mod.int2, 1, 6)
  MSES3["A3->F3 (Int1)",k] <- MSE1(Yhat.F3.int1.A3, YF3.test)
  MSES3["A3->F3 (Int2)",k] <- MSE1(Yhat.F3.int2.A3, YF3.test)
  ## E3->F3
  YE3.test <- YE3[,test.idx[["EF"]]]
  Yhat.F3 <- Prediction(YE3.test, Xs.test[["EF"]], mods[["EF"]], 1, 2)
  Yhat.F3.int1.E3 <- Prediction(YE3.test, Xs.test[["EF"]], mod.int1, 5, 6)
  Yhat.F3.int2.E3 <- Prediction(YE3.test, Xs.test[["EF"]], mod.int2, 5, 6)
  MSES3["E3->F3 (EF)",k] <- MSE1(Yhat.F3, YF3.test)
  MSES3["E3->F3 (Int1)",k] <- MSE1(Yhat.F3.int1.E3, YF3.test)
  MSES3["E3->F3 (Int2)",k] <- MSE1(Yhat.F3.int2.E3, YF3.test)
}
as.matrix(rowMeans(MSES3))

kableExtra::kbl(matrix(rowMeans(MSES3),nrow = 3), digits = 4, format = "latex")

