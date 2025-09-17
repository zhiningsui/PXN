library(MatchMixeR)
library(PPCAXNORM)

load("../data/simulated_dataset.RData")

######################################################################
## Analysis I. Compare the performance of PPCA-XPN with MM, OLS for
## platforms A and B using the first set ("AB") of data in Ys
######################################################################

## This function only produces the MSEs for all pairwise predictions
train_predict <- function(X.train, X.test, Y.train, Y.test, K, method="MM", ...) {
  m <- nrow(Y.train)/K
  ##all pairwise combinations of platforms
  combs <- combn(K, 2)
  MSEs <- rep(0, ncol(combs)); Y.target.hat_list <- vector("list", ncol(combs))
  names(Y.target.hat_list) <- apply(combs, 2, paste, collapse="->")
  for (j in 1:ncol(combs)) {
    k.source <- combs[1, j]; k.target <- combs[2, j]
    ## extract Ysource and Ytarget from Y.test
    Y.train.source <- Y.train[((k.source-1)*m+1):(k.source*m),]
    Y.train.target <- Y.train[((k.target-1)*m+1):(k.target*m),]
    Y.test.source <- Y.test[((k.source-1)*m+1):(k.source*m),]
    Y.test.target <- Y.test[((k.target-1)*m+1):(k.target*m),]
    ##
    if (method=="MM") {
      mod <- MM(Y.train.source, Y.train.target)
      b0 <- mod$betamat[, 1]; b1 <- mod$betamat[, 2]
      Y.target.hat <- sweep(sweep(Y.test.source, 1, b1, "*"), 1, b0, "+")
    } else if (method=="OLS") {
      mod <- OLS(Y.train.source, Y.train.target)
      b0 <- mod$betamat[, 1]; b1 <- mod$betamat[, 2]
      Y.target.hat <- sweep(sweep(Y.test.source, 1, b1, "*"), 1, b0, "+")
    } else if (method=="PPCA-NOX") {
      mod <- InitEst(X=NULL, Y=rbind(Y.train.source, Y.train.target), K=2, ...)
      Y.target.hat <- Prediction(Y.test.source, X=NULL, trained.model=mod, k.source=1, k.target=2)
    } else if (method=="PPCA-NOX-GD") {
      mod0 <- InitEst(X=NULL, Y=rbind(Y.train.source, Y.train.target), K=2, ...)
      mod <- GDfun(X=NULL, Y=rbind(Y.train.source, Y.train.target), K=2, mod0)
      Y.target.hat <- Prediction(Y.test.source, X=NULL, trained.model=mod, k.source=1, k.target=2)
    } else if (method=="PPCA-X") {
      mod <- InitEst(X.train, Y=rbind(Y.train.source, Y.train.target), K=2, ...)
      Y.target.hat <- Prediction(Y.test.source, X=X.test, trained.model=mod, k.source=1, k.target=2)
    } else if (method=="PPCA-X-GD") {
      mod0 <- InitEst(X=X.train, Y=rbind(Y.train.source, Y.train.target), K=2, ...)
      mod <- GDfun(X=X.train, Y=rbind(Y.train.source, Y.train.target), K=2, mod0)
      Y.target.hat <- Prediction(Y.test.source, X=X.test, trained.model=mod, k.source=1, k.target=2)
    } else {
      stop("method can only be: MM, OLS, and PPCA.")
    }
    MSEs[j] <- MSE1(Y.target.hat, Y.test.target)
    Y.target.hat_list[[j]] <- Y.target.hat
  }
  return(list(MSEs = MSEs, Y.target.hat_list = Y.target.hat_list))
}

CV <- function(X, Y, K, folds, method="MM", ...) {
  n <- ncol(Y); m <- nrow(Y)/K
  n.folds <- length(folds)
  ##all pairwise combinations of platforms
  combs <- combn(K, 2); MSEs <- matrix(0, n.folds, ncol(combs))
  rownames(MSEs) <- paste0("Fold", 1:n.folds)
  colnames(MSEs) <- apply(combs, 2, paste, collapse="->")
  ##
  Y.target.hat <- matrix(0, nrow=m, ncol=n)
  colnames(Y.target.hat) <- colnames(Y)
  Y.target.hat_list <- replicate(ncol(combs), Y.target.hat, simplify = FALSE)
  names(Y.target.hat_list) <- apply(combs, 2, paste, collapse="->")

  for (k in 1:n.folds){
    test.idx <- folds[[k]]; train.idx <- setdiff(1:n, test.idx)
    Y.test <- Y[, test.idx, drop=FALSE]
    Y.train <- Y[, train.idx, drop=FALSE]
    if (!is.null(X)) {
      X.test <- X[, test.idx, drop=FALSE]
      X.train <- X[, train.idx, drop=FALSE]
    } else {
      X.test <- X.train <- NULL
    }
    rst <- train_predict(X.train, X.test, Y.train, Y.test, K, method=method, ...)
    MSEs[k,] <- rst$MSEs
    for (i in 1:ncol(combs)) {
      Y.target.hat_list[[names(rst$Y.target.hat_list)[i]]][, test.idx] <- rst$Y.target.hat_list[[i]]
    }
  }
  return(list(MSEs = MSEs, Y.target.hat_list = Y.target.hat_list))
}

set.seed(1234)
YA <- Ys$AB[1:m,]
YB <- Ys$AB[(m+1):(2*m),]
n1 <- ns["AB"]
ss10 <- cv.split(n1, n.folds=10)

## YBhat is a list of predicted YB with PPCA and other two methods.
Ls <- seq(0, 15)
YBhat <- list()
MSEs <- list()
# Without X, Without GD
t1 <- system.time(PPCA.NOX <- lapply(Ls, function(L) PPCA=CV(X = Xs$AB, Y = Ys$AB, K = 2, folds = ss10, method="PPCA-NOX", L=L))); t1
MSEs$PPCA.NOX <- lapply(PPCA.NOX, "[[", 1)
YBhat$PPCA.NOX <- lapply(PPCA.NOX, "[[", 2)
# With X, Without GD
t2 <- system.time(PPCA.X <- lapply(Ls, function(L) PPCA=CV(X = Xs$AB, Y = Ys$AB, K = 2, folds = ss10, method="PPCA-X", L=L))); t2
MSEs$PPCA.X <- lapply(PPCA.X, "[[", 1)
YBhat$PPCA.X <- lapply(PPCA.X, "[[", 2)
# Without X, With GD
t3 <- system.time(PPCA.NOX.GD <- lapply(Ls, function(L) PPCA=CV(X = Xs$AB, Y = Ys$AB, K = 2, folds = ss10, method="PPCA-NOX-GD", L=L))); t3
MSEs$PPCA.NOX.GD <- lapply(PPCA.NOX.GD, "[[", 1)
YBhat$PPCA.NOX.GD <- lapply(PPCA.NOX.GD, "[[", 2)
# With X, With GD
t4 <- system.time(PPCA.X.GD <- lapply(Ls, function(L) PPCA=CV(X = Xs$AB, Y = Ys$AB, K = 2, folds = ss10, method="PPCA-X-GD", L=L))); t4

## Compare with MM and OLS
YBhat.MM <- CV(X = NULL, Y = Ys$AB, K = 2, folds = ss10, method="MM")
YBhat.OLS <- CV(X = NULL, Y = Ys$AB, K = 2, folds = ss10, method="OLS")


# Calculate MSEs
MSEs$PPCA.X.GD <- lapply(PPCA.X.GD, "[[", 1)
YBhat$PPCA.X.GD <- lapply(PPCA.X.GD, "[[", 2)
for (ys in names(YBhat)) names(YBhat[[ys]]) <- paste0("L=", Ls)
for (ys in names(MSEs)) names(MSEs[[ys]]) <- paste0("L=", Ls)

MSEs.PPCA <- sapply(MSEs, function(ys) sapply(ys, colMeans))
MSEs.MM <- colMeans(YBhat.MM$MSEs)
MSEs.OLS <- colMeans(YBhat.OLS$MSEs)

pdf("simulations/AB-CVpred2.pdf", width=8, height=8)
yl <- c(min(MSEs.PPCA)-0.002, max(max(MSEs.PPCA), MSEs.MM, MSEs.OLS))
# matplot(Ls, MSEs.PPCA, type="b", col=1, main="", ylim=yl, ylab="MSE", xlab="L")
matplot(Ls, MSEs.PPCA, type="b", lty=c(1,1,2,2), col=1, main="", ylim=yl, ylab="MSE", xlab="L")
abline(v=cvest.AB$Lstar, lty=4)
abline(h=MSEs.MM, col=2)
abline(h=MSEs.OLS, col=4)
legend("bottomright", legend=c("PPCA-NOX", "PPCA-X", "PPCA-NOX-GD", "PPCA-X-GD", "MM", "OLS"),
       lty=c(1,1,2,2,1,1), col=c(1,1,1,1, 2,4), pch=c("1", "2", "3", "4", "", ""))
dev.off()

pdf("simulations/AB-CVpred3.pdf", width=7, height=7)
MSEs.PPCA1 <- MSEs.PPCA[, c("PPCA.X", "PPCA.X.GD")]
yl <- c(min(MSEs.PPCA1)-0.002, max(max(MSEs.PPCA1), MSEs.MM, MSEs.OLS))
matplot(Ls, MSEs.PPCA1, type="b", lty=c(1,2), col=1, main="", ylim=yl, ylab="MSE", xlab="L")
abline(v=cvest.AB$Lstar, lty=4)
## compare with MM and OLS
abline(h=MSEs.MM, col=2)
abline(h=MSEs.OLS, col=4)
legend("bottomright", legend=c("PXN.noGD", "PXN", "MM", "OLS"),
       lty=c(1,2,1,1), col=c(1,1, 2,4), pch=c("1", "2", "", ""))
dev.off()

## As a comparison, the RMSE YA and YB is *very large*, which means
## that without cross-platform normalization, these two datasets are
## not directly comparable.
MSEs.YA <- MSE1(YA, truth=YB)
MSEs.YA

## To select 100 genes for plot "sim1_prediction.pdf"
set.seed(12345)
genes <- sample(1:m, 100)
YtoPlot <- lapply(YBhat$PPCA.NOX, function(ys) as.vector(ys[[paste0("L=", Lstars[["AB"]])]][genes,]))
## orig. data
YtoPlot$YA <- as.vector(YA[genes,]); YtoPlot$YB <- as.vector(YB[genes,])
## MM and OLS
YtoPlot$OLS <- as.vector(YBhat.OLS[genes,]); YtoPlot$MM <- as.vector(YBhat.MM[genes,])

mses1 <- round(c(YA=MSEs.YA, OLS=MSEs.OLS, MM=MSEs.MM, c(MSEs.PPCA[paste0("L=", Lstars[["AB"]], ".1->2"),])), 4)

colnames(YA) <- paste0(colnames(YA), "_A")
colnames(YB) <- paste0(colnames(YB), "_B")

######################################################################
## save useful analytic results for making plots/tables
######################################################################
save(est1, m, X, Beta, K.all, PFs, DS, ns, b0s, b1s, Ls, t1, t2, t3, t4, mses1, MSEs.PPCA, MSEs.MM, MSEs.OLS, YtoPlot, CVMSEs2, best.Ls, ddmat, Ests, genes, Lstars, Lmean, t7, MSEtab3, BCD.MSEs, ss10.all, mses3.nonorm, MSES3, rr4a, rr4b, tabs4a, tabs4b, file="results/sim1-results.rda")



#  Run pairwise normalization for all three pairs ------------------------

load("../data/simulated_dataset.RData")

m = 1000
Ys <- list(
  YA = Ys$AB[1:m,],
  YB = Ys$AB[(m+1):(2*m),],
  YC = Ys$BCD[(m+1):(2*m),],
  YD = Ys$BCD[(2*m+1):(3*m),],
  YE = Ys$EF[1:m,],
  YF = Ys$EF[(m+1):(2*m),]
)

PFs <- list(
  "AB" = c("YA", "YB"),
  "BC" = c("YB", "YC"),
  "BD" = c("YB", "YD"),
  "CD" = c("YC", "YD")
)

MSEs_Sh2 <- c()
colnames(MSEs_Sh2) <- names(Ys)

for (dp in PFs) {
  dpn <- paste(dp, collapse="_")
  Y.source <- Ys[[dp[1]]]
  Y.target <- Ys[[dp[2]]]
  n <- ncol(Y1)

  output_file <- paste0("Harmonized_", dpn, ".rda")

  Y.target.hat <- Shambhala2(
    InputData = Y.source,
    PData = Y.target,
    QData = Y.target,
    delete_buffer_files = TRUE,
    k = 5
  )
  save(Y.target.hat, file = output_file)

  Y.target.hat <- Y.target.hat %>%
    column_to_rownames("SYMBOL")

  MSEs_Sh2[dpn] <- MSE1(Y.target.hat, Y.target)

}

rownames(MSEs_Sh2) <- paste0("Fold", 1:nrow(MSEs_Sh2))


Shambhala2 <- function( InputData, PData, QData, delete_buffer_files = TRUE, k = 5 ) {

  P = PData # log2
  cat("Range of P:", range(P), "\n")

  MAS0 = as.matrix(InputData) # log2
  cat("Range of MAS0:", range(MAS0), "\n")

  NS = ncol(MAS0)
  SYMBOL = rownames(MAS0)
  CN = colnames(MAS0)

  pool = merge(MAS0, P, by = 0) # log2
  colnames(pool)[1] <- "SYMBOL"

  NS1 = ncol(pool)
  NG1 = nrow(pool)

  pool[, 2:ncol(pool)] <- dplyr::mutate_all(pool[, 2:ncol(pool)], function(x) as.numeric(as.character(x)))

  P1FN = "P_prim.txt"
  write.table(pool, P1FN, row.names = FALSE, col.names = TRUE, sep = "\t")

  NH = ncol(MAS0)
  NP = ncol(P)

  args = c(NH,NP,k)
  AFN = "args.txt"
  write.table(args, AFN, row.names = FALSE, col.names = FALSE)

  system("/Applications/MATLAB_R2024b.app/bin/matlab -nodesktop -nosplash -nodisplay -r \"cd('/Users/zhiningsui/GitHub/PPCAXNORM/simulations');run('Shambhala2.m');exit;\"")

  Q = QData # log2
  cat("Range of Q:", range(Q), "\n")

  Cu2FN = "Cu_bis.txt"
  MAS = read.table(Cu2FN, header = FALSE, sep = " ") # original
  MAS = as.matrix(MAS) # original

  MAS1 = merge(MAS,Q, by.x = 1, by.y = 0) # 2: (NS+1) original, (NS+2):ncol(MAS1) log2
  MAS1[, 2:ncol(MAS1)] <- dplyr::mutate_all(MAS1[, 2:ncol(MAS1)], function(x) as.numeric(as.character(x)))

  MAS3 = MAS1[, (NS+2):ncol(MAS1)] # log2
  cat("Range of MAS3:", range(MAS3), "\n")

  RM = rowMeans(as.matrix(MAS3))
  RS = matrixStats::rowSds(as.matrix(MAS3))

  MAS2 = MAS1[, 2:(NS+1)] # original
  cat("Range of MAS2:", range(MAS2), "\n")
  MAS22 = log2(MAS2+1) # log2
  cat("Range of MAS22:", range(MAS22), "\n")

  NR = nrow(MAS22)
  for ( nr in 1:NR ) {
    MAS22[nr,] = RM[nr] + RS[nr]*MAS22[nr,]
  }

  # MAS23 = exp(MAS22)
  MAS23 = MAS22
  cat("Range of MAS23:", range(MAS23), "\n")

  SYMBOL = as.vector(MAS1[,1])
  MAS33 = cbind(SYMBOL, MAS23)

  for ( j in 2:(NS+1) ) {
    MAS33[,j] = as.vector(as.numeric(MAS33[,j]))
  }

  if ( delete_buffer_files ) {

    if (file.exists(P1FN)) {
      #Delete file if it exists
      file.remove(P1FN)
    }

    if (file.exists(Cu2FN)) {
      #Delete file if it exists
      file.remove(Cu2FN)
    }

    if (file.exists(AFN)) {
      #Delete file if it exists
      file.remove(AFN)
    }

  }

  colnames(MAS33) = c("SYMBOL", CN)

  return(MAS33)

}



YBhat.Sh2 <- Shambhala2(InputData = YA, PData = YB, QData = YB)
