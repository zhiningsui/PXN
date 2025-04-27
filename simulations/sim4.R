library(MatchMixeR)
library(PPCAXNORM)

load("data/simulated_dataset.RData")

######################################################################
## Analysis IVa: Using the AB data (contains subsets A1 and B1) to
## train a PPCA-XPN model. Use it to normalize data A1 to
## \hat{B1}. Use an LMER to regress y~treat+age+(1|ID) for each
## gene. Note that we need the random intercept to account for the
## fact that A1 and B1 are repeated measures of the same biological
## samples.
######################################################################
## these two indicator vectors are required for computing TPR/FPR and AUC
DEGs.treat.vec <- rep(0,m)
DEGs.treat.vec[DEGs.treat] <- 1
DEGs.age.vec <- rep(0,m)
DEGs.age.vec[DEGs.age] <- 1

## H-T function. XAB is the covariate matrix (must include ID for LMER)
HT <- function(YA, YB, XAB, method=c("LMER", "OLS", "pcomb"), alpha=0.05, percentage=TRUE) {
  method <- match.arg(method)
  n <- ncol(YA); m <- nrow(YA)
  YAB <- cbind(YA, YB)
  if (method=="LMER") {
    ## ID is the unique identifier for biological samples
    pp <- t(sapply(1:m, function(i) {
      data.i <- data.frame(y=YAB[i,], XAB)
      mod.i <- suppressMessages(lmer(y~treat+age+(1|ID), data=data.i))
      return(coef(summary(mod.i))[-1,"Pr(>|t|)"])
    }))
  } else if (method=="OLS") {
    pp <- t(sapply(1:m, function(i) {
      data.i <- data.frame(y=YAB[i,], XAB)
      mod.i <- lm(y~treat+age, data=data.i)
      return(coef(summary(mod.i))[-1,"Pr(>|t|)"])
    }))
  } else if (method=="pcomb") {
    pp <- t(sapply(1:m, function(i) {
      data.i <- data.frame(y=YAB[i,], XAB)
      mod.i.a <- lm(y~treat+age, data=data.i[1:n,])
      mod.i.b <- lm(y~treat+age, data=data.i[(n+1):(2*n),])
      pp.a <- coef(summary(mod.i.a))[-1,"Pr(>|t|)"]
      pp.b <- coef(summary(mod.i.b))[-1,"Pr(>|t|)"]
      ## p-value combination based on Fisher's method
      return(sapply(names(pp.a), function(nn) fisher(c(pp.a[nn], pp.b[nn]))))
    }))
  }
  ## summarize the results
  sel.treat <- which(pp[,"treat"]<alpha)
  sel.age <- which(pp[,"age"]<alpha)
  tp.treat <- length(intersect(DEGs.treat, sel.treat))
  fp.treat <- length(sel.treat)-tp.treat
  tp.age <- length(intersect(DEGs.age, sel.age))
  fp.age <- length(sel.age)-tp.age
  tpfp <- c(tp.treat=tp.treat, tp.age=tp.age, fp.treat=fp.treat, fp.age=fp.age)
  aucs <- c(auc.treat=auc(DEGs.treat.vec, -pp[,"treat"]),
            auc.age=auc(DEGs.age.vec, -pp[,"age"]))
  if (percentage) {
    tpfp <- 100*tpfp/c(m1.treat, m1.age, m-m1.treat, m-m1.age)
    aucs <- 100*aucs
  }
  return(c(tpfp, aucs))
}

########## the main loop for H-T ##########
nreps <- 100
## create a covariate matrix for HT()
ID <- rep(1:ns["AB"], 2); XAB <- data.frame(ID=ID, t(Xs$AB)[ID,])
rr4a <- rr4b <- list()
## it takes about 40 minutes
t4 <- system.time( for (i in 1:nreps) {
  ## Generate the overall data (6000x240) and put a selected subset of
  ## into a list of three subsets: (AB, BCD, and EF)
  Y.all <- datagen(X, Beta, A, sigma2s, sigmaU2s, b0.all, b1.all)
  ## select subsets of Y to serve as datasets AB, BCD, and EF
  slist.u <- cumsum(ns); slist.l <- slist.u-ns+1
  slist <- sapply(1:length(ns), function(k) seq(slist.l[k], slist.u[k]))
  names(slist) <- DS
  glist0 <- rep(PF.all, each=m)
  glist <- lapply(DS, function(ds) which(glist0 %in% PFs[[ds]])); names(glist) <- DS
  Ys <- lapply(DS, function(ds) {
    Y.all[glist[[ds]], slist[[ds]]]
  }); names(Ys) <- DS

  # Analysis 4a.
  mod0 <- InitEst(X=Xs$AB, Y=Ys$AB, K=2, L=Lstars[["AB"]])
  mod1 <- GDfun(X=Xs$AB, Y=Ys$AB, K=2, mod0)
  mod0.nox <- InitEst(X=NULL, Y=Ys$AB, K=2, L=Lstars[["AB"]])
  mod1.nox <- GDfun(X=NULL, Y=Ys$AB, K=2, mod0.nox)
  YA1 <- Ys$AB[1:m,]; YB1 <- Ys$AB[(m+1):(2*m),]
  YBhat <- Prediction(YA1, X=Xs$AB, trained.model=mod1, k.source=1, k.target=2)
  YBhat.nox <- Prediction(YA1, X=NULL, trained.model=mod1.nox, k.source=1, k.target=2)
  mod1.MM <- MM(YA1, YB1)
  b0 <- mod1.MM$betamat[, 1]; b1 <- mod1.MM$betamat[, 2]
  YBhat.MM <- sweep(sweep(YA1, 1, b1, "*"), 1, b0, "+")
  # YBhat.Sh2.B <- Shambhala2(input_mtx = YA1, P_mtx = YB1, Q_mtx = YB1)
  # YAB <- cbind(YA1, YB1)
  # colnames(YAB) <- paste0("s",1:ncol(YAB))
  # YBhat.Sh2.AB <- Shambhala2(input_mtx = YA1, P_mtx = YAB, Q_mtx = YB1)

  ## it takes about 2~3 minutes to finish one iteration
  MM4a <- c("OLS", "LMER", "pcomb"); names(MM4a) <- MM4a
  rr4a[[i]] <- list(orig=sapply(MM4a, function(mm) HT(YA1, YB1, XAB, method=mm)),
                    ppca=sapply(MM4a, function(mm) HT(YBhat, YB1, XAB, method=mm)),
                    ppca.nox=sapply(MM4a, function(mm) HT(YBhat.nox, YB1, XAB, method=mm)),
                    # sh2=sapply(MM4a, function(mm) HT(YBhat.Sh2.B, YB1, XAB, method=mm)),
                    # sh2.AB=sapply(MM4a, function(mm) HT(YBhat.Sh2.AB, YB1, XAB, method=mm))
                    MM=sapply(MM4a, function(mm) HT(YBhat.MM, YB1, XAB, method=mm))
  )

  ## Analysis 4b.
  modL <- lapply(names(Xs), function(dn) {
    mod0 <- InitEst(Xs[[dn]], Ys[[dn]], K=Ks[[dn]], L=Lmean)
    mod1 <- GDfun(Xs[[dn]], Ys[[dn]], K=Ks[[dn]], mod0)
    return(mod1)
  })
  ModAll <- ModIntegrate(modL, Xs, Ys, PFs)
  ## nox version
  modL.nox <- lapply(names(Xs), function(dn) {
    mod0 <- InitEst(NULL, Ys[[dn]], K=Ks[[dn]], L=Lmean)
    mod1 <- GDfun(NULL, Ys[[dn]], K=Ks[[dn]], mod0)
    return(mod1)
  })
  ModAll.nox <- ModIntegrate(modL.nox, NULL, Ys, PFs)

  ## Generate new data, YA4 and YB4 (n=30); YA5 and YB5 (n=50).
  Ynew <- datagen(Xnew, Beta, A, sigma2s, sigmaU2s, b0.all, b1.all)
  ## We only observe YA4 and YB5
  YA4 <- Ynew[1:m, 1:n4]; YB5 <- Ynew[(m+1):(2*m), (n4+1):n]
  ## 3. Apply PPCA-XPN to estimate YB4 from YA4
  YB4hat <- Prediction(YA4, X=X4, trained.model=ModAll, k.source=1, k.target=2)
  YB4hat.nox <- Prediction(YA4, X=NULL, trained.model=ModAll.nox, k.source=1, k.target=2)
  YB4hat.MM <- sweep(sweep(YA4, 1, b1, "*"), 1, b0, "+")
  # YB4hat.Sh2.A <- Shambhala2(input_mtx = YA4, P_mtx = YA1, Q_mtx = YB5)

  # Y.tmp <- list()
  # for (mat in Ys) {
  #   Y.split <- lapply(seq(1, nrow(mat), by = 1000), function(i) mat[i:(i+999), ])
  #   Y.tmp[[length(Y.tmp) + 1]] <- do.call(cbind, Y.split)
  # }
  # Y <- do.call(cbind, Y.tmp); colnames(Y) <- paste0("s", 1:ncol(Y))
  # YA.all <- Ys[["AB"]][1:m, ]
  # YB.all <- cbind(Ys[["AB"]][(m+1):(2*m), ], Ys[["BCD"]][1:m, ])
  # YC.all <- Ys[["BCD"]][1001:2000, ]
  # YD.all <- Ys[["BCD"]][2001:3000, ]
  # YE.all <- Ys[["EF"]][1:1000, ]
  # YF.all <- Ys[["EF"]][1001:2000, ]
  # YB4hat.Sh2 <- Shambhala2(input_mtx = YA4, P_mtx = Y, Q_mtx = YB.all)
  # YB4hat.Sh2.B <- Shambhala2(input_mtx = YA4, P_mtx = YB.all, Q_mtx = YB.all)
  # YB4hat.Sh2.C <- Shambhala2(input_mtx = YA4, P_mtx = YC.all, Q_mtx = YB.all)
  # YB4hat.Sh2.D <- Shambhala2(input_mtx = YA4, P_mtx = YD.all, Q_mtx = YB.all)
  # YB4hat.Sh2.E <- Shambhala2(input_mtx = YA4, P_mtx = YE.all, Q_mtx = YB.all)
  # YB4hat.Sh2.F <- Shambhala2(input_mtx = YA4, P_mtx = YF.all, Q_mtx = YB.all)


  ## H-T. For Analysis 4b, no pairing --> no need to use LMER.
  MM4b <- c("OLS", "pcomb"); names(MM4b) <- MM4b
  rr4b[[i]] <- list(orig=sapply(MM4b, function(mm) HT(YA4, YB5, t(Xnew), method=mm)),
                    ppca=sapply(MM4b, function(mm) HT(YB4hat, YB5, t(Xnew), method=mm)),
                    ppca.nox=sapply(MM4b, function(mm) HT(YB4hat.nox, YB5, t(Xnew), method=mm)),
                    # sh2=sapply(MM4b, function(mm) HT(YB4hat.Sh2, YB5, t(Xnew), method=mm)),
                    # sh2=sapply(MM4b, function(mm) HT(YB4hat.Sh2.A, YB5, t(Xnew), method=mm)),
                    # sh2.B=sapply(MM4b, function(mm) HT(YB4hat.Sh2.B, YB5, t(Xnew), method=mm)),
                    # sh2.C=sapply(MM4b, function(mm) HT(YB4hat.Sh2.C, YB5, t(Xnew), method=mm)),
                    # sh2.D=sapply(MM4b, function(mm) HT(YB4hat.Sh2.D, YB5, t(Xnew), method=mm)),
                    # sh2.E=sapply(MM4b, function(mm) HT(YB4hat.Sh2.E, YB5, t(Xnew), method=mm)),
                    # sh2.F=sapply(MM4b, function(mm) HT(YB4hat.Sh2.F, YB5, t(Xnew), method=mm))
                    MM=sapply(MM4b, function(mm) HT(YB4hat.MM, YB5, t(Xnew), method=mm))
  )
  ##
  print(paste0("Done with iteration ", i, "."))
})

## summarize rr4a and rr4b into more readable tables
tabs4a <- lapply(names(rr4a[[1]]), function(mm) {
  Reduce("+", lapply(rr4a, function(rr) rr[[mm]]))/nreps
}); names(tabs4a) <- names(rr4a[[1]])
##
tabs4b <- lapply(names(rr4b[[1]]), function(mm) {
  Reduce("+", lapply(rr4b, function(rr) rr[[mm]]))/nreps
}); names(tabs4b) <- names(rr4b[[1]])

lapply(tabs4a, round, digits=1)
lapply(tabs4b, round, digits=1)

