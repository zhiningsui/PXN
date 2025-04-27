######################################################################
## The X matrix is designed so that x1 (treat) is binary, x2 (age) and
## x3 are continuous.
## beta2 is nonzero for only for a subset of 100 genes.
######################################################################

## Use the estimates from array.seq of MatchMixeR as the "oracle" for simulation
load("data/est.array.seq.rda")

est1 <- ests$array.seq
sigmaU2s <- est1$sigmaU2s
sigma2s <- est1$sigma2s
U <- est1$U
dd <- est1$dd
L <- length(dd)
m <- nrow(est1$b0)
gnames <- rownames(U)

## Compute the oracle A
A <- sweep(U, 2, sqrt(dd), "*")
# Check if AA' +sigmaU2s +sigma2s is a correlation matrix?
rowsums(head(A)^2) + head(sigmaU2s) + head(sigma2s) # Yes

## We will use the PPCA covariance structure estimated from array.seq,
## together with simulated b0, b1, X, and Beta for data generated on 6
## platforms. Note that "B" is shared by the first and second set of
## matched samples. However, samples in the first set (n=40) and the
## second set (n=60) are independent.
PFs <- list("AB"=c("A", "B"),
            "BCD"=c("B", "C", "D"),
            "EF"=c("E", "F"))
DS <- names(PFs) # Names of the datasets
Ks <- sapply(PFs, length)
PF.all <- Reduce("union", PFs) # A...F
K.all <- length(PF.all)   # 6
ns <- c("AB"=60, "BCD"=80, "EF"=100) # Paired samples in each group
N <- sum(ns)

## Generate X matrix
set.seed(11111)
m1.treat <- 100
m1.age <- 100
m1 <- m1.treat + m1.age
m0 <- m-m1
# Define DEGs
DEGs.treat <- seq(m1.treat)
DEGs.age <- seq((m1.treat+1),m1)
# Generate x1
Treats <- list(AB=c(rep(0,20), rep(1, 40)),
               BCD=c(rep(0, 50), rep(1, 30)),
               EF=c(rep(0, 30), rep(1, 70)))
treat <- Reduce("c", Treats)
# Generate x2, x3
age <- round(runif(N, min=18, max=80))
x3 <- rnorm(N)
# Obtain X matrix
X <- rbind(treat, age, x3)
colnames(X) <- paste0("s", 1:N)
p <- nrow(X)
xnames <- rownames(X)
# Re-organize X as a list
Xs <- list()
n.start <- 1
for (l in 1:length(ns)) {
  n.end <- n.start+ns[l]-1
  Xs[[names(ns)[l]]] <- X[,n.start:n.end]
  n.start <- n.end+1
}

## Generate betas
# Beta: U[1/2, 1] for treat; U[-.05, .05] for age; N(0, .5^2) for x3.
beta.treat <- rep(0,m)
beta.treat[DEGs.treat] <- runif(m1.treat, min=.5, max=1.5)
beta.age <- rep(0,m)
beta.age[DEGs.age] <- runif(m1.age, min=-0.1, max=0.1)
beta.x3 <- 1.0*rnorm(m)
Beta <- cbind(treat=beta.treat, age=beta.age, x3=beta.x3)

## Generate platform-dependent b0 and b1. b0 is different on different platforms
b0.all <- sapply(1:K.all, function(k) 8-.5*k+.5*rnorm(m))
# Define all b1 >0.1
b1.fixed <- c(.5, .6, .55, .65, .45, .7)
b1.all <- pmax(sweep(.05*matrix(rnorm(m*K.all), m), 2, b1.fixed, "+"), 0.1)
colnames(b0.all) <- colnames(b1.all) <- LETTERS[1:K.all]
# Group b0.all and b1.all into lists (easier to use)
b0s <- b1s <- list()
for (i in 1:length(PFs)) {
  b0s[[names(PFs)[i]]] <- b0.all[,PFs[[i]]]
  b1s[[names(PFs)[i]]] <- b1.all[,PFs[[i]]]
}

## Encapsulate data generation in a function to help run simulations with repetitions.
datagen <- function(X, Beta, A, sigma2s, sigmaU2s, b0.all, b1.all) {
  L <- ncol(A)
  N <- ncol(X)
  m <- length(sigma2s)
  K <- ncol(b0.all)
  Z <- matrix(rnorm(L*N),L)
  Upsilon <- matrix(sqrt(sigmaU2s)*rnorm(m*N), m)
  Y0 <- Beta %*% X + A %*% Z + Upsilon
  ## generate Y by b1*(Y0 +error)+b0 for all platforms.
  Emat <- sqrt(sigma2s)*matrix(rnorm(m*K*N),m*K)
  Y <- ((rep(1,K)%x%Y0)+Emat)*as.vector(b1.all) +as.vector(b0.all)
  rownames(Y) <- rep(names(sigma2s),K); colnames(Y) <- colnames(X)
  return(Y)
}

## Generate one set of the main data
set.seed(123)
Y.all <- datagen(X, Beta, A, sigma2s, sigmaU2s, b0.all, b1.all)
## Select subsets of Y to serve as datasets AB, BCD, and EF
slist.u <- cumsum(ns)
slist.l <- slist.u-ns+1
slist <- sapply(1:length(ns), function(k) seq(slist.l[k], slist.u[k]))
names(slist) <- DS
glist0 <- rep(PF.all, each=m)
glist <- lapply(DS, function(ds) which(glist0 %in% PFs[[ds]]))
names(glist) <- DS
Ys <- lapply(DS, function(ds) {
  Y.all[glist[[ds]], slist[[ds]]]
})
names(Ys) <- DS

## Apply 5-fold CV once to obtain the optimal Ls for each dataset
ss <- lapply(ns, cv.split, n.folds=5)
cvest.AB <- CV.InitEst(X=Xs$AB, Y=Ys$AB, K=2, folds=ss$AB, k.source=1, k.target=2)
cvest.EF <- CV.InitEst(X=Xs$EF, Y=Ys$EF, K=2, folds=ss$EF, k.source=1, k.target=2)

## BCD is more complicated
cvests.BCD <- list()
cvests.BCD[[1]] <- CV.InitEst(Xs$BCD, Ys$BCD, 3, ss$BCD, k.source=1, k.target=2)
cvests.BCD[[2]] <- CV.InitEst(Xs$BCD, Ys$BCD, 3, ss$BCD, k.source=1, k.target=3)
cvests.BCD[[3]] <- CV.InitEst(Xs$BCD, Ys$BCD, 3, ss$BCD, k.source=2, k.target=3)
CVMSEs2 <- sapply(cvests.BCD, function(ee) rowmeans(ee[["MSE"]]))
rownames(CVMSEs2) <- rownames(cvests.BCD[[1]]$MSE)
colnames(CVMSEs2) <- c("1to2", "1to3", "2to3")
best.Ls <- apply(CVMSEs2, 2, which.min)
best.Ls
Lstars <- c(AB=cvest.AB$Lstar, BCD=round(mean(best.Ls)), EF=cvest.EF$Lstar)

## Generate the covariate matrix just once for Analysis 4b.
# Group4: n=30. Group5: n=50.
n4.c <- 20
n4.t <- 10
n4 <- n4.c+n4.t
n5.c <- 15
n5.t <- 35
n5 <- n5.c+n5.t
n <- n4+n5

treat4 <- c(rep(0, n4.c), rep(1, n4.t))
treat5 <- c(rep(0, n5.c), rep(1, n5.t))

# Intentionally make age in group4 younger than group5
age4 <- round(runif(n4, 18, 60))
age5 <- round(runif(n5, 30, 80))

Xnew <- rbind(treat=c(treat4, treat5), age=c(age4,age5), x3=rnorm(n))
X4 <- Xnew[, 1:n4]
X5 <- Xnew[, (n4+1):n]

save.image("data/simulated_dataset.RData")




