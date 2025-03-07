## 05/25/2024. Main improvements in version 13: \Sigma_{Y} is now
## estimated by gppca(). See Ycircbar2Cov() for more details. Note
## that ModIntegrate() still needs more work.

## 05/12/2024. 1. ModIntegrate() works now -- it uses weighted average
## to estimate all parameters but U, which is estimated via SVD of the
## combined Ycircbar. 2. I tried to modify Ycircbar2Cov() so that the
## denominator in cov(R0sbar) is N-(p+1) instead of N-1, to account
## for the fact R0sbar are the centered residuals of the regression of
## X (with p regressors). However, the results, both in terms of
## parameter estimation and prediction, are slightly worse than using
## the simple denominator N-1. So I changed them back.

## 05/11/2024. Changed Ycirc2Cov() to Ycircbar2Cov() so that it is
## more efficient and reusable for the new ModIntegrate(). A side
## effect is that, as a workaround, I have to disable
## allow.b1.negative in InitEst() and InitEst.LargeSTD() -- that
## option needs more thinking anyway.

## 04/20/2024. Fixed some small issues in InitEst(), Prediction(), and
## CV.InitEst() in "auxilary.r", so that: (a) InitEst() is now robust
## for data with high collinearity, (b)Prediction(), and CV.InitEst()
## now works for LOOCV, in which the test data contains only one
## sample.

## 04/16/2024. Replace sweep(mat, 1, vec, op) by the simpler and faster
## (but more dangerous) direct operation between a matrix and a vector.

## 04/14/2024. All functions are working better and faster than the
## previous version in the new debugging/simulation script,
## debug_smallvar8.r. Only ModIntegrate() needs more tests.

## 04/13/2024. 1. Changed the interface of loglik() and gradLoglik()
## to make them more user-friendly: instead of a long sequence of
## parameters, they now accepts an encapsulated "params" object.

## 04/11/2024. 1. Fixed function gradLoglik().

## 04/10/2024. 1. 3x faster loglik().

## 04/09/2024. I noticed that the gradient descend algorithm, which
## includes GDfun(), loglik(), and gradLoglik(), does not work because
## it was written based on an old model in which sigma2 is a constant
## for all genes.  Changes: 1. Update loglik() so that it works with
## the vector version of sigma2s and fixed an error in
## D2inv. 2. Change Ytilde --> Ys to match Y^{(s)} in the report, and
## Ytildebar --> Ytilde to match the new definition of \tilde{Y},
## defined by eq:Ytilde in the report. The new \tilde{Y} is the
## platform-mean of Y^{(s)}.


## 04/08/2024. 1. Move many auxilary functions to "auxilary.r" to make
## this R script more readable. 2. replace the default svd() with
## fsvd(), a faster version of svd(). 3. Wrote a convenience function
## predict.mm() for making prediction from trained MM (MatchMixeR)
## models. 4. Wrote CV.InitEst() to facilitate the selection of L.

library(Rfast)
# source("auxilary.r")

#' Log-Likelihood Calculation for PPCA-XPN Model
#'
#' Computes the joint log-likelihood for gene expression data across multiple platforms using the PPCA-XPN model.
#'
#' Supports both paired and unpaired platforms, with optional covariates.
#'
#' This function produces log-likelihood without the "-2" multiplier.
#'
#' @param X Optional matrix of covariates (genes x covariates). Set to NULL if covariates are not used.
#' @param Y Gene expression matrix (genes x samples), with columns ordered by platform.
#' @param K Number of platforms.
#' @param params List of parameters containing:
#'   \describe{
#'     \item{b0}{Gene-specific intercepts (vector of length genes).}
#'     \item{b1}{Gene-specific scaling factors (vector of length genes).}
#'     \item{Beta}{Matrix of gene-specific regression coefficients (genes x covariates).}
#'     \item{U}{Matrix of latent factors (genes x L).}
#'     \item{dd}{Vector of latent factor variances (length L).}
#'     \item{sigmaU2s}{Gene-specific shared variance (length genes).}
#'     \item{sigma2s}{Gene-specific residual noise variance (length genes).}
#'   }
#' @param b1_min Minimum allowed value for gene-specific slopes to prevent numerical instability (default = 0.1).
#' @param verbose Logical; if TRUE, returns a list of detailed log-likelihood components instead of just total likelihood.
#' @return Total log-likelihood (numeric), or detailed components if `verbose = TRUE`.
#' @export
loglik <- function(X, Y, K, params, b1_min=.1, verbose=FALSE){
  b0 <- as.vector(params$b0)
  b1 <- as.vector(params$b1)
  # Prevent log(0) or divide by 0 error
  b1[b1 < b1_min] <- b1_min

  Beta <- params$Beta
  U <- params$U
  dd <- params$dd
  sigmaU2s <- params$sigmaU2s
  sigma2s <- params$sigma2s

  N <- ncol(Y)
  L <- length(dd)
  m <- nrow(U)

  if (m*K != nrow(Y)) stop("The rows of input matrix Y does not equal to m x K!")

  ## COMPUTATION STARTS HERE
  Const <- det(diag(x=1+1/dd,nrow=L) -t(U)%*%(U/(1+K*sigmaU2s/sigma2s)))/prod(1+1/dd)

  ## logdet: see eq:SigmaY-logdet
  logdet <- log(Const) +2*sum(log(b1)) +(K-1)*sum(log(sigma2s)) +sum(log(sigma2s+K*sigmaU2s)) +sum(log(1+dd))

  ## Ys: see eq:Ycirc
  Ys <- (Y-b0)/b1
  if (is.null(X) | is.null(Beta)){
    Ycirc <- Ys
  } else {
    Ycirc <- Ys-rep(1,K)%x%(Beta%*%X)
  }
  ## eq:Ycircbar-def
  Ycircbar <- ameans(array(Ycirc, c(m,K,N)),2)

  ## H, M, U2 (\tilde{U}), and W, are defined in between eq:GAAIinv and eq:SigmaY-inv
  hh <- sigma2s/(K*sigmaU2s+sigma2s)
  D2inv <- diag(x=1/dd, nrow=L) # this is a small matrix
  M <- solve(D2inv +diag(nrow=L) -t(U)%*%(U*hh)) ## M: eq:M-def
  U2 <- U*(sqrt(sigmaU2s)/(K*sigmaU2s+sigma2s))
  ## W <- diag(sigmaU2s/(sigma2s*(K*sigmaU2s+sigma2s))) +U2%*%(M%*%t(U2))

  ## Three terms in -2*loglik, defined in eq:joint-likelihood2
  Term1 <- m*N*K*log(2*pi) +N*logdet
  Term2 <- sum(rowsums(Ycirc^2)/rep(sigma2s, K))
  Ycircbar2 <- t(U2)%*%Ycircbar
  Term3a <- K^2*sum((Ycircbar*(sigmaU2s/(sigma2s*(K*sigmaU2s+sigma2s))))*Ycircbar)
  Term3b <- K^2*tr(M%*%tcrossprod(Ycircbar2))
  Term3 <- Term3a+Term3b

  ## loglik
  loglik <- -(Term1+Term2-Term3)/2
  if (verbose) { # for debugging
    return(c(logdet=logdet, Term1=Term1, Term2=Term2, Term3a=Term3a, Term3b=Term3b, Term3=Term3, loglik=loglik))
  } else {
    return(loglik)
  }
}

#' Gradient of Log-Likelihood with Respect to b1
#'
#' Computes the gradient of the log-likelihood with respect to the gene-specific scaling factors (b1).
#'
#' All other parameters can be estimated by closed-form conditional MLE (`b0` and betas)
#' or moment estimator (covariance related parameters, `U`, `dd`, `sigma2s`, and `sigmaU2s`).
#' If there is no covariate in the data, set `X` and `Beta` to NULL.
#'
#' @param X Optional matrix of covariates.
#' @param Y Gene expression matrix.
#' @param K Number of platforms.
#' @param params Parameter list (see \code{loglik} for details).
#' @param b1_min Minimum allowed b1 value (default = 0.1).
#' @return Vector of length equal to number of genes, giving gradient for each gene.
#' @export
gradLoglik <- function(X, Y, K, params, b1_min=.1){
  b0 <- as.vector(params$b0)
  b1 <- as.vector(params$b1)
  Beta <- params$Beta
  U <- params$U
  dd <- params$dd
  sigmaU2s <- params$sigmaU2s
  sigma2s <- params$sigma2s

  N <- ncol(Y)
  L <- length(dd)
  m <- length(b1)/K
  ## Ybar <- rowMeans(Y)
  ## to prevent log(0) or divide by 0 error
  b1[b1 < b1_min] <- b1_min
  Ys <- (Y-b0)/b1
  ## Ytilde <- ameans(array(Ys, c(m,K,N)), 2)
  if (is.null(X) | is.null(Beta)){
    Ycirc <- Ys
  } else {
    Ycirc <- Ys-rep(1,K)%x%(Beta%*%X)
  }
  Ycircbar <- ameans(array(Ycirc, c(m,K,N)), 2)
  ## lambdas <- dd*K; mus <- lambdas/(sigma2+lambdas)
  ## compute M, W, and F
  hh <- sigma2s/(K*sigmaU2s+sigma2s)
  D2inv <- diag(x=1/dd, nrow=L) # this is a small matrix
  ## M: eq:M-def
  M <- solve(D2inv +diag(nrow=L) -t(U)%*%(U*hh))
  U2 <- U*(sqrt(sigmaU2s)/(K*sigmaU2s+sigma2s))
  F <- Ycircbar*(sigmaU2s/(sigma2s*(K*sigmaU2s+sigma2s))) +(U2%*%M)%*%(t(U2)%*%Ycircbar)
  ## only produce the gradients for b1 based on eq:grad-b1. .
  Second <- Ycirc/rep(sigma2s,K) -K*(rep(1,K)%x%F)
  grad.b1 <- -N/b1 +rowsums((Y-b0)*Second)/b1^2
  return(grad.b1)
}

#' Estimate Regression Coefficients (Beta)
#'
#' Estimates gene-specific regression coefficients using expression data and gene-wise covariate slopes (b1).
#'
#' @param X Matrix of covariates.
#' @param Y Gene expression matrix.
#' @param b1 Gene-specific scaling factors.
#' @param K Number of platforms.
#' @param b1_min Minimum allowed b1 value (default = 0.1).
#' @return Matrix of gene-specific regression coefficients (genes x covariates).
BetaEst <- function(X, Y, b1, K, b1_min=.1){
  if (is.null(X)) {
    return(NULL)
  } else {
    b1 <- as.vector(b1)
    m <- length(b1)/K
    N <- ncol(Y)
    b0 <- rowmeans(Y)
    ## to prevent divide by 0 error
    b1[b1 < b1_min] <- b1_min
    Ys <- (Y-b0)/b1
    Ycircbar <- ameans(array(Ys, c(m,K,N)), 2)
    ## XXX <- t(X)%*%solve(tcrossprod(X))
    XXX <- t(X)%*%rsolve2(t(X))
    return(Ycircbar %*% XXX)
  }
}

#' Gradient Descent Optimization for PPCA-XPN Model
#'
#' Uses gradient descent to optimize gene-specific scaling factors (b1) and regression coefficients (Beta) for PPCA-XPN.
#'
#' `b0` is always `rowMeans(Y)`. All other parameters (`U`, `dd`, `sigmaU2s`,
#' `sigma2s`) need to be computed via `InitEst()`.
#'
#' @param X Optional matrix of covariates.
#' @param Y Gene expression matrix.
#' @param K Number of platforms.
#' @param params Initial parameter estimates (see \code{loglik} for details).
#' @param platforms Optional vector of platform names (length K).
#' @param s0 Initial step size for gradient descent.
#' @param smax Maximum allowed step size.
#' @param sigma2.min Minimum residual variance.
#' @param sigma2.max Maximum residual variance.
#' @param b1_min Minimum b1 value.
#' @param b1.max Maximum b1 value.
#' @param beta.min Minimum Beta coefficient value.
#' @param beta.max Maximum Beta coefficient value.
#' @param max.iter Maximum number of iterations.
#' @param tol Convergence tolerance.
#' @param verbose Logical; if TRUE, returns optimization trace.
#' @return List of optimized parameters (same structure as initial params).
#' @export
GDfun <- function(X = NULL, Y, K, params, platforms=NULL,
                  s0=0.01, smax=0.5,
                  sigma2.min=0.01, sigma2.max=2,
                  b1_min=0.1, b1.max=5, beta.min = -5, beta.max = 5,
                  max.iter=30, tol=0.1, verbose=FALSE){
  if((!is.null(X)) & (!is.matrix(X))) { # If there is only one covariate and X is not a matrix.
    X <- matrix(X, nrow = 1)
  }

  b1 <- params$b1
  U <- params$U
  dd <- params$dd
  sigmaU2s <- params$sigmaU2s
  sigma2s <- params$sigma2s
  N <- ncol(Y)
  m <- nrow(Y)/K
  L <- length(dd)

  ## If there are no U and dd (means that L=0), we do not do any
  ## gradient descend, and just return the params.
  if (is.null(U)) {
    if (verbose) {
      return(list(params=params))
    } else{
      return(params)
    }
  }

  ## use genenames names of covariates if available
  gnames <- rownames(Y)[1:m]
  ## assign/compute the initial parameters
  b0 <- rowmeans(Y)
  b1.old <- as.vector(b1)
  if (is.null(X)) {
    Beta.old <- Xnames <- NULL; p <- 0
  } else {
    Xnames <- rownames(X); p <- nrow(X)
    Beta.old <- BetaEst(X, Y, b1.old, K)
  }
  ## as of 04/12/2024, theta only contains b1. We may be able to
  ## expand it to include sigmaU2s and sigma2s in the future.
  theta.old <- b1.old
  ## compute the initial gradient
  params.old <- list(Beta=Beta.old, b0=b0, b1=b1.old, U=U, dd=dd,
                     sigmaU2s=sigmaU2s, sigma2s=sigma2s)
  grad.old <- gradLoglik(X,Y,K, params.old, b1_min=b1_min)
  ## grad.old <- c(g0$grad.sigma2, g0$grad.b1)
  Norm.old <- sqrt(sum(grad.old^2))
  v.old <- grad.old/Norm.old
  ll.old <- loglik(X, Y, K, params.old, b1_min=b1_min)
  ## the main loop
  j <- 1
  s.j <- s0
  loglik.old <- 10^5
  diff.j <- 1e6
  history <- t(c(loglik=ll.old, s=s0, r=1, theta.old))
  while (j < max.iter && diff.j>tol){
    ## compute the new theta/parameters
    theta.j <- theta.old + s.j*v.old
    ## threshold the estimated b1
    b1.j <- pmax(pmin(theta.j, b1.max), b1_min)
    if (is.null(X)) {
      Beta.j <- NULL
    } else {
      Beta.j <- pmax(pmin(BetaEst(X, Y, b1.j, K), beta.max), beta.min)
    }
    ## calculate new log-likelihood
    params.j <- list(Beta=Beta.j, b0=b0, b1=b1.j, U=U, dd=dd,
                     sigmaU2s=sigmaU2s, sigma2s=sigma2s)
    ll.j <- loglik(X, Y, K, params.j)
    ## If log-likelihood decreases (it's supposed to increase),
    ## keep cutting step size s.j in half
    k = 1
    while(ll.j < ll.old && k < max.iter){
      s.j = s.j/2
      ## compute the new theta/parameters
      theta.j <- theta.old + s.j*v.old
      ## threshold the estimated b1
      b1.j <- pmax(pmin(theta.j, b1.max), b1_min)
      ## calculate new log-likelihood
      params.j <- list(Beta=Beta.j, b0=b0, b1=b1.j, U=U, dd=dd,
                       sigmaU2s=sigmaU2s, sigma2s=sigma2s)
      ll.j <- loglik(X, Y, K, params.j)
      k = k+1
    }
    ## compute the new gradients
    params.j <- list(Beta=Beta.j, b0=b0, b1=b1.j, U=U, dd=dd,
                     sigmaU2s=sigmaU2s, sigma2s=sigma2s)
    grad.j <- gradLoglik(X,Y,K, params.j, b1_min=b1_min)
    Norm.j <- sqrt(sum(grad.j^2))
    v.j <- grad.j/Norm.j
    ## update the optimal step for the next iteration
    vv <- sum(v.j*v.old)
    r.j <- Norm.j*vv/(Norm.old*vv -Norm.j)
    s.j <- r.j*s.j
    ## to be on the safe side
    if (s.j < 0 | s.j > smax) s.j <- smax
    ## check the difference for convergence
    diff.j <- abs(ll.j -ll.old)
    ## final steps
    history <- rbind(history, c(ll.j, s.j, r=r.j, theta.j))
    j=j+1
    ll.old <- ll.j
    theta.old=theta.j
    grad.old=grad.j
    Norm.old <- Norm.j
    v.old <- v.j
  }
  b1 <- b1.j
  ## finishing touches
  b0 <- matrix(b0, nrow=m)
  b1 <- matrix(b1, nrow=m)
  if (!is.null(gnames)) {
    rownames(b0) <- rownames(b1) <- rownames(U) <- gnames
  }
  if (!is.null(Beta.j)) {
    if (!is.null(gnames)) rownames(Beta.j) <- gnames
    if (!is.null(Xnames)) colnames(Beta.j) <- Xnames
  }
  if (!is.null(platforms)) {
    if (length(platforms) != K) {
      warning("Length of platforms must be K.")
    } else {
      colnames(b0) <- colnames(b1) <- platforms
    }
  }
  ## the only updated parameters are Beta and b1.
  params <- list(Beta=Beta.j, b0=b0, b1=b1, U=U, dd=dd, sigmaU2s=sigmaU2s, sigma2s=sigma2s, N=N, K=K)
  ##
  if (verbose) {
    Ys <- (Y-as.vector(b0))/as.vector(b1)
    if (is.null(X)){
      Ycirc <- Ys
    } else {
      Ycirc <- Ys-rep(1,K)%x%(Beta.j%*%X)
    }
    Ycircbar <- ameans(array(Ycirc, c(m,K,N)), 2)
    return(list(params=params, m=m, p=p, Ycirc=Ycirc,
                Ycircbar=Ycircbar, loglik=ll.j, loglik.old=ll.old,
                iters=j-1, history=history))
  } else {
    return(params)
  }
}

#' Cross-Platform Gene Expression Prediction
#'
#' Predicts gene expression on a target platform using a trained PPCA-XPN model.
#'
#' The \code{trained.model} should contain:
#' \describe{
#'   \item{b0}{Gene-specific intercepts (length = genes).}
#'   \item{b1}{Gene-specific scaling factors (length = genes).}
#'   \item{Beta}{Regression coefficients w.r.t \code{X} (genes x covariates).}
#'   \item{U}{Latent factors (genes x L).}
#'   \item{dd}{Latent factor variances (length L).}
#'   \item{sigmaU2s}{Gene-specific shared variance (length = genes).}
#'   \item{sigma2s}{Gene-specific residual variances (length = genes).}
#' }
#'
#' @param Ysource Expression matrix for source platform (genes x samples).
#' @param X Optional covariate matrix for target samples.
#' @param trained.model Trained PPCA-XPN model (see details).
#' @param k.source Index of source platform.
#' @param k.target Index of target platform.
#' @param b1_min Minimum allowed b1 value.
#' @param min.sigmaU2 Minimum allowed shared variance.
#' @param verbose Logical; if TRUE, returns intermediate prediction components.
#' @return Predicted gene expression matrix for target platform (genes x samples).
#' @export
Prediction <- function(Ysource, X, trained.model, k.source=1, k.target=2, b1_min=0.02, min.sigmaU2=0.05, verbose=FALSE){

  if((!is.null(X)) & (!is.matrix(X))) {
    X <- matrix(X, nrow = 1)
  }
  ## collect information from trained.model
  Beta <- trained.model$Beta
  b0k <- trained.model$b0[, k.source]
  b0k2 <- trained.model$b0[, k.target]
  b1k <- trained.model$b1[, k.source]
  b1k2 <- trained.model$b1[, k.target]
  U <- trained.model$U; dd <- trained.model$dd
  sigmaU2s <- pmax(trained.model$sigmaU2s, min.sigmaU2)
  sigma2s <- trained.model$sigma2s
  ## stop if NAs are detected
  if (sum(is.na(Ysource))>0) stop("Prediction() does not work with NAs and is not robust to outliers. Consider running winsor() on the input data to impute missing values and remove outliers.")
  ## 1. Compute the "standardized" Yk, which is B_{1,k}^{-1} (Y_k -mu_k)
  Yk <- as.matrix(Ysource)
  m <- nrow(Yk); N <- ncol(Yk); L <- length(dd)
  if (is.null(X) | is.null(Beta)){
    BetaX <- matrix(0, m, N)
  } else {
    BetaX <- Beta%*%as.matrix(X)
  }
  ## use a simple method for those with very small b1k
  ids0 <- which(abs(b1k)<b1_min)
  ## ## to prevent divide by 0 error and keep the +/- sign
  ## b1k <- pmax(abs(b1k), b1_min)*sign(b1k)
  muk <- BetaX*b1k+b0k
  muk2 <- BetaX*b1k2+b0k2
  Ycirc.k <- (Yk-muk)/b1k
  Ycirc.k[ids0, ] <- 0
  ## Wu is a vector
  Wu <- sigmaU2s/(sigmaU2s+sigma2s)
  ## GSI is gene-specific information
  GSI <- Ycirc.k*(b1k2*Wu)
  ## SI is the shared information
  if (is.null(dd)){ #L=0
    SI <- matrix(0, m, N)
  } else {
    Term1 <- b1k2*sqrt(sigmaU2s) #a vector
    PCs <- t(U)%*% (Ycirc.k*(Wu/sqrt(sigmaU2s)))
    UD2 <- sweep(U, 2, dd, "*")
    UsqrtWu <- U*sqrt(Wu); UWuU <- Crossprod(UsqrtWu, UsqrtWu)
    Term2 <- (U%*%solve(diag(1/dd, nrow=L) +UWuU))*Wu
    SI <- (UD2%*%((diag(nrow=L) -t(U)%*%Term2)%*%PCs) -Term2%*%PCs)*Term1
  }
  ## The estimated gene expressions on platform k2
  Yhat <- muk2 +GSI +SI
  ## Now just replace those genes with very small b1k by muk2
  Yhat[ids0,] <- muk2[ids0]%*%t(rep(1,N))
  ## for debugging purpose
  if (verbose) {
    return(list(Yhat=Yhat, muk2=muk2, GSI=GSI, SI=SI))
  } else {
    return(Yhat)
  }
}


#' Initialize Parameters for PPCA-XPN
#'
#' Initializes gene-specific intercepts, slopes, covariances, and latent factors for PPCA-XPN.
#'
#' @param X Optional covariate matrix.
#' @param Y Gene expression matrix.
#' @param K Number of platforms.
#' @param L Number of latent factors.
#' @param min.sigma2 Minimum residual variance.
#' @param min.sigmaU2 Minimum shared variance.
#' @return List of initialized parameters (see \code{loglik} for structure).
#' @export
InitEst <- function(X, Y, K, L, min.sigma2=0.01, min.sigmaU2=0.05){
  if((!is.null(X)) & (!is.matrix(X))) {
    X <- matrix(X, nrow = 1)
  }
  N <- ncol(Y)
  m <- nrow(Y)/K
  ## use gene names of Y if available
  gnames <- rownames(Y)[1:m]
  ## remove genes with very small STD
  varY <- matrix(rowVars(Y), ncol=K)
  ids1 <- rowMins(varY)>=min.sigma2 #genes with large STD
  low.var.genes <- which(!ids1)
  ## the main part
  Y2 <- Y[rep(ids1, K),]
  initEst1 <- InitEst.LargeSTD(X, Y2, K, L, min.sigma2=min.sigma2, min.sigmaU2=min.sigmaU2)
  ## initEst1 <- InitEst.LargeSTD.old(X, Y2, K, L, min.sigmaU2=min.sigmaU2)

  #Combine the simple estimators for genes with very small STD and all other genes
  b0 <- matrix(rowmeans(Y), ncol=K); b1 <- matrix(0, m, K)
  rownames(b0) <- rownames(b1) <- gnames
  colnames(b0) <- colnames(b1) <- paste0("Platform", 1:K)
  b1[ids1,] <- initEst1$b1
  ## estimate beta matrix
  if (is.null(X)) {
    Beta <- NULL
  } else {
    Beta <- matrix(0, m, nrow(X))
    rownames(Beta) <- gnames; colnames(Beta) <- rownames(X)
  }
  Beta[ids1,] <- initEst1[["Beta"]]
  ## other parameters. Note that due to collinearity, U produced by
  ## InitEst.LargeSTD() may have less than L columns
  U0 <- initEst1$U
  if (is.null(U0)) { #this could happen when L==0
    U <- NULL
  } else {
    U <- matrix(0, m, ncol(U0))
    rownames(U) <- gnames; colnames(U) <- colnames(U0)
    U[ids1,] <- initEst1$U
  }
  dd <- initEst1$dd
  varprops <- initEst1$varprops
  sigmaU2s <- rep(1, m); names(sigmaU2s) <- gnames
  sigmaU2s[ids1] <- initEst1$sigmaU2s
  sigma2s <- rep(0, m); names(sigma2s) <- gnames
  sigma2s[ids1] <- initEst1$sigma2s
  return(list(b0=b0, b1=b1, Beta=Beta, U=U, dd=dd, varprops=varprops,
              sigmaU2s=sigmaU2s, sigma2s=sigma2s, N=N, K=K, low.var.genes=low.var.genes))
}

#' Estimate Covariance from Platform-Averaged Expression (Ycircbar)
#'
#' Calculates gene-wise covariance matrix components based on platform-averaged expression residuals.
#'
#' @param Ycircbar Gene expression residual matrix averaged across platforms.
#' @param K Number of platforms.
#' @param L Number of latent factors.
#' @param min.sigma2 Minimum residual variance.
#' @param min.sigmaU2 Minimum shared variance.
#' @param W Weight parameter to control the variance split.
#' @return List with U, dd, sigmaU2s, sigma2s.
#' @export
Ycircbar2Cov <- function(Ycircbar, K, L, min.sigma2=0.01, min.sigmaU2=0.05, W=0.5) {
  N <- ncol(Ycircbar); m <- nrow(Ycircbar)
  if (K==1) { # no pairing
    ## scale=TRUE --> \Sigma_{Y} is a correlation matrix
    gg <- gppca(t(Ycircbar), scale=TRUE, retx=FALSE, L=L, min.sigma2=min.sigma2+min.sigmaU2)
    U <- gg$U
    dd <- gg$dd
    sigmaU2s <- gg$sigma2s*W
    sigma2s <- gg$sigma2s*(1-W)
  } else { # with K>=2 paired platforms
    ## scale=FALSE --> \Sigma_{Y} is NOT a correlation matrix
    gg <- gppca(t(Ycircbar), scale=FALSE, retx=FALSE, L=L, min.sigma2=min.sigma2+min.sigmaU2)
    U <- gg$U
    dd <- gg$dd
    if (is.null(U)) { #no AA'
      varR0bar <- gg$sigma2s
    } else {
      diagAA <- drop((U*U)%*%dd)
      varR0bar <- diagAA+gg$sigma2s
    }
    sigma2s <- pmax(K*(1-varR0bar)/(K-1), min.sigma2)
    sigmaU2s <- pmax(gg$sigma2s-sigma2s/K, min.sigmaU2)
  }
  ## prop. of variance
  if (is.null(dd)) {
    varprops <- NULL
  } else {
    varprops <- cumsum(dd)/sum(dd)
  }
  return(list(U=U, dd=dd, varprops=varprops, sigmaU2s=sigmaU2s, sigma2s=sigma2s))
}

## new implementation as of 04/14/2024. This function only works for
## those genes with STD>>0. The new InitEst() is a wrapper that
## includes Winsorizing and prediction.
InitEst.LargeSTD <- function(X, Y, K, L, min.sigma2=0.01, min.sigmaU2=0.05, W=0.5){
  N <- ncol(Y)
  m <- nrow(Y)/K
  ## use genenames names of Y if available
  gnames <- rownames(Y)[1:m]
  ## b0hat is trivial
  b0hat <- rowmeans(Y)
  ## center both X and Y
  R0 <- Y-b0hat
  if (!is.null(X)) {
    if((!is.matrix(X))) {
      X <- matrix(X, nrow = 1)
    }
    p <- nrow(X); Xnames <- rownames(X)
    Xc <- X-rowmeans(X)
    ## regress out X
    HatMat <- rhat(t(Xc))
    R0 <- R0%*%(diag(N)-HatMat)
    ## After this step, R0 is the set of residuals, which is denoted as
    ## R^{(0)} in my notes. Again, no new object is created to save
    ## memory
  } else {
    p <- 0
    Xnames <- NULL
  }
  ## the initial estimate of b1 is simply sample STDs
  s0 <- sqrt(rowsums(R0*R0)/(N-1))
  b1 <- matrix(s0, nrow=m)
  colnames(b1) <- paste0("Platform", 1:K)
  ## standardize R0 --> R0s (an approximation of Ycirc)
  R0s <- R0/s0
  R0sbar <- ameans(array(R0s, c(m,K,N)), 2)
  CovEsts <- Ycircbar2Cov(R0sbar, K, L, min.sigma2=min.sigma2, min.sigmaU2=min.sigmaU2, W=W)
  U <- CovEsts$U
  ## Disabled "allow.b1.negative" on 05/11/2024
  ## if (allow.b1.negative) {
  ##     ## flip the sign of b1 according to sumcor
  ##     neg.ids <- which(CovEsts$sumcor<0, arr.ind=TRUE)
  ##     ## the first platform always have b1>0 so don't toggle its sign
  ##     neg.ids <- neg.ids[neg.ids[,2] !=1, ]
  ##     b1[neg.ids] <- -b1[neg.ids]
  ## }
  ## Betahat
  Betahat <- BetaEst(X, Y, b1, K)
  ## Return b0 and b1 in matrix format
  b0 <- matrix(b0hat, nrow=m); colnames(b0) <- paste0("Platform", 1:K)
  if (!is.null(gnames)) {
    rownames(b0) <- rownames(b1) <- gnames
    if (!is.null(U)) {
      rownames(U) <- gnames; colnames(U) <- paste0("PC", 1:ncol(U))
    }
  }
  if (!is.null(Betahat)){
    if (!is.null(gnames)) rownames(Betahat) <- gnames
    if (!is.null(Xnames)) colnames(Betahat) <- Xnames
  }
  return(list(b0=b0, b1=b1, Beta=Betahat, U=U, dd=CovEsts$dd,
              varprops=CovEsts$varprops, sigmaU2s=CovEsts$sigmaU2s,
              sigma2s=CovEsts$sigma2s))
}

#' Integration of Multiple PPCA-XPN Models
#'
#' Combines models estimated on different platform pairs (returned by `GDfun()`)
#' into a single unified model.
#'
#' @param modList List of trained PPCA-XPN models.
#' @param Xlist Optional list of covariate matrices (one per platform).
#' @param Ylist List of gene expression matrices (one per platform).
#' @param PairList Matrix (pairs x 2) specifying platform pairs included in each model.
#' @return Combined integrated PPCA-XPN model.
#' @export
ModIntegrate <- function(modList, Xlist=NULL, Ylist, PairList) {
  ## collect information from the inputs
  ## use genenames names of covariates if available
  gnames <- rownames(modList[[1]][["b0"]])
  G <- length(modList)
  platforms <- Reduce("union", PairList)
  nPFs <- length(platforms)
  m <- length(gnames)
  ## the effective sample size for each group is Ng*Kg
  Ns0 <- sapply(modList, function(mod) mod[["N"]])
  Ks <- sapply(modList, function(mod) mod[["K"]])
  ## Ns are the effective samples in each platform
  Ns <- Ns0*Ks
  ## extract the estimates from each model
  sigmaU2mat <- sapply(modList, function(mod) mod[["sigmaU2s"]])
  sigma2mat <- sapply(modList, function(mod) mod[["sigma2s"]])
  Betas <- lapply(modList, function(mod) mod[["Beta"]])
  b0s <- lapply(modList, function(mod) mod[["b0"]])
  b1s <- lapply(modList, function(mod) mod[["b1"]])
  Us <- lapply(modList, function(mod) mod[["U"]])
  dds <- lapply(modList, function(mod) mod[["dd"]])
  ## sanity check
  Ls <- sapply(dds, length); L <- min(Ls)
  if (diff(range(Ls))!=0) {
    warning(paste0("Number of PCs (L) are not equal across all groups. The integrated model will use the minimum number of PCs, Lmin=", L))
  }
  ## use weighted mean for sigma2s, sigmaU2s, and Beta
  ww <- Ns/sum(Ns)
  sigma2s <- drop(sigma2mat%*%ww)
  sigmaU2s <- drop(sigmaU2mat%*%ww)
  ## dd
  ddmat <- sapply(dds, function(ds) ds[1:L])
  dd <- drop(ddmat%*%ww)
  ## Beta
  if (any(sapply(Betas, is.null))) {
    Betahat <- NULL
  } else {
    Betahat <- Reduce("+", lapply(1:G, function(g) ww[g]*Betas[[g]]))
  }
  ## b0 and b1 needs a bit workx
  b0 <- b1 <- matrix(0, m, nPFs)
  colnames(b0) <- colnames(b1) <- platforms
  for (plfm in platforms){
    ids <- lapply(PairList, function(x) which(x==plfm))
    ## sanity check
    LL <- sapply(ids, length)
    if (max(LL>1)) stop("Please check PairList; one platform can only appear once in each group.")
    ## now the main algorithm
    Grps <- which(LL==1)
    if (length(Grps)==1){ #no need to do weighted average
      WithinGrpID <- ids[[Grps]]
      b0[, plfm] <- b0s[[Grps]][((WithinGrpID-1)*m+1):(WithinGrpID*m)]
      b1[, plfm] <- b1s[[Grps]][((WithinGrpID-1)*m+1):(WithinGrpID*m)]
    } else {
      b0[, plfm] <- Reduce("+", lapply(Grps, function(g) {
        WithinGrpID <- ids[[g]]
        return(Ns[g]*b0s[[g]][((WithinGrpID-1)*m+1):(WithinGrpID*m)])
      }))/sum(Ns[Grps])
      b1[, plfm] <- Reduce("+", lapply(Grps, function(g) {
        WithinGrpID <- ids[[g]]
        return(Ns[g]*b1s[[g]][((WithinGrpID-1)*m+1):(WithinGrpID*m)])
      }))/sum(Ns[Grps])
    }
  }
  ## Estimate U, modified from Ycircbar2Cov()
  if (L==0) {                         #do not compute U and dd
    U <- dd <- varprops <- NULL
  } else {
    ## use the combined Ycircbars to estimate U and dd
    Ycircbars <- Reduce("cbind", lapply(1:G, function(g) {
      Yg <- Ylist[[g]]
      platforms.g <- PairList[[g]]; Kg <- length(platforms.g)
      b0g <- as.vector(b0[, platforms.g])
      b1g <- as.vector(b1[, platforms.g])
      Ys <- (Yg-b0g)/b1g
      Ytilde <- ameans(array(Ys, c(m,Kg,ncol(Ys))), 2)
      if (is.null(Xlist)){
        Ycirc <- Ys
      } else {
        Xg <- Xlist[[g]]
        if (is.null(Xg)) {
          Ycirc <- Ys
        } else {
          BetaX <- Ytilde %*% rhat(t(Xg))
          Ycirc <- Ys -rep(1,Kg)%x%BetaX
        }
      }
      Ycircbar <- ameans(array(Ycirc, c(m,Kg,ncol(Ycirc))), 2)
    }))
    ## Using SVD of Ycircbar to estimate U
    ss <- fsvd(Ycircbars, k=L)
    U <- ss$u
    ## varprop is very rough
    varprops <- cumsum(ss$d^2)/sum(Ycircbars^2)
    ## there is a chance that ncol(U) < L, because some singular
    ## values may be numerically zero (especially when L is large)
    if (ncol(U) < L) {
      warning(paste0("Effective rank of Ycircbar is less than Lmin=", L, ", so we reduce the number of PCs to L=", ncol(U), " to ensure the numerical stability of the results."))
      L <- ncol(U)
    }
  }
  ## Final step: add row/column names to the output. Note that both
  ## b0 and b1 are matrices, not long vectors
  rownames(b0) <- rownames(b1) <- rownames(U)  <- names(sigmaU2s)  <- names(sigma2s) <- gnames
  if (!is.null(Betahat)) {
    rownames(Betahat) <- gnames
    colnames(Betahat) <- rownames(Xs[[1]])
  }
  colnames(U) <- paste0("PC", 1:ncol(U))
  return(list(b0=b0, b1=b1, Beta=Betahat, U=U, dd=dd, varprops=varprops,
              sigma2s=sigma2s, sigmaU2s=sigmaU2s, Ns=Ns0, Ks=Ks))
}


## We may consider writing a wrapper for a dummy user that combines
## InitEst(), GDfun(), and Prediction().
## EigenXPNLearn <- function(Xtrain, Ytrain, K, L, ...) {
##   init.est <- InitEst(X, Y, K, L)
##   final.est <- GDfun(X, Y, K, U, dd, sigma2.0, b1, ...)
##   return(...)
## }



## This function takes the output of EigenXPNLearn, a new set of data,
## and produces the predicted (normalized) gene expressions.
## EigenXPNPred <- function(EigenXPNModel, Xtest, Ytest, k1, k2) { ...
## }




#' Train EigenXPN Model
#'
#' A high-level wrapper for learning the PPCA-XPN model.
#'
#' This wrapper:
#' 1. Initializes the parameters using `InitEst()`.
#' 2. Optimizes the model parameters using `GDfun()`.
#'
#' @param Xtrain Optional covariate matrix (genes x covariates) for training data.
#' @param Ytrain Gene expression matrix (genes x samples) for training data.
#' @param K Number of platforms.
#' @param L Number of latent factors.
#' @param ... Additional parameters passed to `GDfun()`.
#' @return List representing the trained EigenXPN model, containing:
#' \describe{
#'   \item{b0}{Gene-specific intercepts.}
#'   \item{b1}{Gene-specific scaling factors.}
#'   \item{Beta}{Gene-specific regression coefficients (if X provided).}
#'   \item{U}{Latent factors.}
#'   \item{dd}{Latent factor variances.}
#'   \item{sigmaU2s}{Shared gene-specific variance.}
#'   \item{sigma2s}{Gene-specific noise variance.}
#' }
#' @export
#' @examples
#' Ytrain <- matrix(rnorm(1000), nrow=100, ncol=10)  # 100 genes, 10 samples
#' model <- EigenXPNLearn(NULL, Ytrain, K=2, L=3)
EigenXPNLearn <- function(Xtrain, Ytrain, K, L, ...) {

  # Step 1: Initialize model parameters.
  init.est <- InitEst(Xtrain, Ytrain, K, L)

  # Step 2: Optimize parameters via gradient descent.
  final.est <- GDfun(Xtrain, Ytrain, K, init.est, ...)

  # Return the final trained model (same structure as params in loglik())
  return(final.est)
}

#' Predict Gene Expression Using Trained EigenXPN Model
#'
#' Predicts gene expression on a new platform using a previously trained EigenXPN model.
#'
#' This wrapper simplifies calling `Prediction()` using the trained model.
#'
#' @param EigenXPNModel Trained model from \code{EigenXPNLearn()}.
#' @param Xtest Optional covariate matrix (genes x covariates) for new data.
#' @param Ytest Gene expression matrix (genes x samples) for new data.
#' @param k1 Index of the source platform.
#' @param k2 Index of the target platform.
#' @param ... Additional arguments passed to \code{Prediction()}.
#' @return Predicted gene expression matrix (genes x samples) for platform k2.
#' @export
#' @examples
#' Ytrain <- matrix(rnorm(1000), nrow=100, ncol=10)  # 100 genes, 10 samples
#' model <- EigenXPNLearn(NULL, Ytrain, K=2, L=3)
#' Ytest <- Ytrain  # Example reuse for simplicity
#' pred <- EigenXPNPred(model, NULL, Ytest, k1=1, k2=2)
EigenXPNPred <- function(EigenXPNModel, Xtest, Ytest, k1, k2, ...) {
  # Slice out source platform samples
  genes <- nrow(Ytest)
  samples <- ncol(Ytest)
  samples.per.platform <- samples / EigenXPNModel$K

  if (samples %% EigenXPNModel$K != 0) {
    stop("Ytest columns must evenly divide across K platforms.")
  }

  # Extract source platform data
  source.start <- (k1 - 1) * samples.per.platform + 1
  source.end <- k1 * samples.per.platform
  Ysource <- Ytest[, source.start:source.end, drop=FALSE]

  # Use Prediction() to predict target platform expression
  predicted <- Prediction(
    Ysource = Ysource,
    X = Xtest,
    trained.model = EigenXPNModel,
    k.source = k1,
    k.target = k2,
    ...
  )

  return(predicted)
}
