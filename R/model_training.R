#' Log-Likelihood Computation for Latent Gene Expression Models
#'
#' Computes the log-likelihood (without the \code{-2} multiplier) for a latent variable model based on
#' observed gene expression data, structured into \code{K} domains or platforms. The model includes
#' affine transformations, latent correlation structure via eigen decomposition, and covariate adjustment if applicable.
#'
#' @param X Optional covariate matrix. If there are no covariates, both \code{X} and \code{Beta} in \code{params} should be set to \code{NULL}.
#' @param Y A numeric matrix of observed data with \code{m*K} rows and \code{N} columns, where \code{m} is the number of genes and \code{K} is the number of platforms.
#' @param K Integer. Number of platforms (or blocks) the rows of \code{Y} are divided into.
#' @param params A list containing the model parameters:
#' \itemize{
#'   \item \code{b0}: Vector of intercepts (length \code{m}).
#'   \item \code{b1}: Vector of slopes (length \code{m}).
#'   \item \code{Beta}: Coefficient matrix for covariates (optional).
#'   \item \code{U}: Eigenvector matrix of the latent correlation structure.
#'   \item \code{dd}: Corresponding eigenvalues (\code{dd[l]} equals squared singular values).
#'   \item \code{sigmaU2s}: Platform-specific variance components.
#'   \item \code{sigma2s}: Residual variances per gene.
#' }
#' @param b1_min Numeric. Minimum allowed value for \code{b1} slopes to avoid numerical instability. Defaults to \code{0.1}.
#' @param verbose Logical. If \code{TRUE}, returns detailed components of the likelihood calculation for debugging. Defaults to \code{FALSE}.
#'
#' @return
#' If \code{verbose = FALSE}, returns a numeric scalar representing the log-likelihood.
#' If \code{verbose = TRUE}, returns a named vector containing intermediate quantities (\code{logdet}, \code{Term1}, \code{Term2}, \code{Term3a}, \code{Term3b}, \code{Term3}, \code{loglik}).
#'
#' @details
#' \itemize{
#'   \item \code{U} and \code{dd} characterize the latent correlation structure \eqn{AA'}.
#'   \item The likelihood accounts for transformations, platform-specific noise variances, and optional covariates.
#'   \item Computations are stabilized by imposing minimum thresholds on \code{b1} values.
#'   \item The function uses \code{ameans()} for fast marginalization across domains and assumes \code{Y} is ordered by stacked platforms.
#' }
#'
#' @examples
#' \dontrun{
#' # Example pseudo-code usage
#' params <- list(
#'   b0 = rnorm(5),
#'   b1 = runif(5, 0.5, 2),
#'   Beta = NULL,
#'   U = matrix(rnorm(25), nrow = 5),
#'   dd = runif(5, 0.1, 1),
#'   sigmaU2s = runif(5, 0.1, 0.5),
#'   sigma2s = runif(5, 0.1, 0.5)
#' )
#' Y <- matrix(rnorm(5 * 2 * 50), nrow = 10)
#' loglik(X = NULL, Y = Y, K = 2, params = params)
#' }
#'
#' @seealso \code{\link{ameans}}, \code{\link{tr}}, \code{\link{solve}}
#'
#' @export
loglik <- function(X, Y, K, params, b1_min=.1, verbose=FALSE){
  b0 <- as.vector(params$b0); b1 <- as.vector(params$b1)
  Beta <- params$Beta; U <- params$U; dd <- params$dd; sigmaU2s <- params$sigmaU2s; sigma2s <- params$sigma2s
  N <- ncol(Y); L <- length(dd); m <- nrow(U)
  ##  to prevent log(0) or divide by 0 error
  b1[b1 < b1_min] <- b1_min
  if (m*K != nrow(Y)) stop("The rows of input matrix Y does not equal to mxK!")
  ## actual computation starts
  Const <- det(diag(x=1+1/dd,nrow=L) -t(U)%*%(U/(1+K*sigmaU2s/sigma2s)))/prod(1+1/dd)
  ## logdet: eq:SigmaY-logdet
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
  ## H, M, U2 (\tilde{U}), and W, are defined in between eq:GAAIinv
  ## and eq:SigmaY-inv
  hh <- sigma2s/(K*sigmaU2s+sigma2s)
  D2inv <- diag(x=1/dd, nrow=L) #this is a small matrix
  ## M: eq:M-def
  M <- solve(D2inv +diag(nrow=L) -t(U)%*%(U*hh))
  U2 <- U*(sqrt(sigmaU2s)/(K*sigmaU2s+sigma2s))
  ## W <- diag(sigmaU2s/(sigma2s*(K*sigmaU2s+sigma2s))) +U2%*%(M%*%t(U2))
  ## three terms in -2*loglik, defined in eq:joint-likelihood2
  Term1 <- m*N*K*log(2*pi) +N*logdet
  Term2 <- sum(rowsums(Ycirc^2)/rep(sigma2s, K))
  Ycircbar2 <- t(U2)%*%Ycircbar
  Term3a <- K^2*sum((Ycircbar*(sigmaU2s/(sigma2s*(K*sigmaU2s+sigma2s))))*Ycircbar)
  Term3b <- K^2*tr(M%*%tcrossprod(Ycircbar2))
  Term3 <- Term3a+Term3b
  ## loglik
  loglik <- -(Term1+Term2-Term3)/2
  if (verbose) { #for debugging
    return(c(logdet=logdet, Term1=Term1, Term2=Term2, Term3a=Term3a, Term3b=Term3b, Term3=Term3, loglik=loglik))
  } else {
    return(loglik)
  }
}

#' Gradient of the Log-Likelihood with Respect to \code{b1}
#'
#' Computes the gradient of the log-likelihood function with respect to the slope parameters \code{b1},
#' assuming all other parameters are either estimated via closed-form conditional MLEs (\code{b0} and \code{Beta})
#' or moment estimators (\code{U}, \code{dd}, \code{sigma2s}, and \code{sigmaU2s}).
#' Designed for models of stacked domain-specific data with optional covariates.
#'
#' @param X Optional covariate matrix. If there are no covariates, set both \code{X} and \code{Beta} in \code{params} to \code{NULL}.
#' @param Y A numeric matrix of observed data with \code{m*K} rows and \code{N} columns, where \code{m} is the number of genes and \code{K} is the number of domains/platforms.
#' @param K Integer. Number of domains (blocks) stacked in \code{Y}.
#' @param params A list containing model parameters:
#' \itemize{
#'   \item \code{b0}: Vector of intercepts (length \code{m}).
#'   \item \code{b1}: Vector of slopes (length \code{m}).
#'   \item \code{Beta}: Coefficient matrix for covariates (optional).
#'   \item \code{U}: Eigenvector matrix for latent structure.
#'   \item \code{dd}: Corresponding eigenvalues.
#'   \item \code{sigmaU2s}: Latent platform-specific variances.
#'   \item \code{sigma2s}: Residual variances.
#' }
#' @param b1_min Numeric. Minimum allowed value for \code{b1} to avoid division by zero or \code{log(0)} errors. Defaults to \code{0.1}.
#'
#' @return A numeric vector of length equal to \code{length(b1)}, containing the gradient of the log-likelihood with respect to each slope parameter.
#'
#' @details
#' \itemize{
#'   \item Only the derivatives with respect to \code{b1} are computed.
#'   \item The structure of \code{Y} is assumed to match the stacking: \code{[Domain1; Domain2; ...; DomainK]}.
#'   \item Uses efficient array manipulations and precomputes cross-term matrices for gradient stability.
#' }
#'
#' @examples
#' \dontrun{
#' # Example pseudo-code
#' params <- list(
#'   b0 = rnorm(5),
#'   b1 = runif(5, 0.5, 2),
#'   Beta = NULL,
#'   U = matrix(rnorm(25), nrow = 5),
#'   dd = runif(5, 0.1, 1),
#'   sigmaU2s = runif(5, 0.1, 0.5),
#'   sigma2s = runif(5, 0.1, 0.5)
#' )
#' Y <- matrix(rnorm(5 * 2 * 50), nrow = 10)
#' grad <- gradLoglik(X = NULL, Y = Y, K = 2, params = params)
#' }
#'
#' @seealso \code{\link{loglik}}, \code{\link{ameans}}, \code{\link{tr}}
#'
#' @export
gradLoglik <- function(X, Y, K, params, b1_min=.1){
  b0 <- as.vector(params$b0); b1 <- as.vector(params$b1)
  Beta <- params$Beta; U <- params$U; dd <- params$dd; sigmaU2s <- params$sigmaU2s; sigma2s <- params$sigma2s
  N <- ncol(Y); L <- length(dd); m <- length(b1)/K
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
  D2inv <- diag(x=1/dd, nrow=L) #this is a small matrix
  ## M: eq:M-def
  M <- solve(D2inv +diag(nrow=L) -t(U)%*%(U*hh))
  U2 <- U*(sqrt(sigmaU2s)/(K*sigmaU2s+sigma2s))
  F <- Ycircbar*(sigmaU2s/(sigma2s*(K*sigmaU2s+sigma2s))) + (U2%*%M)%*%(t(U2)%*%Ycircbar)
  ## only produce the gradients for b1 based on eq:grad-b1. .
  Second <- Ycirc/rep(sigma2s,K) -K*(rep(1,K)%x%F)
  grad.b1 <- -N/b1 +rowsums((Y-b0)*Second)/b1^2
  return(grad.b1)
}

#' Quasi-MLE Estimation of Covariate Effects \code{Beta}
#'
#' Computes the quasi-maximum likelihood estimator (quasi-MLE) of the covariate effect matrix \code{Beta}
#' in a latent expression model with stacked domain data. If no covariates are provided (\code{X = NULL}), returns \code{NULL}.
#'
#' @param X A numeric covariate matrix. If \code{NULL}, no covariate adjustment is performed.
#' @param Y A numeric matrix of observed stacked outcomes, with \code{m*K} rows and \code{N} columns.
#' @param b1 A numeric vector of slope parameters used to scale \code{Y}. Length must be divisible by \code{K}.
#' @param K Integer. Number of domains (blocks) stacked in \code{Y}.
#' @param b1_min Numeric. Minimum allowed value for \code{b1} to avoid division by zero. Defaults to \code{0.1}.
#'
#' @return
#' If \code{X} is provided, returns a numeric matrix estimating \code{Beta} (effect of covariates on latent variables).
#' Otherwise, returns \code{NULL}.
#'
#' @details
#' \itemize{
#'   \item If \code{X} is not \code{NULL}, the function scales and centers \code{Y}, computes averaged residuals across domains,
#'         and solves a stabilized linear system to estimate \code{Beta}.
#'   \item Stability is ensured by bounding the minimum \code{b1} values and using a robust pseudo-inverse \code{rsolve2()}.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' Y <- matrix(rnorm(1000), nrow = 20)  # stacked domains
#' X <- matrix(rnorm(50 * 3), nrow = 50)  # 3 covariates
#' b1 <- runif(20, 0.5, 2)
#' Beta <- BetaEst(X = X, Y = Y, b1 = b1, K = 4)
#' }
#'
#' @seealso \code{\link{rsolve2}}, \code{\link{ameans}}
#'
#' @export
BetaEst <- function(X, Y, b1, K, b1_min=.1){
  if (is.null(X)) {
    return(NULL)
  } else {
    b1 <- as.vector(b1)
    m <- length(b1)/K; N <- ncol(Y)
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

#' Gradient Descent for Updating \code{b1} and \code{Beta}
#'
#' Performs gradient descent optimization to update the slope parameters \code{b1} and covariate coefficients \code{Beta}
#' in a latent expression model. The intercepts \code{b0} are recomputed as row means of \code{Y} at each step.
#' The eigenstructure parameters (\code{U}, \code{dd}, \code{sigma2s}, \code{sigmaU2s}) are assumed to be fixed and
#' estimated separately via \code{InitEst()}.
#'
#' @param X Optional covariate matrix. If not \code{NULL}, must be a numeric matrix with \code{ncol(X) = ncol(Y)}.
#' @param Y A numeric matrix of observed data with \code{m*K} rows and \code{N} columns, where \code{m} is the number of genes/features.
#' @param K Integer. Number of domains/platforms.
#' @param params A list containing the current parameter estimates:
#' \itemize{
#'   \item \code{b1}: Initial slope vector.
#'   \item \code{U}, \code{dd}: Eigenvectors and eigenvalues of the latent structure.
#'   \item \code{sigmaU2s}, \code{sigma2s}: Variance components.
#' }
#' @param platforms Optional character vector of platform/domain names of length \code{K}.
#' @param s0 Numeric. Initial step size for gradient descent. Defaults to \code{0.01}.
#' @param smax Numeric. Maximum allowed step size. Defaults to \code{0.5}.
#' @param sigma2.min,sigma2.max Numeric. (Currently unused) Reserved for future extensions to variance parameter updates.
#' @param b1_min,b1.max Numeric. Lower and upper bounds for \code{b1}. Defaults are \code{0.1} and \code{5}.
#' @param beta.min,beta.max Numeric. Lower and upper bounds for \code{Beta} entries. Defaults are \code{-5} and \code{5}.
#' @param max.iter Integer. Maximum number of gradient descent iterations. Defaults to \code{30}.
#' @param tol Numeric. Tolerance for convergence based on log-likelihood difference. Defaults to \code{0.1}.
#' @param verbose Logical. If \code{TRUE}, returns intermediate results including \code{Ycirc} and \code{Ycircbar}. Defaults to \code{FALSE}.
#'
#' @return
#' If \code{verbose = FALSE}, returns a list containing updated \code{params} with components:
#' \itemize{
#'   \item \code{b0}, \code{b1}, \code{Beta}, \code{U}, \code{dd}, \code{sigma2s}, \code{sigmaU2s}, \code{N}, \code{K}.
#' }
#' If \code{verbose = TRUE}, additionally returns:
#' \itemize{
#'   \item \code{Ycirc}: Adjusted residual matrix.
#'   \item \code{Ycircbar}: Averaged residual matrix across domains.
#'   \item \code{loglik}: Final log-likelihood value.
#'   \item \code{loglik.old}: Initial log-likelihood value.
#'   \item \code{iters}: Number of iterations performed.
#'   \item \code{history}: Trajectory of log-likelihood and step size per iteration.
#' }
#'
#' @details
#' The function uses a stabilized gradient descent approach with adaptive step size control.
#' If the log-likelihood decreases after a proposed update, the step size is halved until improvement is achieved or
#' the maximum number of retries is reached. All updates to \code{b1} and \code{Beta} are thresholded to specified bounds
#' for numerical stability.
#'
#' If no eigenstructure (\code{U} or \code{dd}) is available (i.e., \code{NULL}), the function skips optimization and
#' returns the input parameters.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' Y <- matrix(rnorm(1000), nrow = 20)
#' params <- list(
#'   b1 = runif(20, 0.5, 2),
#'   U = matrix(rnorm(100), nrow = 20),
#'   dd = runif(5, 0.1, 1),
#'   sigmaU2s = runif(20, 0.1, 0.5),
#'   sigma2s = runif(20, 0.1, 0.5)
#' )
#' res <- GDfun(X = NULL, Y = Y, K = 4, params = params, verbose = TRUE)
#' res$params$b1
#' }
#'
#' @seealso \code{\link{gradLoglik}}, \code{\link{loglik}}, \code{\link{BetaEst}}, \code{\link{ameans}}
#'
#' @export
GDfun <- function(X, Y, K, params, platforms=NULL,
                  s0=0.01, smax=0.5,
                  sigma2.min=0.01, sigma2.max=2,
                  b1_min=0.1, b1.max=5, beta.min = -5, beta.max = 5,
                  max.iter=30, tol=0.1, verbose=FALSE){
  # Edited by Zhining. In the case when there is only one covariate and X is not a matrix.
  if((!is.null(X)) & (!is.matrix(X))) {
    X <- matrix(X, nrow = 1)
  }

  b1 <- params$b1; U <- params$U; dd <- params$dd
  sigmaU2s <- params$sigmaU2s; sigma2s <- params$sigma2s
  gnames_params <- rownames(b1)

  # Added by Zhining 02/27/2026 ---
  idx <- which(rownames(Y) %in% gnames_params)
  Y <- Y[idx, , drop = FALSE]
  # ---
  N <- ncol(Y); m <- nrow(Y)/K;  L <- length(dd)

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
  stopifnot(all.equal(gnames_params, gnames))

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
  Norm.old <- sqrt(sum(grad.old^2)); v.old <- grad.old/Norm.old
  ll.old <- loglik(X, Y, K, params.old, b1_min=b1_min)
  ## the main loop
  j <- 1; s.j <- s0; loglik.old <- 10^5; diff.j <- 1e6
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
    Norm.j <- sqrt(sum(grad.j^2)); v.j <- grad.j/Norm.j
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
    j=j+1; ll.old <- ll.j; theta.old=theta.j; grad.old=grad.j
    Norm.old <- Norm.j; v.old <- v.j
  }
  b1 <- b1.j
  ## finishing touches
  b0 <- matrix(b0, nrow=m); b1 <- matrix(b1, nrow=m)
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
    ## 08/27/2021. I decide to output Ycirc and Ycircbar when
    ## verbose=TRUE. Ycircbar can be used to update U and dd, and is
    ## required for model integration.
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
