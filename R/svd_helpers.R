#' Fast SVD with Automatic Transposition
#'
#' A fast singular value decomposition for high-dimensional matrices.
#' Automatically transposes if rows < columns to optimize computation.
#'
#' @param X A numeric matrix.
#' @param k Number of singular vectors/values to retain. Defaults to \code{min(nrow(X), ncol(X))}.
#' @param tol Tolerance threshold for filtering small eigenvalues. Defaults to \code{1e-6}.
#'
#' @return A list with components \code{d} (singular values), \code{u}, and \code{v}.
#' @export
fsvd <- function(X, k=min(nrow(X), ncol(X)), tol=1e-6) {
  n <- nrow(X); p <- ncol(X)
  if (n<p) {
    o <- fsvd0(t(X), k=k, tol=tol)
    return(list(d=o$d, u=o$v, v=o$u))
  } else {
    return(fsvd0(X, k=k, tol=tol))
  }
}

#' Core Fast SVD Computation
#'
#' Performs fast SVD using eigendecomposition of X'X. Assumes input matrix is wide.
#'
#' @param X A numeric matrix.
#' @param k Number of singular vectors/values to retain. Defaults to \code{min(nrow(X), ncol(X))}.
#' @param tol Tolerance threshold for filtering small eigenvalues. Defaults to \code{1e-6}.
#'
#' @return A list with components \code{d} (singular values), \code{u}, and \code{v}.
#'
#' @importFrom Rfast Crossprod submatrix
#' @export
fsvd0 <- function(X, k=min(nrow(X), ncol(X)), tol=1e-6) {
  xx <- Rfast::Crossprod(X, X); n <- nrow(xx)
  a <- eigen(xx, symmetric=TRUE)
  ## remove very small (even negative) values in d due to numerical
  ## errors and/or collinearity in X
  l <- a$values[1:k]; d <- sqrt(l[l>tol]); k <- length(d)
  ## compute u and v
  v <- Rfast::submatrix(a$vectors, 1, n, 1, k)
  u <- tcrossprod(X, t(v)/d)
  return(list(d = d, u = u, v=v))
}

#' Efficient Principal Component Analysis (EPCA)
#'
#' A fast and memory-efficient PCA routine compatible with and faster than \code{hd.eigen()} in the \pkg{Rfast} package.
#' Automatically handles wide and tall matrices. Returns scaled eigenvalues, principal components, scores (optional),
#' and cumulative variance proportions. The function is especially useful for high-dimensional data where standard PCA
#' methods are too slow or memory-intensive.
#'
#' @param x A numeric matrix or data frame. Rows represent observations and columns represent variables.
#' @param center Logical. Should the variables be centered to have zero mean? Defaults to \code{TRUE}.
#' @param scale Logical. Should the variables be scaled to have unit standard deviation? Defaults to \code{FALSE}.
#' @param retx Logical. Should the principal component scores be returned? Defaults to \code{TRUE}.
#' @param k Integer. Number of principal components to compute. Defaults to \code{min(nrow(x), ncol(x))}.
#' @param tol Numeric. Eigenvalue tolerance. Components with eigenvalues below this threshold are discarded. Defaults to \code{1e-6}.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{values}}{A numeric vector of eigenvalues scaled by \(1 / (n - 1)\).}
#'   \item{\code{vectors}}{A matrix of principal component directions (eigenvectors). Columns are named \code{"PC1"}, \code{"PC2"}, etc.}
#'   \item{\code{x}}{(If \code{retx = TRUE}) The principal component scores matrix.}
#'   \item{\code{var.props}}{A numeric vector of cumulative variance proportions explained by the components.}
#' }
#'
#' @details This implementation computes the SVD or eigendecomposition of the (possibly transposed) covariance matrix,
#' depending on the dimensionality of \code{x}. For wide matrices (\code{n < m}), it avoids computing the full covariance matrix
#' to save time and memory. Eigen decomposition is preferred for numerical stability and speed.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' X <- matrix(rnorm(100 * 20), nrow = 100)
#' result <- epca(X)
#' plot(result$values, type = "b", main = "Eigenvalues")
#' }
#'
#' @seealso \code{\link[Rfast]{hd.eigen}}, \code{\link[stats]{prcomp}}, \code{\link[stats]{eigen}}
#'
#' @importFrom Rfast colVars transpose Crossprod submatrix mat.mult
#' @export
epca <- function(x, center=TRUE, scale=FALSE, retx=TRUE, k=min(nrow(x), ncol(x)), tol=1e-6) {
  n <- nrow(x); m <- ncol(x)
  var.names <- colnames(x); sample.names <- rownames(x)
  if (center) x <- eachrow(x, colmeans(x), oper="-")
  if (scale) {
    s <- Rfast::colVars(x, std = TRUE)
    x <- eachrow(x, s, oper="/")
  }
  ##
  if (n<m) { #wide matrix
    y <- Rfast::transpose(x)
    yy <- Rfast::Crossprod(y, y)
    ## As of 04/27/2024, eigen() is faster than any alternative
    ## implementations such as eigs_sym() or eigen.sym().
    a <- eigen(yy, symmetric=TRUE)
    L <- a$values[1:k]; L <- L[L>tol]; k2 <- length(L)
    if (k2<k) warning("Number of nonzero eigenvalues are less than k.")
    U <- Rfast::submatrix(a$vectors, 1, n, 1, k2)
    if (retx) X <- eachrow(U, sqrt(L))
    ## need to compute eigenvectors
    V <- tcrossprod(y, Rfast::transpose(U)*L^(-0.5))
  } else { #thin matrix, n>m
    xx <- Rfast::Crossprod(x, x)
    a <- eigen(xx, symmetric=TRUE)
    L <- a$values[1:k]; L <- L[L>tol]; k2 <- length(L)
    if (k2<k) warning("Number of nonzero eigenvalues are less than k.")
    V <- submatrix(a$vectors, 1, m, 1, k2)
    if (retx) X <- Rfast::mat.mult(x, V) #PC scores
  }
  ## final steps
  var.props <- cumsum(a$values)/sum(a$values)
  colnames(V) <- paste0("PC", 1:k2)
  rownames(V) <- var.names
  if (retx) {
    colnames(X) <- paste0("PC", 1:k2)
    rownames(X) <- sample.names
    return(list(values=L/(n-1), vectors=V, x=X, var.props=var.props))
  } else {
    return(list(values=L/(n-1), vectors=V, var.props=var.props))
  }
}


#' Generalized Probabilistic PCA (gPPCA) with Non-Scalar IID Noise
#'
#' Performs a generalized form of probabilistic PCA allowing for heterogeneous (non-scalar) noise variances across variables.
#' The method estimates latent components and variable-specific noise variances using an iterative EM-like algorithm. It is compatible
#' with \code{svd} and \code{prcomp} data formats and automatically handles centering and scaling. When the number of latent components
#' is zero, the function returns noise variances only.
#'
#' @param Y A numeric matrix of size \eqn{N \times m}, where rows are samples and columns are features.
#' @param center Logical. Should the columns of \code{Y} be centered? Defaults to \code{TRUE}.
#' @param scale Logical. Should the columns of \code{Y} be standardized? Defaults to \code{FALSE}.
#' @param retx Logical. Currently unused. Kept for compatibility. Defaults to \code{TRUE}.
#' @param L Integer. Number of latent dimensions (principal components) to extract. Defaults to \code{min(nrow(Y), ncol(Y))}.
#' @param min.sigma2 Minimum allowed residual variance per variable. Used for numerical stability. Defaults to \code{0.01}.
#' @param pca.tol Tolerance for PCA eigenvalue thresholding. Used during initialization. Defaults to \code{1e-6}.
#' @param iter.tol Convergence tolerance for the iterative updates. Defaults to \code{1e-6}.
#' @param max.iter Maximum number of iterations allowed for convergence. Defaults to \code{20}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{U}}{Matrix of estimated loading directions (eigenvectors of the generalized correlation matrix).}
#'   \item{\code{dd}}{Estimated eigenvalues of the latent space (diagonal elements of the covariance matrix).}
#'   \item{\code{sigma2s}}{Estimated residual variances for each feature (non-scalar noise).}
#'   \item{\code{history}}{A numeric vector recording the convergence error at each iteration.}
#' }
#'
#' @details The algorithm starts with an initial PCA estimate, infers variable-specific noise variances,
#' and iteratively refines both the loading matrix and noise model using generalized singular value decomposition.
#' If \code{L = 0}, the function estimates only per-variable noise variances. For standardized input,
#' the model ensures estimated components correspond to a correlation structure.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' Y <- matrix(rnorm(100 * 20), nrow = 100)
#' res <- gppca(Y, L = 5)
#' plot(res$dd, type = "b", main = "Generalized Eigenvalues")
#' }
#'
#' @seealso \code{\link{epca}}, \code{\link[stats]{prcomp}}, \code{\link[stats]{svd}}
#'
#' @importFrom Rfast colVars transpose Crossprod submatrix mat.mult
#' @export
gppca <- function(Y, center=TRUE, scale=FALSE, retx=TRUE, L=min(nrow(Y), ncol(Y)), min.sigma2=0.01, pca.tol=1e-6, iter.tol=1e-6, max.iter=20) {
  N <- nrow(Y); m <- ncol(Y); m2 <- min(N,m)
  var.names <- colnames(Y); sample.names <- rownames(Y)
  ## 1. Data pre-processing
  if (center) {
    Y <- eachrow(Y, colmeans(Y), oper="-")
    N <- N-1 #to reduce bias
  }
  if (scale) {
    s <- Rfast::colVars(Y, std = TRUE)
    Y <- eachrow(Y, s, oper="/")
  }
  dimnames(Y) <- list(sample.names, var.names)
  ## 1b. A very special case, L=0
  if (L==0) {
    if (scale) {
      sigma2s <- rep(1,m)
    } else {
      sigma2s <- colsums(Y^2)/N
    }
    U <- dd <- history <- NULL
  } else { #L>0
    ## 2. Initial PPCA of Y, starting from the standard PCA
    ee <- epca(Y, center=FALSE, scale=FALSE, retx=FALSE, k=L, tol=pca.tol)
    lambdasY <- ee$values; U <- ee$vectors[, 1:L, drop=FALSE]
    ## it is possible that length(ee$values) is smaller than the
    ## specified k, so we need to update it just in case.
    if (L != length(lambdasY)) {
      L <- length(lambdasY)
      warning(paste0("The maximum effective number of PCs is less than the specified value. This parameter is set to be L=", L, "."))
    }
    ## 2b. Using the remaining variance to compute sigma2
    vp <- ee$var.props[L]
    if (L==m2) { #no remaining variance and m2-k=0
      sigma2 <- min.sigma2
    } else {
      v.remain <- sum(lambdasY)*(1-vp)/vp
      sigma2 <- max(v.remain/(m-L), min.sigma2)
    }
    dd <- lambdasY-sigma2
    ## 3. Optional: standardization
    if (scale) {
      s0 <- drop((U*U)%*%dd)+sigma2
      U <- U/sqrt(s0); sigma2s <- sigma2/s0
      ## sample variance of Y
      sigmaY2 <- rep(1, m)
    } else {
      sigma2s <- rep(sigma2, m)
      sigmaY2 <- colsums(Y^2)/N
    }
    ## Start of the iteration.
    i <- 0; err <- Inf; history <- c()
    while (i <= max.iter & err > iter.tol) {
      ## 4. Update sigma2s
      diagAA <- drop((U*U)%*%dd)
      sigma2s.new <- pmax(sigmaY2-diagAA, min.sigma2)
      ## 5. Update AA' (U and dd)
      Ytilde <- sweep(Y, 2, sqrt(sigma2s.new), "/")
      ss <- fsvd(Ytilde, k=L, tol=pca.tol)
      ## ss$v is U_{\tilde{Y},L}^{(k)}
      U.new <- ss$v*sqrt(sigma2s.new)
      ## ss$d^2/N is \lambda_{AA',l}^{(k)}
      dd.new <- pmax(ss$d^2/N, 0)
      ## 6. Stopping rule
      Bkk <- tcrossprod(sqrt(dd.new))*crossprod(U.new)
      Bkminus1 <- tcrossprod(sqrt(dd))*crossprod(U)
      Bkkminus1 <- (sqrt(dd.new)%*%t(sqrt(dd)))*(t(U.new)%*%U)
      Term1 <- sum(Bkk^2)+sum(Bkminus1^2)-2*sum(Bkkminus1^2)
      err <- (Term1 +sum((sigma2s.new-sigma2s)^2))/m^2
      U <- U.new; dd <- dd.new; sigma2s <- sigma2s.new
      i <- i+1
      history <- c(history, err)
    }
    ## 7. Final steps
    Astar <- sweep(U, 2, sqrt(dd), "*")
    ## 7b. Standardize Astar so that the outcome is correlation matrix
    if (scale) {
      s0 <- rowsums(Astar^2)+sigma2s
      Astar <- Astar/sqrt(s0); sigma2s <- sigma2s/s0
    }
    ss <- fsvd(Astar)
    U <- ss$u; dd <- ss$d^2
    rownames(U) <- names(sigma2s) <- var.names
    colnames(U) <- names(dd) <- paste0("PC", 1:ncol(U))
    if (sum(dd<=0)>0) warning("Some estimated eigenvalues are zero or negative. Please consider using a smaller L (number of PCs)!")
  }
  return(list(U=U, dd=dd, sigma2s=sigma2s, history=history))
}

