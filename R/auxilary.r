## useful auxiliary functions for the PPCA-XPN project.

#' Matrix Visualization Plot
#'
#' Visualizes a numeric matrix using a grayscale heatmap.
#'
#' @param mat Numeric matrix to plot.
#' @param ... Additional parameters passed to \code{image}.
#' @return None. Produces a heatmap plot.
#' @export
#' @examples
#' mat <- matrix(runif(100), 10, 10)
#' mplot(mat)
mplot <- function(mat, ...) {
  # Plot matrix with reversed y-axis so row 1 is at the top.
  image(t(mat[nrow(mat):1,]),
        col=grey(seq(1, 0, length.out=101)),
        xaxt="n",
        yaxt="n", ...)
}

#' Matrix Trace
#'
#' Calculates the sum of the diagonal elements of a square matrix.
#'
#' @param M A numeric matrix.
#' @return Numeric scalar, the matrix trace.
#' @export
#' @examples
#' M <- matrix(1:9, 3, 3)
#' tr(M)
tr <- function(M){
  sum(diag(as.matrix(M)))
}

#' Flatten Matrix to Vector
#'
#' Converts a matrix into a single vector in column-major order.
#'
#' @param x Matrix or array.
#' @return Numeric vector.
#' @export
#' @examples
#' mat <- matrix(1:6, 2, 3)
#' vec(mat)
vec <- as.vector

#' Fast Singular Value Decomposition (SVD)
#'
#' Efficient computation of thin SVD, optimized for tall and wide matrices.
#' Modeled after hd.eigen() in Rfast() but even faster and more convenient.
#'
#' Automatically handles cases where rows > columns or vice versa.
#'
#' @param X Numeric matrix.
#' @param k Integer number of components to retain.
#' @param tol Tolerance threshold for small eigenvalues.
#' @return List with components \code{d} (singular values), \code{u} (left vectors), \code{v} (right vectors).
#' @export
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' fsvd(X, k=5)
fsvd <- function(X, k=min(nrow(X), ncol(X)), tol=1e-6) {
  n <- nrow(X)
  p <- ncol(X)
  if (n<p) {
    o <- fsvd0(t(X), k=k, tol=tol)
    return(list(d=o$d, u=o$v, v=o$u))
  } else {
    return(fsvd0(X, k=k, tol=tol))
  }
}

#' Core Fast SVD for Thin Matrix
#'
#' Performs SVD on thin matrices using cross-product and eigen decomposition.
#'
#' @param X Numeric matrix.
#' @param k Number of singular values to retain.
#' @param tol Tolerance for small values.
#' @return List with components \code{d}, \code{u}, \code{v}.
#' @export
fsvd0 <- function(X, k=min(nrow(X), ncol(X)), tol=1e-6) {
  xx <- Rfast::Crossprod(X, X)
  n <- nrow(xx)
  a <- eigen(xx, symmetric=TRUE)
  # Remove very small (even negative) values in d due to numerical errors and/or collinearity in X
  l <- a$values[1:k]
  d <- sqrt(l[l>tol])
  k <- length(d)
  # Compute u and v
  v <- Rfast::submatrix(a$vectors, 1, n, 1, k)
  u <- tcrossprod(X, t(v)/d)
  return(list(d = d, u = u, v=v))
}

#' Efficient Principal Component Analysis (PCA)
#'
#' Fast PCA for wide or tall matrices using eigen decomposition.
#'
#' Automatically handles thin and wide matrices and optionally returns scores.
#'
#' Compatible and faster than hd.eigen() in Rfast. Note that var.props provides
#' the full list of variance proportions, typically much longer than k.
#'
#' @param x Data matrix.
#' @param center Logical, whether to center columns.
#' @param scale Logical, whether to scale columns.
#' @param retx Logical, return PC scores.
#' @param k Number of principal components to retain.
#' @param tol Eigenvalue tolerance threshold.
#' @return List containing eigenvalues, eigenvectors, variance proportions, and optionally scores.
#' @export
#' @examples
#' x <- matrix(rnorm(100), 10, 10)
#' epca(x, k=3)
epca <- function(x, center=TRUE, scale=FALSE, retx=TRUE, k=min(nrow(x), ncol(x)), tol=1e-6) {
  n <- nrow(x)
  m <- ncol(x)
  var.names <- colnames(x)
  sample.names <- rownames(x)

  if (center) x <- eachrow(x, colmeans(x), oper="-")
  if (scale) {
    s <- Rfast::colVars(x, std = TRUE)
    x <- eachrow(x, s, oper="/")
  }

  if (n<m) { # Wide matrix
    y <- Rfast::transpose(x)
    yy <- Rfast::Crossprod(y, y)
    a <- eigen(yy, symmetric=TRUE)
    L <- a$values[1:k]
    L <- L[L>tol]
    k2 <- length(L)
    if (k2<k) {
      warning("Number of nonzero eigenvalues are less than k.")
    }
    U <- Rfast::submatrix(a$vectors, 1, n, 1, k2)
    if (retx) {
      X <- eachrow(U, sqrt(L))
    }
    ## need to compute eigenvectors
    V <- tcrossprod(y, Rfast::transpose(U)*L^(-0.5))
  } else { # Thin matrix, n > m
    xx <- Rfast::Crossprod(x, x)
    a <- eigen(xx, symmetric=TRUE)
    L <- a$values[1:k]
    L <- L[L>tol]
    k2 <- length(L)
    if (k2<k) {
      warning("Number of nonzero eigenvalues are less than k.")
    }
    V <- submatrix(a$vectors, 1, m, 1, k2)
    if (retx) {
      X <- Rfast::mat.mult(x, V) # PC scores
    }
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

#' Generalized Probabilistic PCA (GPPCA)
#'
#' Extended probabilistic PCA for data with non-scalar iid noise.
#'
#' U == U_{AA'}, dd == d_{l}^{2}. Y is an (N x m)-dim matrix (compatible with svd/prcomp)
#'
#' @param Y Data matrix.
#' @param center Center columns.
#' @param scale Scale columns.
#' @param retx Return scores.
#' @param L Number of latent factors.
#' @param min.sigma2 Minimum noise variance.
#' @param pca.tol Tolerance for PCA convergence.
#' @param iter.tol Iteration tolerance for updates.
#' @param max.iter Maximum iterations.
#' @return List with U, dd, sigma2s, history (tracking error convergence).
#' @export
#' @examples
#' Y <- matrix(rnorm(100), 10, 10)
#' gppca(Y, L=3)
gppca <- function(Y, center=TRUE, scale=FALSE, retx=TRUE, L=min(nrow(Y), ncol(Y)), min.sigma2=0.01, pca.tol=1e-6, iter.tol=1e-6, max.iter=20) {
  N <- nrow(Y)
  m <- ncol(Y)
  m2 <- min(N,m)
  var.names <- colnames(Y)
  sample.names <- rownames(Y)

  ## 1. Data pre-processing
  if (center) {
    Y <- eachrow(Y, colmeans(Y), oper="-")
    N <- N-1 # to reduce bias
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
  } else { # L>0
    ## 2. Initial PPCA of Y, starting from the standard PCA
    ee <- epca(Y, center=FALSE, scale=FALSE, retx=FALSE, k=L, tol=pca.tol)
    lambdasY <- ee$values
    U <- ee$vectors[, 1:L, drop=FALSE]

    ## It is possible that length(ee$values) is smaller than the
    ## specified k, so we need to update it just in case.
    if (L != length(lambdasY)) {
      L <- length(lambdasY)
      warning(paste0("The maximum effective number of PCs is less than the specified value. This parameter is set to be L=", L, "."))
    }

    ## 2b. Using the remaining variance to compute sigma2
    vp <- ee$var.props[L]
    if (L==m2) { # no remaining variance and m2-k=0
      sigma2 <- min.sigma2
    } else {
      v.remain <- sum(lambdasY)*(1-vp)/vp
      sigma2 <- max(v.remain/(m-L), min.sigma2)
    }
    dd <- lambdasY-sigma2

    ## 3. Optional: standardization
    if (scale) {
      s0 <- drop((U*U)%*%dd)+sigma2
      U <- U/sqrt(s0)
      sigma2s <- sigma2/s0

      ## sample variance of Y
      sigmaY2 <- rep(1, m)
    } else {
      sigma2s <- rep(sigma2, m)
      sigmaY2 <- colsums(Y^2)/N
    }

    ## Start of the iteration.
    i <- 0
    err <- Inf
    history <- c()
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
      U <- U.new
      dd <- dd.new
      sigma2s <- sigma2s.new
      i <- i+1
      history <- c(history, err)
    }

    ## 7. Final steps
    Astar <- sweep(U, 2, sqrt(dd), "*")

    ## 7b. Standardize Astar so that the outcome is correlation matrix
    if (scale) {
      s0 <- rowsums(Astar^2)+sigma2s
      Astar <- Astar/sqrt(s0)
      sigma2s <- sigma2s/s0
    }
    ss <- fsvd(Astar)
    U <- ss$u
    dd <- ss$d^2
    rownames(U) <- names(sigma2s) <- var.names
    colnames(U) <- names(dd) <- paste0("PC", 1:ncol(U))
    if (sum(dd<=0)>0) warning("Some estimated eigenvalues are zero or negative. Please consider using a smaller L (number of PCs)!")
  }
  return(list(U=U, dd=dd, sigma2s=sigma2s, history=history))
}

#' Robust Hat Matrix Computation
#'
#' Computes a robust projection matrix (hat matrix) with SVD truncation.
#'
#' A singular value threshold is defined. Only those eigenvalues that passed this
#' threshold are used in the hat matrix.
#'
#' @param X Design matrix.
#' @param d.prop Proportion cutoff for small singular values.
#' @param dmin Minimum allowable singular value.
#' @return Hat matrix.
#' @export
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' rhat(X)
rhat <- function(X, d.prop=1e-6, dmin=1e-9) {
  N <- nrow(X)
  o <- svd(X)
  thresh <- max(sum(o$d) * d.prop, dmin)
  idx <- which(o$d>thresh)
  if (length(idx)==0) { # nothing is left
    return(matrix(0, N, N))
  } else {
    return(tcrossprod(o$u[, idx]))
  }
}

#' Robust Inverse of X'X for Thin Matrices
#'
#' Computes a robust regularized inverse of \code{X'X}, handling near-singular cases.
#'
#' X is a thin matrix with dimension n x p (n >> p).
#'
#' @param X Numeric matrix (tall and thin ideally).
#' @param d.prop Proportion cutoff for small singular values.
#' @param dmin Minimum singular value.
#' @param dmax Maximum allowed singular value (cap for stability).
#' @return Inverted matrix (p x p).
#' @export
#' @examples
#' X <- matrix(rnorm(100), 10, 10)
#' rsolve2(X)
rsolve2 <- function(X, d.prop=1e-6, dmin=1e-9, dmax=1e9){
  if (!is.matrix(X)) stop("X must be a matrix, not a dataframe or vector.")

  N <- nrow(X)
  p <- ncol(X)

  if (N <= p) warning("N is less or equal to p in rsolve2(); the results may be numericall unstable.")

  if (N==1){ # Special case for 1 row matrix (single sample)
    d <- sqrt(sum(X^2))
    v <- as.vector(X)/d
    thresh <- min(max(d * d.prop, dmin), dmax)
    inv.mat <- (1/d^2 -1/thresh^2)*tcrossprod(v) +(1/thresh^2)*diag(nrow=p)
  } else if (p==1) { # Special case for single-column X (column vector)
    d <- sqrt(sum(X^2))
    thresh <- min(max(d * d.prop, dmin), dmax)
    inv.mat <- as.matrix(1/max(d,thresh)^2)
  } else { ## general case, N,p>1
    o <- svd(X)
    # Singular value regularization
    thresh <- min(max(sum(o$d) * d.prop, dmin), dmax)
    inv.mat <- o$v %*% diag(1/(pmax(o$d, thresh))^2) %*% t(o$v)
  }
  return(inv.mat)
}

#' Fast Apply Mean Along Array Dimension
#'
#' Efficiently calculates the mean along a specified dimension of an array.
#'
#' This is a faster alternative to `apply(an_array, c(1,3), mean)`.
#'
#' @param a Numeric array.
#' @param i Integer, dimension index along which to compute the mean.
#' @return Array with means along dimension \code{i}.
#' @export
#' @examples
#' a <- array(runif(24), dim=c(3,4,2))
#' ameans(a, 3)
ameans <- function(a, i) {
  n <- length(dim(a))
  # Permute array dimensions so 'i' is last
  b <- aperm(a, c(seq_len(n)[-i], i))
  # Take rowMeans over permuted array (last dimension gets reduced)
  rowMeans(b, dims = n - 1)
}

#' Calculate Mean Squared Error (MSE)
#'
#' Computes MSE between two vectors or matrices.
#'
#' @param est Estimated values.
#' @param truth True values.
#' @param relative Logical, if TRUE compute relative MSE.
#' @return Numeric MSE.
#' @export
#' @examples
#' MSE1(c(1,2,3), c(1,2,2.5))
MSE1  <-  function(est, truth, relative=FALSE){
  total.mse <- sum((est - truth)^2, na.rm=TRUE)
  if (relative) {
    return(total.mse/sum(truth^2, an.rm=TRUE))
  } else {
    return(total.mse/length(as.vector(truth)))
  }
}

#' Winsorize Matrix with Median Imputation for NAs
#'
#' Imputes missing values and applies winsorization to control outliers.
#'
#' Median imputation for NAs; then winsorize outliers by the 5STD rule.
#' Note that this function only works for numeric, matrix-like data
#'
#' @param mydata Numeric matrix.
#' @param n.std Number of standard deviations for winsorization cutoff.
#' @param trim Proportion of values trimmed when computing mean/std.
#' @param verbose Logical, if TRUE prints summary.
#' @return Winsorized matrix.
#' @export
#' @examples
#' mat <- matrix(rnorm(100), 10, 10)
#' mat[1, 1] <- NA
#' winsor(mat)
winsor <- function(mydata, n.std=5, trim=0.1, verbose=FALSE) {
  mydata <- as.matrix(mydata)

  # Median imputation for missing values
  n.na <- sum(is.na(mydata))
  if (n.na>0) {
    na.rows <- which(rowSums(is.na(mydata))>0)
    for (i in na.rows) mydata[i,is.na(mydata[i,])] <- median(mydata[i,], na.rm=TRUE)
  }

  # Compute trimmed mean and standard deviation for each row and then lower/upper bounds
  ysort <- Rfast::rowSort(mydata)
  n <- ncol(mydata)
  nL <- ceiling(n*trim)
  nU <- floor(n*(1-trim))
  ytrimmed <- ysort[, seq(nL, nU)]
  mus <- rowMeans(ytrimmed)
  sds <- rowVars(ytrimmed, std=TRUE)

  # Winsorization bounds
  Ls <- mus-n.std*sds
  Us <- mus+n.std*sds

  # Apply winsorization row-wise
  y <- t(sapply(1:nrow(mydata), function(i) {
    x <- mydata[i,]
    L <- Ls[i]
    U <- Us[i]
    normalRange <- range(x[x>=L & x<=U])

    # Clamp values outside range
    x[x<L] <- normalRange[1]
    x[x>U] <- normalRange[2]
    return(x)
  }))
  rownames(y) <- rownames(mydata)

  # Verbose summary if required
  if (verbose) {
    n.out <- sum(y-mydata != 0)
    N <- nrow(mydata)*ncol(mydata)
    print(paste0(n.na, " (", round(100*n.na/N,2), "%) NAs were replaced by gene-specific trimmed means. ",
                 n.out, " (", round(100*n.out/N, 2), "%) outliers were replaced by trimmed mean +/-", n.std, "*STDs." ))
  }
  return(y)
}

#' Predict Target Expression From MatchMixeR/OLS Model
#'
#' Generates predicted expression values using estimated slope and intercept.
#'
#' @param Ysource Matrix of source expression data.
#' @param trained.model A trained model object containing \code{betamat}.
#' @return Matrix of predicted target expression values.
#' @export
#' @examples
#' trained.model <- list(betamat = matrix(c(0.5, 1.2), ncol=2, dimnames=list(NULL, c("Intercept", "Slope"))))
#' Ysource <- matrix(runif(20), 5, 4)
#' predict.mm(Ysource, trained.model)
predict.mm <- function(Ysource, trained.model) {
  bb <- trained.model$betamat
  b0 <- bb[, "Intercept"]
  b1 <- bb[, "Slope"]
  Ytarget <- b0 +sweep(Ysource, 1, b1, "*")
  return(Ytarget)
}

#' Create Cross-Validation Splits
#'
#' Divides data indices into approximately equal-sized subsets for k-fold cross-validation.
#'
#' @param N Total number of samples.
#' @param n.folds Number of cross-validation folds.
#' @return List of length \code{n.folds}, each containing sample indices.
#' @export
#' @examples
#' cv.split(100, 5)
cv.split <- function(N, n.folds=5) {
  n0 <- N %/% n.folds # n0 is (approx) the sample size for each subset
  r <- N %% n.folds   # remainder of the division
  ns <- c(rep(n0+1, r), rep(n0, n.folds-r)) # samples in each subset

  # Start and end indices for each fold
  i.start <- c(1, cumsum(ns[-n.folds])+1)
  i.end <- cumsum(ns)

  # Shuffle and split indices
  shuffle <- sample(1:N)
  return(lapply(1:n.folds, function(k) shuffle[i.start[k]:i.end[k]]))
}

#' Cross-Validation Based Initial Estimation
#'
#' Performs k-fold cross-validation to select optimal L for latent factor models.
#'
#' @param X Optional predictor matrix (NULL for no covariate case).
#' @param Y Response matrix.
#' @param K Number of source-target pairs.
#' @param folds List of folds created by \code{cv.split}.
#' @param k.source Index for source samples.
#' @param k.target Index for target samples.
#' @param Ls Vector of candidate latent factor counts.
#' @param ... Additional arguments passed to \code{InitEst}.
#' @return List containing cross-validation MSE matrix and optimal L.
#' @export
CV.InitEst <- function(X = NULL, Y, K, folds, k.source, k.target, Ls=1:10, ...) {
  n <- ncol(Y)
  m <- nrow(Y)/K
  nFolds <- length(folds)
  nL <- length(Ls)

  MSEs <- matrix(-1, nL, nFolds)
  rownames(MSEs) <- paste0("L=", Ls)
  colnames(MSEs) <- paste0("fold", 1:nFolds)

  # Fold evaluation function
  fun <- function(k) {
    test.idx <- folds[[k]]
    train.idx <- setdiff(1:n, test.idx)
    Y.test <- Y[, test.idx, drop=FALSE]
    Y.train <- Y[, train.idx, drop=FALSE]
    if (!is.null(X)) {
      X.test <- X[, test.idx, drop=FALSE]
      X.train <- X[, train.idx, drop=FALSE]
    } else {
      X.test <- X.train <- NULL
    }
    rr <- InitEst(X=X.train, Y=Y.train, K=K, L=L, ...)

    ## extract Ysource and Ytarget from Y
    Ysource.test <- Y.test[((k.source-1)*m+1):(k.source*m),]
    Ytarget.test <- Y.test[((k.target-1)*m+1):(k.target*m),]
    Ypred <- Prediction(Ysource=Ysource.test, X=X.test, trained.model=rr,
                        k.source=k.source, k.target=k.target)
    return(MSE1(Ypred, Ytarget.test))
  }

  for (l in 1:nL) {
    L <- Ls[l]
    MSEs[l,] <- sapply(1:nFolds, fun)
  }

  Lstar <- Ls[which.min(rowMeans(MSEs))]
  return(list(MSE=MSEs, Lstar=Lstar))
}

#' Extract Estimated Covariance Matrix Sigma_Y
#'
#' @param est Estimated model object from gppca or similar.
#' @param var.only If TRUE, return only diagonal elements.
#' @return Estimated covariance matrix.
#' @export
getSigmaY <- function(est, var.only=FALSE) {
  U <- est$U
  dd <- est$dd
  sigmaU2s <- est$sigmaU2s
  sigma2s <- est$sigma2s

  # Make it work with the results produced by gppca()
  if (is.null(sigmaU2s)) sigmaU2s <- 0

  ## var.only ==> only need to produce the diagonal elements
  if (var.only) {
    if (is.null(dd)) {
      diagAA <- 0
    } else {
      diagAA <- drop((U*U)%*%dd)
    }
    return(sigmaU2s+sigma2s+diagAA)
  } else {
    if (is.null(dd)) { # L=0
      AA <- 0
    } else { # L>0
      AA <- U %*% diag(x=dd,nrow=length(dd)) %*% t(U)
    }
    return(diag(x=sigmaU2s+sigma2s)+AA)
  }
}

#' Mean Squared Error for Sigma_Y
#'
#' Calculates MSE between estimated and oracle covariance matrices.
#'
#' This function is efficient for very large A.
#'
#' @param est Estimated model.
#' @param oracle Oracle model.
#' @return Numeric MSE value.
#' @export
MSE.SigmaY <- function(est, oracle) {
  A <- sweep(est$U, 2, sqrt(est$dd), "*")
  m <- nrow(A)
  B <- sweep(oracle$U, 2, sqrt(oracle$dd), "*")
  ss.A <- est$sigmaU2s+est$sigma2s
  ss.B <- oracle$sigmaU2s+oracle$sigma2s
  ##
  Term1 <- sum(crossprod(A)^2)+sum(crossprod(B)^2)-2*sum((crossprod(A, B))^2)
  Term2 <- sum((ss.A-ss.B)^2)
  return((Term1+Term2)/m^2)
}

#' Differential Expression Analysis via Limma
#'
#' Fits linear models and performs empirical Bayes moderation using limma.
#'
#' @param gdata Gene expression matrix.
#' @param v Design matrix or data frame.
#' @param padj.method Multiple testing correction method.
#' @return Data frame with coefficients, t-stats, raw and adjusted p-values.
#' @export
limma <- function(gdata, v, padj.method="BH"){
  v <- as.matrix(v)
  vn <- colnames(v)
  if (is.null(vn)) {
    vn <- paste0("X", 1:ncol(v))
    colnames(v) <- vn
  }
  ## remove missing samples
  na.id <- apply(v, 1, function(x) any(is.na(x)))
  gdata <- gdata[, !na.id]
  v <- v[!na.id,,drop=FALSE]
  design.mat <- cbind(Intercept=1, v)
  efit <- eBayes(lmFit(gdata, design.mat))
  betahat <- efit$coefficients
  colnames(betahat) <- paste0("betahat.", colnames(betahat))
  tstats=efit$t[,-1, drop=FALSE]
  colnames(tstats) <- paste0("tstat.", colnames(tstats))
  pvals <- efit$p.value[,-1, drop=FALSE]
  colnames(pvals) <- paste0("pvals.", vn)
  adjP <- matrix(p.adjust(pvals, padj.method), nrow=nrow(pvals))
  colnames(adjP) <- paste0("adjP.", vn)
  return(data.frame(betahat, tstats, pvals, adjP))
}

#' Fisher's Combined P-value Test
#'
#' Combines independent p-values using Fisher's method.
#'
#' @param ps Vector of p-values.
#' @param pmin Minimum allowable p-value.
#' @return Combined p-value.
#' @export
#' @examples
#' fisher(c(0.01, 0.02, 0.05))
fisher <- function(ps, pmin=1e-6) {
  S <- -2*sum(log(pmax(ps,pmin)))
  return(1-pchisq(S, 2*length(ps)))
}


