#' Generate Cross-Validation Splits with Approximately Equal Fold Sizes
#'
#' Creates a list of indices that partition \code{1:N} into approximately equal-sized subsets for use in K-fold cross-validation.
#' The function ensures that the folds are randomly ordered and differ in size by at most one.
#'
#' @param N Integer. Total number of samples to be split.
#' @param n.folds Integer. Number of cross-validation folds. Defaults to \code{5}.
#'
#' @return A list of length \code{n.folds}, where each element is an integer vector of indices for one fold.
#'
#' @details
#' The function ensures:
#' \itemize{
#'   \item Fold sizes differ by at most one.
#'   \item Samples are randomly shuffled before assignment to folds.
#'   \item Useful in K-fold cross-validation for supervised learning tasks.
#' }
#'
#' @examples
#' set.seed(123)
#' folds <- cv.split(23, n.folds = 5)
#' sapply(folds, length)  # check size distribution
#' unlist(folds)          # ensure full coverage of 1:23
#'
#' @seealso \code{\link[stats]{sample}}, \code{\link[caret]{createFolds}} in other ML packages
#'
#' @export
cv.split <- function(N, n.folds=5) {
  n0 <- N %/% n.folds #n0 is (approx) the sample size for each subset
  r <- N %% n.folds   #remainder of the division
  ns <- c(rep(n0+1, r), rep(n0, n.folds-r)) #samples in each subset
  i.start <- c(1, cumsum(ns[-n.folds])+1)
  i.end <- cumsum(ns)
  ## now generate the shuffled data
  shuffle <- sample(1:N)
  return(lapply(1:n.folds, function(k) shuffle[i.start[k]:i.end[k]]))
}

#' Cross-Validation for Selecting Optimal Latent Dimension \code{L}
#'
#' Performs K-fold cross-validation across a grid of candidate latent dimensions (\code{L})
#' to select the best number of latent factors for initializing the PPCA-Xnorm model using \code{InitEst()}.
#'
#' @param X Optional covariate matrix. If no covariates are used, set \code{X = NULL}.
#' @param Y A numeric matrix of observed stacked outcomes, with \code{m*K} rows and \code{n} columns.
#' @param K Integer. Number of domains/platforms.
#' @param folds A list of test-set column indices for each CV fold (e.g., output of \code{cv.split()}).
#' @param k.source Integer. Index (1-based) of the source platform within the stacked rows of \code{Y}.
#' @param k.target Integer. Index (1-based) of the target platform within the stacked rows of \code{Y}.
#' @param Ls Integer vector. Candidate grid of latent dimensions \code{L} to evaluate. Defaults to \code{1:10}.
#' @param ... Additional arguments passed to \code{InitEst()} and \code{Prediction()}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{MSE}: A matrix of mean squared errors, with rows corresponding to each candidate \code{L} and columns to CV folds.
#'   \item \code{Lstar}: The value of \code{L} that minimizes the average MSE across folds.
#' }
#'
#' @details
#' \itemize{
#'   \item For each \code{L} and fold, \code{InitEst()} is fitted on the training set, predictions are made on the test set using \code{Prediction()},
#'         and mean squared error (MSE) is evaluated.
#'   \item The optimal \code{Lstar} is chosen to minimize the mean MSE across folds.
#'   \item If covariates are present, they are split accordingly between training and testing.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' Y <- matrix(rnorm(20 * 50), nrow = 20)
#' folds <- cv.split(ncol(Y), n.folds = 5)
#' res <- CV.InitEst(X = NULL, Y = Y, K = 4, folds = folds, k.source = 1, k.target = 2)
#' res$Lstar
#' }
#'
#' @seealso \code{\link{InitEst}}, \code{\link{Prediction}}, \code{\link{cv.split}}, \code{\link{MSE1}}
#'
#' @export
CV.InitEst <- function(X, Y, K, folds, k.source, k.target, Ls=1:10, ...) {
  n <- ncol(Y); m <- nrow(Y)/K
  nFolds <- length(folds); nL <- length(Ls)
  MSEs <- matrix(-1, nL, nFolds)
  rownames(MSEs) <- paste0("L=", Ls)
  colnames(MSEs) <- paste0("fold", 1:nFolds)
  ##
  fun <- function(k) {
    test.idx <- folds[[k]]; train.idx <- setdiff(1:n, test.idx)
    Y.test <- Y[, test.idx, drop=FALSE]; Y.train <- Y[, train.idx, drop=FALSE]
    if (!is.null(X)) {
      X.test <- X[, test.idx, drop=FALSE]; X.train <- X[, train.idx, drop=FALSE]
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
  ##
  for (l in 1:nL) {
    L <- Ls[l]
    MSEs[l,] <- sapply(1:nFolds, fun)
  }
  Lstar <- Ls[which.min(rowMeans(MSEs))]
  return(list(MSE=MSEs, Lstar=Lstar))
}
