#' Cross-Platform Prediction of Gene Expression
#'
#' Predicts gene expression levels on a target platform based on source platform data using a trained MatchMixeR model.
#' Accounts for individualized scaling, covariate adjustment, shared latent structure, and platform-specific noise.
#'
#' @param Ysource A numeric matrix of source platform expression data. Rows are genes/features, and columns are samples.
#' @param X Optional covariate matrix for the source data. If no covariates are available, set \code{X = NULL}.
#' @param trained.model A list containing trained model parameters, typically produced by MatchMixeR gradient descent or estimation functions:
#' \itemize{
#'   \item \code{b0}: Matrix of intercepts (genes × platforms).
#'   \item \code{b1}: Matrix of slopes (genes × platforms).
#'   \item \code{Beta}: Covariate coefficients (optional).
#'   \item \code{U}: Eigenvectors of shared latent structure.
#'   \item \code{dd}: Corresponding eigenvalues (singular values squared).
#'   \item \code{sigmaU2s}: Gene-specific latent variances.
#'   \item \code{sigma2s}: Gene-specific residual variances.
#' }
#' @param k.source Integer. Index (1-based) of the source platform in \code{b0} and \code{b1}. Defaults to \code{1}.
#' @param k.target Integer. Index (1-based) of the target platform. Defaults to \code{2}.
#' @param b1_min Numeric. Minimum threshold for slope parameters \code{b1}. Genes with slopes smaller than this value are treated with a simplified prediction. Defaults to \code{0.02}.
#' @param min.sigmaU2 Numeric. Minimum threshold for latent variance \code{sigmaU2s}. Defaults to \code{0.05}.
#' @param verbose Logical. If \code{TRUE}, returns detailed prediction components. Defaults to \code{FALSE}.
#'
#' @return
#' If \code{verbose = FALSE}, returns a numeric matrix \code{Yhat} of predicted expressions on the target platform.
#' If \code{verbose = TRUE}, returns a list containing:
#' \itemize{
#'   \item \code{Yhat}: Predicted target expression matrix.
#'   \item \code{muk2}: Adjusted covariate contribution on the target platform.
#'   \item \code{GSI}: Gene-specific information contribution.
#'   \item \code{SI}: Shared latent information contribution.
#'   \item \code{ids0}: Indices of genes where simplified prediction was applied.
#' }
#'
#' @details
#' \itemize{
#'   \item Genes with very small source slopes (\code{b1k}) are predicted purely based on covariate effects and intercepts.
#'   \item Otherwise, the prediction combines gene-specific information (GSI) and shared information (SI) derived from the latent structure (\code{U} and \code{dd}).
#'   \item If missing values (\code{NA}) are present in \code{Ysource}, the function will stop. Preprocessing with \code{winsor()} is recommended.
#' }
#'
#' @examples
#' \dontrun{
#' # Example pseudo-usage
#' set.seed(123)
#' Ysource <- matrix(rnorm(20 * 5), nrow = 20)  # 20 genes × 5 samples
#' trained.model <- list(
#'   b0 = matrix(rnorm(40), nrow = 20),
#'   b1 = matrix(runif(40, 0.5, 2), nrow = 20),
#'   Beta = matrix(rnorm(20 * 2), nrow = 20),
#'   U = matrix(rnorm(20 * 5), nrow = 20),
#'   dd = runif(5, 0.1, 1),
#'   sigmaU2s = runif(20, 0.05, 0.2),
#'   sigma2s = runif(20, 0.05, 0.2)
#' )
#' X <- matrix(rnorm(2 * 5), nrow = 2)
#' Ypred <- Prediction(Ysource, X, trained.model, k.source = 1, k.target = 2)
#' }
#'
#' @seealso \code{\link{GDfun}}, \code{\link{winsor}}
#'
#' @export
Prediction <- function(Ysource, X, trained.model, k.source=1, k.target=2, b1_min=0.02,
                       min.sigmaU2=0.05, verbose=FALSE){
  # Ensure X is a matrix when provided (supports single covariate)
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
  Yhat <- muk2 + GSI +SI
  ## Now just replace those genes with very small b1k by muk2
  Yhat[ids0,] <- muk2[ids0]%*%t(rep(1,N))
  ## for debugging purpose
  if (verbose) {
    return(list(Yhat=Yhat, muk2=muk2, GSI=GSI, SI=SI, ids0=ids0))
  } else {
    return(Yhat)
  }
}

#' Predict Target Values Using a Trained MatchMixeR Model
#'
#' Applies a trained multiplicative mapping model from the \code{MatchMixeR} framework to source data \code{Ysource}
#' to generate predicted target values. Each feature (row) in \code{Ysource} is transformed by a learned slope and intercept.
#'
#' @param Ysource A numeric matrix of source data. Rows correspond to features (e.g., genes), and columns to samples.
#' @param trained.model A fitted model object from \code{MatchMixeR}, which must contain a matrix \code{betamat}
#'   with named columns \code{"Intercept"} and \code{"Slope"} representing the learned coefficients per row.
#'
#' @return A numeric matrix of the same dimensions as \code{Ysource}, representing the predicted target values.
#'
#' @details
#' For each row \eqn{i}, the predicted output is computed as:
#' \deqn{Y^{\text{target}}_{i, \cdot} = \text{Intercept}_i + \text{Slope}_i \times Y^{\text{source}}_{i, \cdot}}
#' The function uses \code{sweep()} to efficiently apply row-wise multiplicative transformations.
#'
#' @examples
#' set.seed(1)
#' Ysource <- matrix(rnorm(20), nrow = 4)
#' trained.model <- list(betamat = cbind(Intercept = rnorm(4), Slope = runif(4, 0.5, 1.5)))
#' Ytarget <- predict_MM(Ysource, trained.model)
#'
#' @seealso \code{\link[stats]{lm}}
#'
#' @export
predict_MM <- function(Ysource, trained.model) {
  bb <- trained.model$betamat
  b0 <- bb[, "Intercept"]; b1 <- bb[, "Slope"]
  Ytarget <- b0 +sweep(Ysource, 1, b1, "*")
  return(Ytarget)
}
