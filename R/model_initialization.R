
#' Initial Parameter Estimation for Cross-Platform Normalization
#'
#' Provides initial estimates of model parameters for cross-platform normalization or mapping models.
#' Handles genes with both sufficient and insufficient variance separately. Acts as a wrapper around
#' \code{InitEst.LargeSTD()} for genes with sufficiently large variance, and fills simple initializations for others.
#'
#' @param X Optional covariate matrix. If no covariates are used, set \code{X = NULL}.
#' @param Y A numeric matrix of stacked domain-specific observed outcomes, with \code{m*K} rows and \code{N} columns.
#' @param K Integer. Number of platforms (domains).
#' @param L Integer. Number of latent factors to extract for the latent structure (\code{U}, \code{dd}).
#' @param min.sigma2 Numeric. Minimum variance threshold for retaining genes. Genes with lower variance are handled separately. Defaults to \code{0.01}.
#' @param min.sigmaU2 Numeric. Minimum threshold for latent variance estimates. Defaults to \code{0.05}.
#' @param p.correct Logical. If \code{TRUE}, performs variance proportion correction for latent structure. Defaults to \code{TRUE}.
#'
#' @return A list containing the initial parameter estimates:
#' \itemize{
#'   \item \code{b0}: Gene-by-platform matrix of intercepts.
#'   \item \code{b1}: Gene-by-platform matrix of slopes.
#'   \item \code{Beta}: Covariate coefficient matrix (if \code{X} is provided).
#'   \item \code{U}: Estimated eigenvector matrix for latent structure.
#'   \item \code{dd}: Estimated eigenvalues (squared singular values).
#'   \item \code{varprops}: Cumulative variance proportions explained by latent factors.
#'   \item \code{sigmaU2s}: Gene-specific latent variances.
#'   \item \code{sigma2s}: Gene-specific residual variances.
#'   \item \code{N}: Number of samples.
#'   \item \code{K}: Number of platforms.
#'   \item \code{low.var.genes}: Indices of genes removed due to low variance.
#' }
#'
#' @details
#' \itemize{
#'   \item Genes with low variance (\code{rowVars(Y)} below \code{min.sigma2} across all platforms) are assigned simple defaults.
#'   \item For genes with sufficient variance, full initialization is performed using \code{InitEst.LargeSTD()}.
#'   \item Covariate matrix \code{X} is automatically coerced to a matrix if necessary.
#'   \item If latent factors are not identified (e.g., \code{L = 0}), the returned \code{U} will be \code{NULL}.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' Y <- matrix(rnorm(1000), nrow = 20)
#' res <- InitEst(X = NULL, Y = Y, K = 4, L = 3)
#' names(res)
#' }
#'
#' @seealso \code{\link{InitEst.LargeSTD}}, \code{\link{BetaEst}}, \code{\link{rowVars}}
#'
#' @importFrom Rfast rowVars
#' @importFrom Rfast rowMins
#' @importFrom Rfast rowmeans
#' @importFrom Rfast rowsums
#' @export
InitEst <- function(X, Y, K, L, min.sigma2=0.01, min.sigmaU2=0.05, p.correct=TRUE){
  # In the case when there is only one covariate and X is not a matrix.
  if (!is.null(X) && !is.matrix(X)) {
    X <- matrix(X, nrow = 1)
  }

  N <- ncol(Y)
  m <- nrow(Y) / K

  ## use gene names of Y if available
  gnames <- rownames(Y)[1:m]
  ## remove genes with very small STD
  varY <- matrix(rowVars(Y), ncol = K)
  ids1 <- rowMins(varY, value = TRUE) >= min.sigma2  # genes with large STD
  low.var.genes <- which(!ids1)

  ## the main part
  # keep only high-variance genes in all K stacked submatrices
  Y2 <- Y[rep(ids1, K), , drop = FALSE]

  # Added by Zhining 02/27/2026 ---
  gnames_keep <- gnames[ids1]
  m_keep <- sum(ids1)
  # ---

  initEst1 <- InitEst.LargeSTD(X, Y2, K, L, min.sigma2=min.sigma2, min.sigmaU2=min.sigmaU2, p.correct=p.correct)
  ## initEst1 <- InitEst.LargeSTD.old(X, Y2, K, L, min.sigmaU2=min.sigmaU2)

  # Edited by Zhining 02/27/2026 ---
  # ## Combine the simple estimators for genes with very small STD and
  # ## all other genes
  # b0 <- matrix(rowmeans(Y), ncol=K); b1 <- matrix(0, m, K)
  # rownames(b0) <- rownames(b1) <- gnames
  # colnames(b0) <- colnames(b1) <- paste0("Platform", 1:K)
  # b1[ids1,] <- initEst1$b1
  # ## estimate beta matrix
  # if (is.null(X)) {
  #   Beta <- NULL
  # } else {
  #   Beta <- matrix(0, m, nrow(X))
  #   rownames(Beta) <- gnames; colnames(Beta) <- rownames(X)
  # }
  # Beta[ids1,] <- initEst1[["Beta"]]
  # ## other parameters. Note that due to collinearity, U produced by
  # ## InitEst.LargeSTD() may have less than L columns
  # U0 <- initEst1$U
  # if (is.null(U0)) { #this could happen when L==0
  #   U <- NULL
  # } else {
  #   U <- matrix(0, m, ncol(U0))
  #   rownames(U) <- gnames; colnames(U) <- colnames(U0)
  #   U[ids1,] <- initEst1$U
  # }
  # dd <- initEst1$dd
  # varprops <- initEst1$varprops
  # sigmaU2s <- rep(1, m); names(sigmaU2s) <- gnames
  # sigmaU2s[ids1] <- initEst1$sigmaU2s
  # sigma2s <- rep(0, m); names(sigma2s) <- gnames
  # sigma2s[ids1] <- initEst1$sigma2s

  ## combine estimators for retained genes only
  b0 <- matrix(rowmeans(Y2), ncol = K)
  b1 <- matrix(0, m_keep, K)
  rownames(b0) <- rownames(b1) <- gnames_keep
  colnames(b0) <- colnames(b1) <- paste0("Platform", 1:K)
  b1[,] <- initEst1$b1

  # b0 <- matrix(rowmeans(Y), ncol = K)
  # b1 <- matrix(0, m, K)
  # rownames(b0) <- rownames(b1) <- gnames
  # colnames(b0) <- colnames(b1) <- paste0("Platform", 1:K)
  # gnames_ids1 = gnames[ids1]
  # m_ids1 = sum(ids1)
  # b0 <- b0[gnames_ids1,]; b1 <- b1[gnames_ids1,]
  # b1[gnames_ids1,] <- initEst1$b1

  ## estimate beta matrix
  if (is.null(X)) {
    Beta <- NULL
  } else {
    Beta <- matrix(0, m_keep, nrow(X))
    rownames(Beta) <- gnames_keep
    colnames(Beta) <- rownames(X)
    Beta[,] <- initEst1$Beta
  }
  ## other parameters. Note that due to collinearity, U produced by
  ## InitEst.LargeSTD() may have less than L columns
  U0 <- initEst1$U
  if (is.null(U0)) { #this could happen when L==0
    U <- NULL
  } else {
    U <- matrix(0, m_keep, ncol(U0))
    rownames(U) <- gnames_keep
    colnames(U) <- colnames(U0)
    U[,] <- U0
  }
  dd <- initEst1$dd
  varprops <- initEst1$varprops

  sigmaU2s <- rep(NA, m_keep)
  names(sigmaU2s) <- gnames_keep
  sigmaU2s[] <- initEst1$sigmaU2s

  sigma2s <- rep(NA, m_keep)
  names(sigma2s) <- gnames_keep
  sigma2s[] <- initEst1$sigma2s

  # ---

  return(list(b0=b0, b1=b1, Beta=Beta, U=U, dd=dd, varprops=varprops,
              sigmaU2s=sigmaU2s, sigma2s=sigma2s, N=N, K=K, low.var.genes=low.var.genes))
}


#' Initial Parameter Estimation for Genes with Sufficient Variance
#'
#' Computes initial estimates of model parameters (\code{b0}, \code{b1}, \code{Beta}, \code{U}, \code{dd}, \code{sigmaU2s}, \code{sigma2s})
#' for genes with sufficiently large variance using a stabilized two-step procedure:
#' centering and standardization of residuals, followed by generalized probabilistic PCA (\code{gppca}) for covariance estimation.
#' This function is intended for use inside \code{InitEst()} and does not handle low-variance genes.
#'
#' @param X Optional covariate matrix. If no covariates are available, set \code{X = NULL}.
#' @param Y A numeric matrix of observed stacked outcomes, with \code{m*K} rows and \code{N} columns.
#' @param K Integer. Number of domains/platforms.
#' @param L Integer. Number of latent factors to extract.
#' @param min.sigma2 Numeric. Minimum allowed residual variance per gene. Defaults to \code{0.01}.
#' @param min.sigmaU2 Numeric. Minimum allowed latent variance per gene. Defaults to \code{0.05}.
#' @param W Numeric. Weight (between 0 and 1) splitting variance between latent (\code{sigmaU2s}) and residual (\code{sigma2s}) components. Defaults to \code{0.5}.
#' @param p.correct Logical. Whether to apply degrees of freedom correction for \code{R0sbar} after covariate regression. Defaults to \code{TRUE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{b0}: Gene-by-platform matrix of intercepts.
#'   \item \code{b1}: Gene-by-platform matrix of slopes.
#'   \item \code{Beta}: Covariate coefficient matrix (if \code{X} is provided).
#'   \item \code{U}: Estimated eigenvector matrix for shared latent structure.
#'   \item \code{dd}: Estimated eigenvalues.
#'   \item \code{varprops}: Cumulative variance proportions explained by the latent structure.
#'   \item \code{sigmaU2s}: Gene-specific latent variances.
#'   \item \code{sigma2s}: Gene-specific residual variances.
#' }
#'
#' @details
#' \itemize{
#'   \item Genes are assumed to have sufficiently large variance (standard deviation well above \code{min.sigma2}).
#'   \item Residuals are computed after optionally regressing out covariates using a robust hat matrix (\code{rhat()}).
#'   \item Standardized residuals (\code{R0s}) are averaged across domains to form \code{R0sbar}, which is used to estimate the latent structure.
#'   \item Uses \code{gppca()} to extract low-dimensional latent components.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' Y <- matrix(rnorm(1000), nrow = 20)
#' res <- InitEst.LargeSTD(X = NULL, Y = Y, K = 4, L = 3)
#' names(res)
#' }
#'
#' @seealso \code{\link{InitEst}}, \code{\link{gppca}}, \code{\link{BetaEst}}, \code{\link{Ycircbar2Cov}}
#'
#' @importFrom Rfast rowmeans
#' @export
InitEst.LargeSTD <- function(X, Y, K, L, min.sigma2=0.01, min.sigmaU2=0.05, W=0.5, p.correct=TRUE){
  N <- ncol(Y); m <- nrow(Y)/K
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
    R0 <- R0 %*% (diag(N)-HatMat)
    ## After this step, R0 is the set of residuals, which is denoted as
    ## R^{(0)} in my notes. Again, no new object is created to save
    ## memory
  } else {
    p <- 0; Xnames <- NULL
  }
  ## the initial estimate of b1 is simply sample STDs
  s0 <- sqrt(rowsums(R0*R0)/(N-1))
  b1 <- matrix(s0, nrow=m); colnames(b1) <- paste0("Platform", 1:K)
  ## standardize R0 --> R0s (an approximation of Ycirc)
  R0s <- R0/s0
  R0sbar <- ameans(array(R0s, c(m,K,N)), 2)
  ## added 04/10/2025: To correct for the lost DFs due to regression
  if (p.correct) R0sbar <- sqrt((N-1)/(N-1-p))*R0sbar
  CovEsts <- Ycircbar2Cov(R0sbar, K, L, min.sigma2=min.sigma2, min.sigmaU2=min.sigmaU2, W=W)
  U <- CovEsts$U; dd <- CovEsts$dd
  sigmaU2s=CovEsts$sigmaU2s; sigma2s=CovEsts$sigma2s
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
  ## 08/25/2021 Return b0 and b1 in matrix format
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
  return(list(b0=b0, b1=b1, Beta=Betahat, U=U, dd=dd,
              varprops=CovEsts$varprops, sigmaU2s=sigmaU2s,
              sigma2s=sigma2s))
}

#' Covariance Parameter Estimation from Averaged Residuals
#'
#' Estimates the latent covariance structure (\code{U}, \code{dd}) and variance components (\code{sigmaU2s}, \code{sigma2s})
#' from the matrix of domain-averaged residuals \code{Ycircbar}. Uses generalized probabilistic PCA (\code{gppca()})
#' to separate shared structure and individual noise contributions.
#'
#' @param Ycircbar A numeric matrix of residuals averaged across domains, with rows as genes/features and columns as samples.
#' @param K Integer. Number of domains/platforms.
#' @param L Integer. Number of latent factors to extract.
#' @param min.sigma2 Numeric. Minimum allowed value for residual variances. Defaults to \code{0.01}.
#' @param min.sigmaU2 Numeric. Minimum allowed value for latent variances. Defaults to \code{0.05}.
#' @param W Numeric. Weight parameter (between 0 and 1) controlling the division of variance between \code{sigmaU2s} and \code{sigma2s}. Defaults to \code{0.5}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{U}: Estimated eigenvector matrix for shared latent structure.
#'   \item \code{dd}: Corresponding eigenvalues (latent variances).
#'   \item \code{varprops}: Cumulative variance proportions explained by the latent factors.
#'   \item \code{sigmaU2s}: Estimated gene-specific latent variances.
#'   \item \code{sigma2s}: Estimated gene-specific residual variances.
#' }
#'
#' @details
#' \itemize{
#'   \item If \code{K = 1} (no platform pairing), \code{Ycircbar} is treated as standardized and \code{scale = TRUE} is used in \code{gppca()}.
#'   \item If \code{K} \eqn{\geq} \code{2}, platform-paired adjustment is performed, and \code{scale = FALSE} is used.
#'   \item Variance decomposition ensures stability by enforcing minimum thresholds \code{min.sigma2} and \code{min.sigmaU2}.
#'   \item If no latent structure is detected (\code{U = NULL}), a simplified variance model is returned.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' Ycircbar <- matrix(rnorm(20 * 50), nrow = 20)
#' cov.params <- Ycircbar2Cov(Ycircbar, K = 2, L = 3)
#' names(cov.params)
#' }
#'
#' @seealso \code{\link{gppca}}
#'
#' @export
Ycircbar2Cov <- function(Ycircbar, K, L, min.sigma2=0.01, min.sigmaU2=0.05, W=0.5) {
  N <- ncol(Ycircbar); m <- nrow(Ycircbar)
  if (K==1) { #no pairing
    ## scale=TRUE --> \Sigma_{Y} is a correlation matrix
    gg <- gppca(t(Ycircbar), scale=TRUE, retx=FALSE, L=L, min.sigma2=min.sigma2+min.sigmaU2); U <- gg$U; dd <- gg$dd
    sigmaU2s <- gg$sigma2s*W; sigma2s <- gg$sigma2s*(1-W)
  } else { #with K>=2 paired platforms
    ## scale=FALSE --> \Sigma_{Y} is NOT a correlation matrix
    gg <- gppca(t(Ycircbar), scale=FALSE, retx=FALSE, L=L, min.sigma2=min.sigma2+min.sigmaU2); U <- gg$U; dd <- gg$dd
    if (is.null(U)) { #no AA'
      varR0bar <- gg$sigma2s
    } else {
      diagAA <- drop((U*U)%*%dd); varR0bar <- diagAA+gg$sigma2s
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


#' Compute the Covariance Matrix \eqn{\Sigma_Y} from Estimated Parameters
#'
#' Computes the full covariance matrix or only the diagonal (variance terms) of \eqn{\Sigma_Y}
#' based on the estimated model parameters from PPCA-Xnorm or generalized PPCA models.
#'
#' @param est A fitted model object (e.g., output from \code{InitEst()}, \code{gppca()}, or related functions)
#'   containing at least \code{U}, \code{dd}, \code{sigmaU2s}, and \code{sigma2s}.
#' @param var.only Logical. If \code{TRUE}, returns only the vector of marginal variances (diagonal elements of \eqn{\Sigma_Y}). Defaults to \code{FALSE}.
#'
#' @return
#' \itemize{
#'   \item If \code{var.only = TRUE}, returns a numeric vector of variances for each gene/feature.
#'   \item If \code{var.only = FALSE}, returns a full covariance matrix \eqn{\Sigma_Y}.
#' }
#'
#' @details
#' \itemize{
#'   \item \eqn{\Sigma_Y = \text{diag}(\sigma^2 + \sigma_U^2) + AA'}, where \eqn{AA'} captures the shared latent structure.
#'   \item If no latent factors are available (\code{U = NULL} or \code{dd = NULL}), only the diagonal variance terms are used.
#'   \item This function automatically adapts to simple cases (e.g., pure noise models) and complex cases (low-rank plus diagonal structure).
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' est <- list(
#'   U = matrix(rnorm(100), nrow = 10),
#'   dd = runif(10, 0.1, 1),
#'   sigmaU2s = runif(10, 0.05, 0.1),
#'   sigma2s = runif(10, 0.05, 0.1)
#' )
#' SigmaY_full <- getSigmaY(est)
#' SigmaY_diag <- getSigmaY(est, var.only = TRUE)
#' }
#'
#' @seealso \code{\link{InitEst}}, \code{\link{gppca}}
#'
#' @export
getSigmaY <- function(est, var.only=FALSE) {
  U <- est$U; dd <- est$dd; sigmaU2s <- est$sigmaU2s; sigma2s <- est$sigma2s
  ## to make it work with the results produced by gppca()
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
    if (is.null(dd)) { #L=0
      AA <- 0
    } else { #L>0
      AA <- U%*%diag(x=dd,nrow=length(dd))%*%t(U)
    }
    return(diag(x=sigmaU2s+sigma2s)+AA)
  }
}

