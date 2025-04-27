#' Linear Modeling with Empirical Bayes Moderation (limma Wrapper)
#'
#' Fits gene-wise linear models using covariates \code{v} and applies empirical Bayes moderation
#' to stabilize variance estimates. Returns estimated coefficients, t-statistics, p-values, and adjusted p-values
#' for each covariate of interest.
#'
#' @param gdata A numeric gene expression matrix. Rows are genes/features, columns are samples.
#' @param v A covariate matrix or vector, where rows correspond to samples and columns to covariates.
#' @param padj.method Character string specifying the method for p-value adjustment. Defaults to \code{"BH"} (Benjamini-Hochberg).
#'
#' @return A \code{data.frame} containing:
#' \itemize{
#'   \item Estimated coefficients (\code{betahat.*}) for each covariate.
#'   \item Moderated t-statistics (\code{tstat.*}) for each covariate.
#'   \item Raw p-values (\code{pvals.*}) for each covariate.
#'   \item Adjusted p-values (\code{adjP.*}) for each covariate.
#' }
#'
#' @details
#' \itemize{
#'   \item Samples with missing values in \code{v} are removed prior to model fitting.
#'   \item The model fitted is: \eqn{Y = \text{Intercept} + v_1 + v_2 + \cdots + \epsilon}.
#'   \item \code{eBayes()} from the \pkg{limma} package is used to perform empirical Bayes shrinkage of variance estimates.
#' }
#'
#' @examples
#' \dontrun{
#' library(limma)
#' set.seed(123)
#' gdata <- matrix(rnorm(1000), nrow = 100)  # 100 genes × 10 samples
#' v <- matrix(rnorm(10 * 2), ncol = 2)       # 2 covariates
#' colnames(v) <- c("Treatment", "Batch")
#' res <- limma(gdata, v)
#' head(res)
#' }
#'
#' @seealso \code{\link[limma]{lmFit}}, \code{\link[limma]{eBayes}}
#'
#' @importFrom limma lmFit eBayes
#' @importFrom stats p.adjust
#' @export
limma <- function(gdata, v, padj.method="BH"){
  v <- as.matrix(v); vn <- colnames(v)
  if (is.null(vn)) {
    vn <- paste0("X", 1:ncol(v)); colnames(v) <- vn
  }
  ## remove missing samples
  na.id <- apply(v, 1, function(x) any(is.na(x)))
  gdata <- gdata[, !na.id]; v <- v[!na.id,,drop=FALSE]
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


#' Row-wise Linear Modeling (Affine-Invariant Alternative to \code{limma})
#'
#' Fits ordinary least squares (OLS) linear models row-by-row across high-throughput data (e.g., gene expression),
#' manually computing coefficients, standard errors, t-statistics, p-values, and adjusted p-values.
#' Unlike \code{limma()}, this implementation is strictly invariant and equivariant to affine transformations of \code{Y}.
#'
#' @param gdata A numeric matrix of observations. Rows are features (e.g., genes), columns are samples.
#' @param v A covariate matrix or vector. Rows correspond to samples, columns correspond to covariates.
#' @param padj.method Character string specifying the method for p-value adjustment. Defaults to \code{"BH"} (Benjamini-Hochberg).
#'
#' @return A \code{data.frame} containing:
#' \itemize{
#'   \item Estimated coefficients (\code{betahat.*}) for each covariate.
#'   \item Raw t-statistics (\code{tstat.*}) for each covariate.
#'   \item Raw two-sided p-values (\code{pvals.*}) for each covariate.
#'   \item Adjusted p-values (\code{adjP.*}) for each covariate.
#' }
#'
#' @details
#' \itemize{
#'   \item Samples with missing covariate values are removed before model fitting.
#'   \item An intercept is automatically included in the design matrix.
#'   \item Manual matrix-based calculation ensures strict affine invariance of the fitted statistics.
#'   \item P-values are based on the classical t-distribution with \code{n - p} degrees of freedom.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' gdata <- matrix(rnorm(1000), nrow = 100)  # 100 genes × 10 samples
#' v <- matrix(rnorm(10 * 2), ncol = 2)       # 2 covariates
#' colnames(v) <- c("Treatment", "Batch")
#' res <- rowlm(gdata, v)
#' head(res)
#' }
#'
#' @seealso \code{\link{limma}}, \code{\link[stats]{lm}}, \code{\link[stats]{p.adjust}}
#'
#' @importFrom stats pt p.adjust
#' @export
rowlm <- function(gdata, v, padj.method="BH"){
  v <- as.matrix(v); vn <- colnames(v)
  if (is.null(vn)) {
    vn <- paste0("X", 1:ncol(v)); colnames(v) <- vn
  }
  ## remove missing samples
  na.id <- apply(v, 1, function(x) any(is.na(x)))
  if (sum(na.id)>0) {
    gdata <- gdata[, !na.id]; v <- v[!na.id,,drop=FALSE]
  }
  n <- ncol(gdata)
  ## X1 is the design matrix (with the intercept)
  X1 <- cbind(Intercept=1, v); p <- ncol(X1)
  ## manual calculation of betas and se(beta)
  W <- solve(crossprod(X1))
  betahat <- (W%*%t(X1))%*%t(gdata)
  sigma2s <- colSums( (t(gdata)-X1%*%betahat)^2) / (n-p)
  se.betahat <- sqrt(diag(W)%*%t(sigma2s))
  ## we don't care about the intercept
  tstats <- t(betahat/se.betahat)[, -1, drop=FALSE]
  colnames(tstats) <- paste0("tstat.", vn)
  pvals <- 2*pt(-abs(tstats), df=n-p)
  colnames(pvals) <- paste0("pvals.", vn)
  adjP <- matrix(p.adjust(pvals, padj.method), nrow=nrow(pvals))
  colnames(adjP) <- paste0("adjP.", vn)
  return(data.frame(t(betahat), tstats, pvals, adjP))
}


#' Fisher's Method for Combining P-values
#'
#' Combines multiple independent p-values into a single global p-value
#' using Fisher's combination test.
#'
#' @param ps A numeric vector of individual p-values to combine.
#' @param pmin Numeric. Minimum threshold for individual p-values to avoid taking \code{log(0)}. Defaults to \code{1e-6}.
#'
#' @return A numeric scalar representing the combined global p-value.
#'
#' @details
#' \itemize{
#'   \item The test statistic is computed as \eqn{S = -2 \sum \log(p_i)}, where each \eqn{p_i} is thresholded by \code{pmin}.
#'   \item The combined p-value is obtained by comparing \eqn{S} to a chi-squared distribution with \code{2 × length(ps)} degrees of freedom.
#'   \item Useful for aggregating evidence across multiple tests while controlling type I error.
#' }
#'
#' @examples
#' \dontrun{
#' ps <- c(0.01, 0.05, 0.10)
#' fisher(ps)
#' }
#'
#' @seealso \code{\link[stats]{pchisq}}
#'
#' @importFrom stats pchisq
#' @export
fisher <- function(ps, pmin=1e-6) {
  S <- -2*sum(log(pmax(ps,pmin)))
  return(1-pchisq(S, 2*length(ps)))
}
