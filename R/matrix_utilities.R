#' Robust Hat Matrix Computation via Truncated SVD
#'
#' Computes a numerically stable approximation of the projection (hat) matrix from a design matrix \code{X}
#' using truncated singular value decomposition (SVD). This method avoids instability due to near-singular
#' or collinear columns by discarding small singular values below a specified threshold.
#'
#' @param X A numeric matrix (design matrix), typically with rows representing observations and columns representing predictors.
#' @param d.prop Numeric. Relative threshold for singular value truncation. Defaults to \code{1e-6}. Threshold is computed as \code{max(sum(d) * d.prop, dmin)}.
#' @param dmin Numeric. Absolute minimum singular value threshold. Defaults to \code{1e-9}.
#'
#' @return A numeric matrix of size \code{nrow(X)} by \code{nrow(X)} representing the approximate hat (projection) matrix.
#' If all singular values fall below the threshold, returns a zero matrix.
#'
#' @details The function uses the left singular vectors \eqn{U} corresponding to retained singular values \eqn{d}
#' to compute the hat matrix as \eqn{H = U U^\top}. This approach avoids explicit inversion and stabilizes
#' under collinearity or rank deficiency.
#'
#' @examples
#' X <- matrix(rnorm(100 * 10), nrow = 100)
#' H <- rhat(X)
#' image(H)
#'
#' @seealso \code{\link[stats]{lm}}, \code{\link[stats]{svd}}
#' @export
rhat <- function(X, d.prop=1e-6, dmin=1e-9) {
  N <- nrow(X); o <- svd(X)
  ## singular value threshold. Only those eigenvalues that passed this
  ## threshold are used in the hat matrix.
  thresh <- max(sum(o$d) * d.prop, dmin)
  idx <- which(o$d>thresh)
  if (length(idx)==0) { #nothing is left
    return(matrix(0, N, N))
  } else {
    return(tcrossprod(o$u[, idx]))
  }
}

#' Robust Inversion of \eqn{X^\top X} for Tall Matrices
#'
#' Computes a regularized inverse of the Gram matrix \eqn{X^\top X} using truncated singular value decomposition (SVD).
#' Designed for use with tall, thin matrices (\eqn{n \gg p}), the method stabilizes inversion in the presence of near-collinearity
#' or small singular values, which commonly arise in regression problems.
#'
#' @param X A numeric matrix of dimension \eqn{n \times p}, where \eqn{n \gg p}. Must not be a vector or data frame.
#' @param d.prop Numeric. Proportional threshold for regularizing small singular values. Defaults to \code{1e-6}.
#' @param dmin Numeric. Minimum allowed singular value threshold. Defaults to \code{1e-9}.
#' @param dmax Numeric. Maximum allowed singular value threshold (for stability). Defaults to \code{1e9}.
#'
#' @return A numeric matrix of dimension \eqn{p \times p}, representing the robust inverse of \eqn{X^\top X}.
#' If the matrix is nearly singular, regularization ensures a stable estimate.
#'
#' @details
#' \itemize{
#'   \item If \code{X} is a row or column vector, closed-form stabilized solutions are computed directly.
#'   \item For full-rank matrices, SVD-based inversion with regularized singular values is applied.
#'   \item The function avoids division by small or zero singular values by applying a threshold capped between \code{dmin} and \code{dmax}.
#' }
#'
#' @examples
#' X <- matrix(rnorm(100 * 5), nrow = 100)
#' inv_xtx <- rsolve2(X)
#' all.equal(inv_xtx, solve(crossprod(X)), tolerance = 1e-3)  # Should be close
#'
#' @seealso \code{\link[base]{solve}}, \code{\link{rhat}}, \code{\link[stats]{svd}}
#'
#' @export
rsolve2 <- function(X, d.prop=1e-6, dmin=1e-9, dmax=1e9){
  ## A special case for vector-valued
  if (!is.matrix(X)) stop("X must be a matrix, not a dataframe or vector.")
  N <- nrow(X); p <- ncol(X)
  if (N <= p) warning("N is less or equal to p in rsolve2(); the results may be numericall unstable.")
  ## A special (very bad) case if the input is a row vector
  if (N==1){
    d <- sqrt(sum(X^2)); v <- as.vector(X)/d
    thresh <- min(max(d * d.prop, dmin), dmax)
    inv.mat <- (1/d^2 -1/thresh^2)*tcrossprod(v) +(1/thresh^2)*diag(nrow=p)
  } else if (p==1) {
    ## A special case if N>1 and the input is a column vector (good
    ## case)
    d <- sqrt(sum(X^2))
    thresh <- min(max(d * d.prop, dmin), dmax)
    inv.mat <- as.matrix(1/max(d,thresh)^2)
  } else {
    ## general case, N,p>1
    o <- svd(X)
    ## singular value threshold
    thresh <- min(max(sum(o$d) * d.prop, dmin), dmax)
    inv.mat <- o$v %*% diag(1/(pmax(o$d, thresh))^2) %*% t(o$v)
  }
  return(inv.mat)
}

#' Fast Marginal Mean over an Array Dimension
#'
#' Computes the marginal means over a specified dimension of a multi-dimensional array.
#' This is a faster alternative to \code{apply(an_array, i, mean)} when \code{i} is a single dimension.
#'
#' @param a A numeric array of arbitrary dimension.
#' @param i Integer. The dimension over which to compute the marginal means.
#'
#' @return A numeric array with the \code{i}-th dimension collapsed via averaging.
#'
#' @details This function reorders the array dimensions using \code{aperm()} to bring the target dimension to the end,
#' and then efficiently applies \code{rowMeans()} to compute marginal means. It is especially efficient for large arrays
#' compared to the base \code{apply()} function.
#'
#' @examples
#' x <- array(rnorm(2 * 3 * 4), dim = c(2, 3, 4))
#' apply(x, c(1, 3), mean)  # base R version
#' ameans(x, 2)             # faster version over the second dimension
#'
#' @seealso \code{\link[base]{apply}}, \code{\link[base]{aperm}}, \code{\link[base]{rowMeans}}
#' @export
ameans <- function(a, i) {
  n <- length(dim(a))
  b <- aperm(a, c(seq_len(n)[-i], i))
  rowMeans(b, dims = n - 1)
}

#' Matrix Image Plot
#'
#' Displays a grayscale image of a matrix by flipping the vertical axis,
#' commonly used for visualizing matrices like heatmaps.
#'
#' @param mat A numeric matrix to visualize.
#' @param ... Additional arguments passed to \code{image()}.
#'
#' @return None. Produces a plot.
#' @export
mplot <- function(mat, ...) image(t(mat[nrow(mat):1,]), col=grey(seq(1, 0, length.out=101)), xaxt="n", yaxt="n", ...)

#' Trace of a Matrix
#'
#' Computes the trace of a square matrix (sum of diagonal elements).
#'
#' @param M A square numeric matrix.
#'
#' @return A numeric scalar representing the trace.
#' @export
tr <- function(M) sum(diag(as.matrix(M)))

#' Vectorize a Matrix
#'
#' Converts a matrix to a column-major vector.
#'
#' @param M A numeric matrix.
#'
#' @return A numeric vector.
#' @export
vec <- function(x, mode = "any") {
  as.vector(x, mode)
}
