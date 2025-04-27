#' Mean Squared Error (MSE) Between Two Vectors
#'
#' Computes the mean squared error (MSE) between estimated and true values, with an option for relative MSE.
#'
#' @param est A numeric vector of estimated values.
#' @param truth A numeric vector of true values, of the same length as \code{est}.
#' @param relative Logical. If \code{TRUE}, returns the relative MSE (normalized by \code{sum(truth^2)}). Defaults to \code{FALSE}.
#'
#' @return A numeric scalar representing the MSE or relative MSE.
#'
#' @details
#' \itemize{
#'   \item If \code{relative = FALSE}, computes \eqn{\frac{1}{n} \sum_i (\hat{y}_i - y_i)^2}.
#'   \item If \code{relative = TRUE}, computes \eqn{\frac{\sum_i (\hat{y}_i - y_i)^2}{\sum_i y_i^2}}, which is scale-invariant.
#'   \item Missing values (\code{NA}) are ignored in the sum.
#' }
#'
#' @examples
#' set.seed(42)
#' truth <- rnorm(100)
#' est <- truth + rnorm(100, sd = 0.1)
#' MSE1(est, truth)
#' MSE1(est, truth, relative = TRUE)
#'
#' @seealso \code{\link[base]{mean}}, \code{\link[metrics]{rmse}} from other packages
#' @export
MSE1  <-  function(est, truth, relative=FALSE){
  total.mse <- sum((est - truth)^2, na.rm=TRUE)
  if (relative) {
    return(total.mse/sum(truth^2, an.rm=TRUE))
  } else {
    return(total.mse/length(as.vector(truth)))
  }
}

#' Mean Squared Error (MSE) Between Estimated and Oracle \eqn{\Sigma_Y}
#'
#' Computes the mean squared error (MSE) between two covariance models \eqn{\Sigma_Y},
#' based on their low-rank and diagonal components. Designed to efficiently handle large matrices by avoiding explicit matrix reconstruction.
#'
#' @param est A fitted model object (e.g., output from \code{InitEst()} or \code{gppca()}) containing \code{U}, \code{dd}, \code{sigmaU2s}, and \code{sigma2s}.
#' @param oracle The true model object (oracle) with the same structure (\code{U}, \code{dd}, \code{sigmaU2s}, and \code{sigma2s}).
#'
#' @return A numeric scalar representing the mean squared error between the estimated and oracle \eqn{\Sigma_Y} matrices, scaled by \eqn{1/m^2}.
#'
#' @details
#' \itemize{
#'   \item The low-rank term is computed as \eqn{AA'} where \eqn{A = U \sqrt{d}}.
#'   \item The full covariance matrix \eqn{\Sigma_Y} is modeled as the sum of a low-rank term (\eqn{AA'}) and a diagonal variance term.
#'   \item The MSE is computed as the sum of squared differences between the two low-rank structures and the two sets of variances, divided by \eqn{m^2}.
#'   \item This method is memory-efficient and avoids explicitly forming large \eqn{\Sigma_Y} matrices.
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
#' oracle <- list(
#'   U = matrix(rnorm(100), nrow = 10),
#'   dd = runif(10, 0.1, 1),
#'   sigmaU2s = runif(10, 0.05, 0.1),
#'   sigma2s = runif(10, 0.05, 0.1)
#' )
#' MSE.SigmaY(est, oracle)
#' }
#'
#' @seealso \code{\link{getSigmaY}}, \code{\link{InitEst}}, \code{\link{gppca}}
#'
#' @export
MSE.SigmaY <- function(est, oracle) {
  A <- sweep(est$U, 2, sqrt(est$dd), "*"); m <- nrow(A)
  B <- sweep(oracle$U, 2, sqrt(oracle$dd), "*")
  ss.A <- est$sigmaU2s+est$sigma2s
  ss.B <- oracle$sigmaU2s+oracle$sigma2s
  ##
  Term1 <- sum(crossprod(A)^2)+sum(crossprod(B)^2)-2*sum((crossprod(A, B))^2)
  Term2 <- sum((ss.A-ss.B)^2)
  return((Term1+Term2)/m^2)
}
