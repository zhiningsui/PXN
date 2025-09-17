#' Simulation expression matrices for PXN
#'
#' @description
#' `sim_expres` contains simulated gene-expression matrices used in the
#' PPCAXNORM vignette and examples. Each element is a matrix whose rows are
#' genes and columns are samples, with platforms stacked by rows.
#'
#' @format
#' A named list of length 3:
#' \describe{
#'   \item{AB}{Numeric matrix of dimension \eqn{2m \times n}.
#'             Stacked as \code{rbind(YA, YB)} (top = A, bottom = B).}
#'   \item{BCD}{Numeric matrix of dimension \eqn{3m \times n}.
#'              Stacked as \code{rbind(YB, YC, YD)}.}
#'   \item{EF}{Numeric matrix of dimension \eqn{2m \times n}.
#'             Stacked as \code{rbind(YE, YF)}.}
#' }
#'
#' @details
#' Within each element, all platforms share the same gene set and row order.
#' Columns (samples) are specific to that dataset and align with the matching
#' covariate matrix in \code{\link{sim_covars}} (i.e., \code{ncol(sim_covars$AB)
#' == ncol(sim_expres$AB)}, etc.).
#'
#' @usage data(sim_expres)
#'
#' @examples
#' data(sim_expres)
#' str(sim_expres, max.level = 1)
#' dim(sim_expres$AB)
#'
#' ## Split AB into source (A) and target (B) blocks
#' m <- nrow(sim_expres$AB) / 2
#' YA <- sim_expres$AB[1:m, ]
#' YB <- sim_expres$AB[(m + 1):(2 * m), ]
#'
#' @seealso \code{\link{sim_covars}},
#'   \code{\link{InitEst}}, \code{\link{GDfun}},
#'   \code{\link{Prediction}}, \code{\link{ModIntegrate}}
"sim_expres"

#' Simulation covariate matrices aligned to \code{sim_expres}
#'
#' @description
#' `sim_covars` contains covariate matrices that column-align with the
#' expression matrices in \code{\link{sim_expres}}. Each element is a small
#' \eqn{p \times n} numeric matrix (in the shipped data \eqn{p = 3}) whose rows
#' are covariates (see \code{rownames()}) and whose columns are samples in the
#' matching expression matrix.
#'
#' @format
#' A named list of length 3:
#' \describe{
#'   \item{AB}{Numeric matrix \eqn{p \times n}; columns match \code{sim_expres$AB}.}
#'   \item{BCD}{Numeric matrix \eqn{p \times n}; columns match \code{sim_expres$BCD}.}
#'   \item{EF}{Numeric matrix \eqn{p \times n}; columns match \code{sim_expres$EF}.}
#' }
#'
#' @details
#' The first row is often an intercept and remaining rows are simulated
#' covariates (e.g., clinical or technical variables). Use
#' \code{rownames(sim_covars$AB)} to see the covariate labels. Before model
#' fitting, ensure there are no missing values and that the column order matches
#' the corresponding \code{sim_expres} matrix.
#'
#' @usage data(sim_covars)
#'
#' @examples
#' data(sim_expres); data(sim_covars)
#' stopifnot(ncol(sim_covars$AB) == ncol(sim_expres$AB))
#' rownames(sim_covars$AB)        # covariate names
#'
#' ## Minimal dimensions check for all sets
#' lapply(sim_covars, dim)
#'
#' @seealso \code{\link{sim_expres}},
#'   \code{\link{InitEst}}, \code{\link{GDfun}},
#'   \code{\link{Prediction}}, \code{\link{ModIntegrate}}
"sim_covars"

#' Median Imputation and Winsorization by 5-Standard-Deviation Rule
#'
#' Performs robust row-wise preprocessing on numeric matrix-like data by first imputing missing values using row-wise medians,
#' and then winsorizing outliers based on a trimmed mean and standard deviation threshold. Useful for preprocessing expression
#' data, intensity profiles, or any other high-dimensional numeric matrix.
#'
#' @param data A numeric matrix or data frame. Each row is treated as a feature (e.g., gene), and columns are samples.
#' @param n.std Numeric. The number of standard deviations from the trimmed mean to use for winsorizing. Defaults to \code{5}.
#' @param trim Numeric. The trimming proportion (between 0 and 0.5) used to compute the trimmed mean and standard deviation. Defaults to \code{0.1}.
#' @param verbose Logical. Whether to print summary information about imputation and winsorization. Defaults to \code{FALSE}.
#'
#' @return A numeric matrix of the same dimensions as \code{data}, with missing values imputed and outliers winsorized.
#'
#' @details
#' \itemize{
#'   \item Missing values are imputed using the row-wise median (ignoring \code{NA}s).
#'   \item For each row, a trimmed mean and standard deviation are computed using central values (controlled by \code{trim}).
#'   \item Values outside the \code{mean ± n.std × sd} range are winsorized to the nearest in-range value.
#'   \item Uses fast row-wise operations from the \pkg{Rfast} package.
#' }
#'
#' @examples
#' set.seed(42)
#' X <- matrix(rnorm(1000), nrow = 20)
#' X[sample(length(X), 10)] <- NA  # introduce some NAs
#' X[1, 1] <- 1000  # introduce an outlier
#' Y <- winsor(X, verbose = TRUE)
#'
#' @seealso \code{\link[Rfast]{rowSort}}, \code{\link[Rfast]{rowVars}}
#'
#' @importFrom Rfast rowSort rowVars
#' @export
winsor <- function(data, n.std = 5, trim = 0.1, verbose = FALSE) {
  data <- as.matrix(data)
  n.na <- sum(is.na(data))
  if (n.na>0) {
    na.rows <- which(rowSums(is.na(data))>0)
    for (i in na.rows) data[i,is.na(data[i,])] <- median(data[i,], na.rm=TRUE)
  }
  ## pre-compute trimmed mean/stds and then lower/upper bounds
  ysort <- rowSort(data) #Rfast package
  n <- ncol(data); nL <- ceiling(n*trim); nU <- floor(n*(1-trim))
  ytrimmed <- ysort[, seq(nL, nU)]
  mus <- rowMeans(ytrimmed); sds <- rowVars(ytrimmed, std=TRUE)
  Ls <- mus-n.std*sds; Us <- mus+n.std*sds
  ## Winsorization
  y <- t(sapply(1:nrow(data), function(i) {
    x <- data[i,]; L <- Ls[i]; U <- Us[i]
    normalRange <- range(x[x>=L & x<=U])
    x[x<L] <- normalRange[1]; x[x>U] <- normalRange[2]
    return(x)
  }))
  rownames(y) <- rownames(data)
  if (verbose) {
    n.out <- sum(y-data != 0)
    N <- nrow(data)*ncol(data)
    print(paste0(n.na, " (", round(100*n.na/N,2), "%) NAs were replaced by gene-specific trimmed means. ",
                 n.out, " (", round(100*n.out/N, 2), "%) outliers were replaced by trimmed mean +/-", n.std, "*STDs." ))
  }
  return(y)
}
