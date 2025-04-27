#' Median Imputation and Winsorization by 5-Standard-Deviation Rule
#'
#' Performs robust row-wise preprocessing on numeric matrix-like data by first imputing missing values using row-wise medians,
#' and then winsorizing outliers based on a trimmed mean and standard deviation threshold. Useful for preprocessing expression
#' data, intensity profiles, or any other high-dimensional numeric matrix.
#'
#' @param mydata A numeric matrix or data frame. Each row is treated as a feature (e.g., gene), and columns are samples.
#' @param n.std Numeric. The number of standard deviations from the trimmed mean to use for winsorizing. Defaults to \code{5}.
#' @param trim Numeric. The trimming proportion (between 0 and 0.5) used to compute the trimmed mean and standard deviation. Defaults to \code{0.1}.
#' @param verbose Logical. Whether to print summary information about imputation and winsorization. Defaults to \code{FALSE}.
#'
#' @return A numeric matrix of the same dimensions as \code{mydata}, with missing values imputed and outliers winsorized.
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
winsor <- function(mydata, n.std=5, trim=0.1, verbose=FALSE) {
  mydata <- as.matrix(mydata)
  ## median imputation of missing values for convenience
  n.na <- sum(is.na(mydata))
  if (n.na>0) {
    na.rows <- which(rowSums(is.na(mydata))>0)
    for (i in na.rows) mydata[i,is.na(mydata[i,])] <- median(mydata[i,], na.rm=TRUE)
  }
  ## pre-compute trimmed mean/stds and then lower/upper bounds
  ysort <- rowSort(mydata) #Rfast package
  n <- ncol(mydata); nL <- ceiling(n*trim); nU <- floor(n*(1-trim))
  ytrimmed <- ysort[, seq(nL, nU)]
  mus <- rowMeans(ytrimmed); sds <- rowVars(ytrimmed, std=TRUE)
  Ls <- mus-n.std*sds; Us <- mus+n.std*sds
  ## Winsorization
  y <- t(sapply(1:nrow(mydata), function(i) {
    x <- mydata[i,]; L <- Ls[i]; U <- Us[i]
    normalRange <- range(x[x>=L & x<=U])
    x[x<L] <- normalRange[1]; x[x>U] <- normalRange[2]
    return(x)
  }))
  rownames(y) <- rownames(mydata)
  if (verbose) {
    n.out <- sum(y-mydata != 0)
    N <- nrow(mydata)*ncol(mydata)
    print(paste0(n.na, " (", round(100*n.na/N,2), "%) NAs were replaced by gene-specific trimmed means. ",
                 n.out, " (", round(100*n.out/N, 2), "%) outliers were replaced by trimmed mean +/-", n.std, "*STDs." ))
  }
  return(y)
}
