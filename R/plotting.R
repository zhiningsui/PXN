#' Enhanced Boxplot with Overlaid Jittered Points
#'
#' Creates a standard boxplot with jittered individual data points overlaid for better visualization
#' of sample distributions. Outlier symbols from the boxplot itself are suppressed to avoid redundancy.
#'
#' @param formula A formula specifying the grouping structure, e.g., \code{y ~ group}.
#' @param data Optional data frame containing the variables in the formula.
#' @param pch Integer or character specifying the plotting symbol for the jittered points. Defaults to \code{19} (solid circles).
#' @param ... Additional arguments passed to \code{boxplot()}.
#'
#' @return An object of class \code{"boxplot"}, invisibly returned (the same as \code{boxplot()}).
#'
#' @details
#' \itemize{
#'   \item Boxplot outliers are suppressed (\code{outpch = NA}) to focus on the jittered points.
#'   \item Jittering is applied using \code{method = "jitter"} to prevent overplotting of identical or similar values.
#'   \item Useful for small- to medium-sized datasets where individual points are informative.
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' group <- rep(c("A", "B"), each = 20)
#' value <- c(rnorm(20, mean = 0), rnorm(20, mean = 1))
#' df <- data.frame(group, value)
#' myboxplot(value ~ group, data = df)
#' }
#'
#' @seealso \code{\link[graphics]{boxplot}}, \code{\link[graphics]{stripchart}}
#'
#' @importFrom graphics boxplot stripchart
#' @export
myboxplot <- function(formula, data=NULL, pch=19, ...) {
  bx <- boxplot(formula, outpch = NA, data=data, ...)
  stripchart(formula, vertical=TRUE, method="jitter", pch=pch, add=TRUE, data=data)
  return(bx)
}
