#' Model Integration Across Multiple PPCA-XNORM Trained Models
#'
#' Integrates a list of separately trained PPCA-Xnorm models (typically from different datasets or study groups)
#' into a single unified model by combining intercepts, slopes, covariate effects, and latent structures.
#' Shared covariance structure is re-estimated based on the averaged residuals (\code{Ycircbar}) from all groups.
#'
#' @param modList A list of trained PPCA-Xnorm models, where each model is a list output from \code{\link{GDfun}}.
#' @param Xlist Optional list of covariate matrices corresponding to each group. If no covariates are used, set \code{Xlist = NULL}.
#' @param Ylist A list of observed outcome matrices (one for each group), stacked by platform/domain order.
#' @param PairList A list or matrix indicating the platform identities for each group, used to align and integrate \code{b0} and \code{b1}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{b0}: Gene-by-platform matrix of integrated intercepts.
#'   \item \code{b1}: Gene-by-platform matrix of integrated slopes.
#'   \item \code{Beta}: Integrated covariate coefficient matrix (if covariates were provided).
#'   \item \code{U}: Integrated eigenvector matrix of the shared latent structure.
#'   \item \code{dd}: Integrated eigenvalues (latent variances).
#'   \item \code{varprops}: Cumulative variance proportions explained by the integrated latent structure.
#'   \item \code{sigma2s}: Integrated gene-specific residual variances.
#'   \item \code{sigmaU2s}: Integrated gene-specific latent variances.
#'   \item \code{Ns}: Vector of sample sizes per group.
#'   \item \code{Ks}: Vector of number of platforms per group.
#' }
#'
#' @details
#' \itemize{
#'   \item Intercepts (\code{b0}) and slopes (\code{b1}) are combined via sample-size weighted averages across groups sharing the same platform.
#'   \item Covariate effects (\code{Beta}) are combined via weighted averages if available.
#'   \item Shared latent structure (\code{U}, \code{dd}) is re-estimated by applying SVD to concatenated \code{Ycircbar} matrices.
#'   \item Effective number of PCs (\code{L}) may decrease if singular values are numerically zero during SVD.
#' }
#'
#' @examples
#' \dontrun{
#' # Example pseudo-usage
#' set.seed(123)
#' mod1 <- GDfun(NULL, matrix(rnorm(20 * 50), nrow = 20), K = 2, params = list(b1 = rep(1, 20), U = NULL, dd = NULL, sigma2s = rep(1, 20), sigmaU2s = rep(1, 20)))
#' mod2 <- GDfun(NULL, matrix(rnorm(20 * 60), nrow = 20), K = 2, params = list(b1 = rep(1, 20), U = NULL, dd = NULL, sigma2s = rep(1, 20), sigmaU2s = rep(1, 20)))
#' modList <- list(mod1$params, mod2$params)
#' PairList <- list(c("Platform1", "Platform2"), c("Platform3", "Platform4"))
#' Ylist <- list(matrix(rnorm(20 * 50), nrow = 20), matrix(rnorm(20 * 60), nrow = 20))
#' res <- ModIntegrate(modList, Xlist = NULL, Ylist = Ylist, PairList = PairList)
#' names(res)
#' }
#'
#' @seealso \code{\link{GDfun}}, \code{\link{InitEst}}, \code{\link{fsvd}}, \code{\link{ameans}}
#'
#' @export
ModIntegrate <- function(modList, Xlist=NULL, Ylist, PairList) {
  #############################################################################
  #############################################################################
  ## Added by Zhining: Harmonize genes across groups (use shared/intersection)
  gene_list <- lapply(modList, function(m) rownames(m[["b0"]]))
  shared_genes <- Reduce(intersect, gene_list)
  if (length(shared_genes) == 0L)
    stop("No shared genes across groups.")
  if (any(lengths(gene_list) != length(shared_genes)))
    warning(sprintf("Dropping %d–%d non-shared genes to align groups.",
                    min(lengths(gene_list)) - length(shared_genes),
                    max(lengths(gene_list)) - length(shared_genes)))

  ## Reindex all model pieces by shared genes
  modList <- lapply(seq_along(modList), function(g) {
    m <- modList[[g]]
    ord <- match(shared_genes, rownames(m[["b0"]]))
    m$b0 <- m$b0[ord, , drop = FALSE]
    m$b1 <- m$b1[ord, , drop = FALSE]
    if (!is.null(m$Beta)) m$Beta <- m$Beta[ord, , drop = FALSE]
    if (!is.null(m$U)) m$U <- m$U[ord, , drop = FALSE]
    m$sigma2s <- m$sigma2s[shared_genes]
    m$sigmaU2s <- m$sigmaU2s[shared_genes]
    m
  })

  if (!is.null(Ylist)) {
    Ylist <- lapply(seq_along(Ylist), function(g) {
      Yg <- Ylist[[g]]
      Kg <- length(PairList[[g]])
      m_old <- nrow(Yg) / Kg
      if (m_old != length(gene_list[[g]]))
        stop("Ylist[[", g, "]] row count does not match gene count × Kg.")
      arr <- array(Yg, c(m_old, Kg, ncol(Yg)))
      ord <- match(shared_genes, gene_list[[g]])
      arr2 <- arr[ord, , , drop = FALSE]
      out <- matrix(arr2, nrow = length(shared_genes) * Kg, ncol = ncol(Yg))
      colnames(out) <- colnames(Yg)
      out
    })
  }

  gnames <- shared_genes
  G  <- length(modList)
  platforms <- Reduce("union", PairList)
  nPFs <- length(platforms)
  m <- length(gnames)

  Ns0 <- sapply(modList, function(mod) mod[["N"]])
  Ks <- sapply(modList, function(mod) mod[["K"]])
  Ns <- Ns0 * Ks

  #############################################################################
  #############################################################################

  ## collect information from the inputs
  ## use genenames names of covariates if available
  # gnames <- rownames(modList[[1]][["b0"]])
  # G <- length(modList); platforms <- Reduce("union", PairList)
  # nPFs <- length(platforms); m <- length(gnames)
  # ## the effective sample size for each group is Ng*Kg
  # Ns0 <- sapply(modList, function(mod) mod[["N"]])
  # Ks <- sapply(modList, function(mod) mod[["K"]])
  # ## Ns are the effective samples in each platform
  # Ns <- Ns0*Ks

  ## extract the estimates from each model
  sigmaU2mat <- sapply(modList, function(mod) mod[["sigmaU2s"]])
  sigma2mat <- sapply(modList, function(mod) mod[["sigma2s"]])
  Betas <- lapply(modList, function(mod) mod[["Beta"]])
  b0s <- lapply(modList, function(mod) mod[["b0"]])
  b1s <- lapply(modList, function(mod) mod[["b1"]])
  Us <- lapply(modList, function(mod) mod[["U"]])
  dds <- lapply(modList, function(mod) mod[["dd"]])

  ## sanity check
  Ls <- sapply(dds, length); L <- min(Ls)
  if (diff(range(Ls))!=0) {
    warning(paste0("Number of PCs (L) are not equal across all groups. The integrated model will use the minimum number of PCs, Lmin=", L))
  }
  ## use weighted mean for sigma2s, sigmaU2s, and Beta
  ww <- Ns/sum(Ns)
  sigma2s <- drop(sigma2mat%*%ww)
  sigmaU2s <- drop(sigmaU2mat%*%ww)
  ## dd
  ddmat <- sapply(dds, function(ds) ds[1:L])
  dd <- drop(ddmat%*%ww)
  ## Beta
  if (any(sapply(Betas, is.null))) {
    Betahat <- NULL
  } else {
    Betahat <- Reduce("+", lapply(1:G, function(g) ww[g]*Betas[[g]]))
  }
  ## b0 and b1 needs a bit workx
  b0 <- b1 <- matrix(0, m, nPFs)
  colnames(b0) <- colnames(b1) <- platforms
  for (plfm in platforms){
    ids <- lapply(PairList, function(x) which(x==plfm))
    ## sanity check
    LL <- sapply(ids, length)
    if (max(LL>1)) stop("Please check PairList; one platform can only appear once in each group.")
    ## now the main algorithm
    Grps <- which(LL==1)
    if (length(Grps)==1){ #no need to do weighted average
      WithinGrpID <- ids[[Grps]]
      b0[, plfm] <- b0s[[Grps]][((WithinGrpID-1)*m+1):(WithinGrpID*m)]
      b1[, plfm] <- b1s[[Grps]][((WithinGrpID-1)*m+1):(WithinGrpID*m)]
    } else {
      b0[, plfm] <- Reduce("+", lapply(Grps, function(g) {
        WithinGrpID <- ids[[g]]
        return(Ns[g]*b0s[[g]][((WithinGrpID-1)*m+1):(WithinGrpID*m)])
      }))/sum(Ns[Grps])
      b1[, plfm] <- Reduce("+", lapply(Grps, function(g) {
        WithinGrpID <- ids[[g]]
        return(Ns[g]*b1s[[g]][((WithinGrpID-1)*m+1):(WithinGrpID*m)])
      }))/sum(Ns[Grps])
    }
  }
  ## Estimate U, modified from Ycircbar2Cov()
  if (L==0) {                         #do not compute U and dd
    U <- dd <- varprops <- NULL
  } else {
    ## use the combined Ycircbars to estimate U and dd
    Ycircbars <- Reduce("cbind", lapply(1:G, function(g) {
      Yg <- Ylist[[g]]
      platforms.g <- PairList[[g]]; Kg <- length(platforms.g)
      b0g <- as.vector(b0[, platforms.g])
      b1g <- as.vector(b1[, platforms.g])
      Ys <- (Yg-b0g)/b1g
      Ytilde <- ameans(array(Ys, c(m,Kg,ncol(Ys))), 2)
      if (is.null(Xlist)){
        Ycirc <- Ys
      } else {
        Xg <- Xlist[[g]]
        if (is.null(Xg)) {
          Ycirc <- Ys
        } else {
          BetaX <- Ytilde %*% rhat(t(Xg))
          Ycirc <- Ys -rep(1,Kg)%x%BetaX
        }
      }
      Ycircbar <- ameans(array(Ycirc, c(m,Kg,ncol(Ycirc))), 2)
    }))
    ## Using SVD of Ycircbar to estimate U
    ss <- fsvd(Ycircbars, k=L); U <- ss$u
    ## varprop is very rough
    varprops <- cumsum(ss$d^2)/sum(Ycircbars^2)
    ## there is a chance that ncol(U) < L, because some singular
    ## values may be numerically zero (especially when L is large)
    if (ncol(U) < L) {
      warning(paste0("Effective rank of Ycircbar is less than Lmin=", L, ", so we reduce the number of PCs to L=", ncol(U), " to ensure the numerical stability of the results."))
      L <- ncol(U)
    }
  }
  ## Final step: add row/column names to the output. Note that both
  ## b0 and b1 are matrices, not long vectors
  rownames(b0) <- rownames(b1) <- rownames(U)  <- names(sigmaU2s)  <- names(sigma2s) <- gnames
  if (!is.null(Betahat)) {
    rownames(Betahat) <- gnames; colnames(Betahat) <- rownames(Xs[[1]])
  }
  colnames(U) <- paste0("PC", 1:ncol(U))
  return(list(b0=b0, b1=b1, Beta=Betahat, U=U, dd=dd, varprops=varprops,
              sigma2s=sigma2s, sigmaU2s=sigmaU2s, Ns=Ns0, Ks=Ks))
}

