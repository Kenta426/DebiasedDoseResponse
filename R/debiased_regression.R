.pseudo.outcomes <- function(y, a, w, mu, g){
  n <- length(a)
  W.new <- data.frame(w[rep(1:n, n), ])  # repeat a matrix over rows
  w <- data.frame(w); colnames(w) <- colnames(W.new)
  A.new <- rep(a, each=n) # repeat
  muhat.mat <- matrix(mu(A.new, W.new), byrow=TRUE, nrow=n, ncol=n)
  muhat.obs <- diag(muhat.mat)
  mhat.obs <- rowMeans(muhat.mat)
  ghat.obs <- g(a, w)
  list(pseudo.outcome=(y - muhat.obs)/(ghat.obs) + mhat.obs,
       muhat.mat=muhat.mat,
       mhat.obs=mhat.obs)
}


#' Compute debiased nonparametric inference for a covariate-adjusted regression
#' function.
#'
#' This function performs nonparametric inference on the causal dose-response
#' curve associated with continuous exposure. It performs both
#' pointwise and uniform inference. The resulting confidence intervals and
#' confidence bands attain asymptotically valid coverage without undersmoothing
#' subject to technical conditions.
#' See the accompanying paper for details.
#'
#' @param Y \code{n x 1} numeric vector of observed outcome values.
#' @param A \code{n x 1} numeric vector of observed exposure values.
#' @param W \code{n x p} data.frame of potential confounders.
#' @param mu outcome regression function (takes A and W and returns R)
#' @param g standardized propensity function (takes A and W and returns R)
#' @param tau the ratio of two parameters h/b
#' @param eval.pts the list of points to run inference on. the default is
#' 30 equi-distant points from the 5th and 95th percentile of the exposure values
#' @param kernel.type the type of kernel functions. The options are "uni" (Uniform),
#' "tri" (Triangular), "gau" (Gaussian), and "epa" (epanechnikov).
#' The default is "epa".
#' @param bw.seq the choice of bandwidths. suppressed unless leave-one-out
#' cross-validation is selected.
#'
#' @return
#'
#'
#' @examples TODO
#' @export

debiased_inference <- function(Y, A, W, mu, g, tau=1, eval.pts=NULL,
                               kernel.type="epa", bw.seq=NULL){
  require(nprobust)
  # compute an estimated pseudo-outcome sequence
  ord <- order(A)
  Y <- Y[ord]; W <- W[ord,]; A <- A[ord]
  n <- length(A)
  pseudo.res <- .pseudo.outcomes(Y, A, W, mu, g)
  pseudo.out <- pseudo.res$pseudo.outcome
  muhat.mat <- pseudo.res$muhat.mat
  mhat.obs <- pseudo.res$mhat.obs

  kern <- function(t){.kern(t, kernel=kernel.type)}

  if (is.null(eval.pts)){
    eval.pts <- seq(quantile(A, 0.05), quantile(A, 0.95),  length.out=30)
  }

  # compute bandwidth
  h.dpi <- nprobust::lpbwselect(pseudo.out, A, eval=eval.pts,
                                bwselect="imse-dpi")$bws[,2]
  b.dpi <- h.dpi/tau
  # fit regression
  est.res <- .lprobust(A, pseudo.out, h.dpi, b.dpi,
                       eval=eval.pts, kernel.type="epa")

  # compute an estimated influence function sequence and confidence intervals
  rinf.fns <- mapply(function(a, h.val, b.val){
    .compute.rinfl.func(pseudo.out, A, a, h.val, b.val, kern,
                        muhat.mat, mhat.obs)
  }, eval.pts, h.dpi, b.dpi)
  rif.se <- apply(rinf.fns, 2, sd)
  ci.ll.p <- est.res[,"theta.hat"] - 1.96 * rif.se / sqrt(n)
  ci.ul.p <- est.res[,"theta.hat"] + 1.96 * rif.se / sqrt(n)

  # compute uniform bands
  get.unif.ep <- function(alpha){
    std.inf.vals <- scale(rinf.fns)
    ep.maxes <- replicate(1e4,
                          max(abs(rbind(rnorm(n) / sqrt(n)) %*% std.inf.vals)))
    quantile(ep.maxes, alpha)
  }
  ep.unif.quant <- get.unif.ep(0.95) * rif.se / sqrt(n)
  ci.ll.u <- est.res[,"theta.hat"] - ep.unif.quant
  ci.ul.u <- est.res[,"theta.hat"] + ep.unif.quant

  # output data frame
  res <- data.frame(
    eval.pts=eval.pts,
    theta=est.res[,"theta.hat"],
    ci.ul.pts=ci.ul.p,
    ci.ll.pts=ci.ll.p,
    ci.ul.unif=ci.ul.u,
    ci.ll.unif=ci.ll.u,
    bias=est.res[,"b.hat"],
    if.val=rif.se)
  res$h <- h.dpi; res$b <- b.dpi
  return(res)
}
