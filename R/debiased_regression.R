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
#' @param tau the ratio of two parameters h/b
#' @param eval.pts the list of points to run inference on. the default is 30
#' equi-distant points from the 5th and 95th percentile of the exposure values
#' @param ... Additional control parameters. See the \strong{Details} section.
#' @details The following is a list of additional control parameters.
#'
#'\describe{
#'   \item{\code{mu}}{outcome regression function (takes A and W and returns R)}
#'   \item{\code{g}}{standardized propensity function (takes A and W and returns R)}
#'   \item{\code{kernel.type}}{The choice of a kernel function. It must be one
#'   of the following: \code{epa} (epanechnikov), \code{tri} (triangular),
#'   \code{uni} (uniform) and \code{gau} (Gaussian). The default is \code{epa}.}
#'   \item{\code{bandwidth.method}}{The choice of bandwidth selection method.
#'   Is must be one of the following: \code{DPI} (Direct Plug-in),
#'   \code{LOOCV} (Leave-one-out cross-validation), and
#'   \code{LOOCV(h=b) (Leave-one-out cross-validation with the constraint h=b)}.
#'   The default is \code{DPI}.}
#'   \item{\code{alpha.pts}}{The alpha value for the point-wise confidence
#'   intervals. The default value is 0.05.}
#'   \item{\code{unif}}{If TRUE, the function performs uniform inference.
#'   The default value is TRUE.}
#'   \item{\code{alpha.unif}}{The alpha value for the uniform confidence bands.
#'   The default value is 0.05.}
#'   \item{\code{bw.seq}}{The set of bandwidth parameters to choose from. Only
#'   used when \code{bandwidth.method} is \code{LOOCV} or \code{LOOCV(h=b)}.
#'   The default value is `seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 10)` for
#'   \code{LOOCV} and `seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 50)` for
#'   \code{LOOCV(h=b)}}
#'   \item{\code{bootstrap}}{The number of bootstrap samples used to simulate
#'   Gaussian Process in the uniform bands. The default value is 1e4.}
#'}
#' @return
#'
#' @importFrom nprobust lpbwselect
#' @import SuperLearner earth gam ranger rpart Rsolnp nnls
#' @examples
#' # Sample data
#' n <- 200; cols <- 3
#' W <- matrix(runif(n*cols), ncol=cols) # 200 * 3 matrix of covariates
#' A <- rnorm(n, mean=W%*%rnorm(cols))   # 200 * 1 vector of treatment variable
#' Y <- rnorm(n, mean=sin(A)+W%*%rnorm(cols))  # 200 * 1 vector of response variable
#' res <- debiased_inference(Y, A, W)
#' @export

debiased_inference <- function(Y, A, W, tau=1, eval.pts=NULL, ...){
  # Parse control inputs ------------------------------------------------------
  control <- .parse.debiased_inference(...)
  kernel.type <- control$kernel.type

  # Compute a nuisance estimators if not provided -----------------------------
  if(is.null(control$mu)){
    mu <- .fit.regression(Y, A, W)
  }
  else{
    mu <- control$mu
  }
  if(is.null(control$mu)){
    g <- .fit.density(A, W)
  }
  else{
    g <- control$g
  }

  # Compute an estimated pseudo-outcome sequence ------------------------------
  ord <- order(A)
  Y <- Y[ord]; W <- W[ord,]; A <- A[ord]; n <- length(A)
  pseudo.res <- .pseudo.outcomes(Y, A, W, mu, g)
  pseudo.out <- pseudo.res$pseudo.outcome
  muhat.mat <- pseudo.res$muhat.mat
  mhat.obs <- pseudo.res$mhat.obs
  kern <- function(t){.kern(t, kernel=kernel.type)}
  if (is.null(eval.pts)){
    eval.pts <- seq(quantile(A, 0.05), quantile(A, 0.95),  length.out=30)
  }

  # Compute bandwidth ---------------------------------------------------------
  if(control$bandwidth.method == "LOOCV"){
    if (is.null(control$bw.seq)){
      bw.seq <- seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 10)
    }
    else{
      bw.seq <- control$bw.seq
    }
    bw.seq.h <- rep(bw.seq, length(bw.seq))
    bw.seq.b <- rep(bw.seq, each=length(bw.seq))
    risk <- mapply(function(h, b){
      .robust.loocv(A, pseudo.out, h, b, eval.pt=eval.pts,
                    kernel.type=kernel.type)}, bw.seq.h, bw.seq.b)
    h.opt <- bw.seq.h[which.min(risk)]
    b.opt <- bw.seq.b[which.min(risk)]
  }
  else if(control$bandwidth.method == "LOOCV(h=b)"){
    if (is.null(control$bw.seq)){
      bw.seq <- seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 50)
    }
    else{
      bw.seq <- control$bw.seq
    }
    risk <- mapply(function(h, b){
      .robust.loocv(A, pseudo.out, h, b, eval.pt=eval.pts,
                    kernel.type=kernel.type)}, bw.seq, bw.seq)
    h.opt <- bw.seq[which.min(risk)]
    b.opt <- bw.seq[which.min(risk)]
  }
  else{
    h.opt <- lpbwselect(pseudo.out, A, eval=eval.pts,
                                  bwselect="imse-dpi")$bws[,2]
    b.opt <- h.opt/tau
  }

  # Fit debiased local linear regression --------------------------------------
  est.res <- .lprobust(A, pseudo.out, h.opt, b.opt,
                       eval=eval.pts, kernel.type=kernel.type)

  # Estimate influence function sequence --------------------------------------
  rinf.fns <- mapply(function(a, h.val, b.val){
    .compute.rinfl.func(pseudo.out, A, a, h.val, b.val, kern,
                        muhat.mat, mhat.obs)
  }, eval.pts, h.opt, b.opt)
  rif.se <- apply(rinf.fns, 2, sd)/sqrt(n)
  alpha <- control$alpha.pts
  ci.ll.p <- est.res[,"theta.hat"]-qnorm(1-alpha/2)*rif.se
  ci.ul.p <- est.res[,"theta.hat"]+qnorm(1-alpha/2)*rif.se

  # Compute uniform bands by simulating GP ------------------------------------
  if(control$unif){
    get.unif.ep <- function(alpha){
      std.inf.vals <- scale(rinf.fns)
      boot.samples <- control$bootstrap
      ep.maxes <- replicate(boot.samples,
                            max(abs(rbind(rnorm(n)/sqrt(n))%*%std.inf.vals)))
      quantile(ep.maxes, alpha)
    }
    alpha <- control$alpha.unif
    ep.unif.quant <- get.unif.ep(1-alpha)*rif.se
    ci.ll.u <- est.res[,"theta.hat"]-ep.unif.quant
    ci.ul.u <- est.res[,"theta.hat"]+ep.unif.quant
  }
  else{
    ci.ll.u <- ci.ul.u <- 0
  }

  # Output data.frame ---------------------------------------------------------
  res <- data.frame(
    eval.pts=eval.pts,
    theta=est.res[,"theta.hat"],
    ci.ul.pts=ci.ul.p,
    ci.ll.pts=ci.ll.p,
    ci.ul.unif=ci.ul.u,
    ci.ll.unif=ci.ll.u,
    bias=est.res[,"b.hat"],
    if.val=rif.se,
    unif.if.val=ep.unif.quant)
  res$h <- h.opt; res$b <- b.opt
  return(res)
}

#' Compute debiased nonparametric inference for the difference of
#' covariate-adjusted regression function values.
#'
#' This function performs nonparametric inference on the difference of
#' causal dose-response curve values associated with continuous exposure.
#' The resulting confidence intervals attain asymptotically valid coverage
#' without undersmoothing subject to technical conditions.
#' See the accompanying paper for details.
#'
#' @param Y \code{n x 1} numeric vector of observed outcome values.
#' @param A \code{n x 1} numeric vector of observed exposure values.
#' @param W \code{n x p} data.frame of potential confounders.
#' @param tau the ratio of two parameters h/b
#' @param eval.pts.1 the list of points to run inference on theta(a1)-theta(a2),
#' where this entry corresponds to a1. the default is 30 equi-distant points
#' from the 50th and 95th percentile of the exposure values
#' @param eval.pts.2 the list of points to run inference on theta(a1)-theta(a2),
#' where this entry corresponds to a2. the default is the 5th percentile
#' of the exposure values
#' @param ... Additional control parameters. See the \strong{Details} section.
#' @details The following is a list of additional control parameters.
#'
#'\describe{
#'   \item{\code{mu}}{outcome regression function (takes A and W and returns R)}
#'   \item{\code{g}}{standardized propensity function
#'   (takes A and W and returns R)}.
#'   \item{\code{kernel.type}}{The choice of a kernel function. It must be one
#'   of the following: \code{epa} (epanechnikov), \code{tri} (triangular),
#'   \code{uni} (uniform) and \code{gau} (Gaussian). The default is \code{epa}.}
#'   \item{\code{bandwidth.method}}{The choice of bandwidth selection method.
#'   Is must be one of the following: \code{DPI} (Direct Plug-in),
#'   \code{LOOCV} (Leave-one-out cross-validation), and
#'   \code{LOOCV(h=b) (Leave-one-out cross-validation with the constraint h=b)}.
#'   The default is \code{DPI}.}
#'   \item{\code{alpha.pts}}{The alpha value for the point-wise confidence
#'   intervals. The default value is 0.05.}
#'   \item{\code{bw.seq}}{The set of bandwidth parameters to choose from. Only
#'   used when \code{bandwidth.method} is \code{LOOCV} or \code{LOOCV(h=b)}.
#'   The default value is `seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 10)` for
#'   \code{LOOCV} and `seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 50)` for
#'   \code{LOOCV(h=b)}}
#'}
#' @return
#' @export
#' @importFrom nprobust lpbwselect
#' @import SuperLearner earth gam ranger rpart Rsolnp nnls
#' @examples
#' # Sample data
#' n <- 200; cols <- 3
#' W <- matrix(runif(n*cols), ncol=cols) # 200 * 3 matrix of covariates
#' A <- rnorm(n, mean=W%*%rnorm(cols))   # 200 * 1 vector of treatment variable
#' Y <- rnorm(n, mean = sin(A)+W%*%rnorm(cols))  # 200 * 1 vector of response variable
#' res <- debiased_ate_inference(Y, A, W)
debiased_ate_inference <- function(Y, A, W, tau=1,
                                   eval.pts.1=NULL, eval.pts.2=NULL, ...){
  # Parse control inputs ------------------------------------------------------
  control <- .parse.debiased_inference(...)
  kernel.type <- control$kernel.type
  if (is.null(eval.pts.1)){
    eval.pts.1 <- seq(quantile(A, 0.5), quantile(A, 0.95),  length.out=30)
    eval.pts.2 <- rep(quantile(A, 0.05),  length.out=30)
  }
  if (length(eval.pts.1) == 1){
    eval.pts.1 <- rep(eval.pts.1, length(eval.pts.2))
  }
  if (length(eval.pts.2) == 1){
    eval.pts.2 <- rep(eval.pts.2, length(eval.pts.1))
  }

  # Compute a nuisance estimators if not provided -----------------------------
  if(is.null(control$mu)){
    mu <- .fit.regression(Y, A, W)
  }
  else{
    mu <- control$mu
  }
  if(is.null(control$mu)){
    g <- .fit.density(A, W)
  }
  else{
    g <- control$g
  }

  # Estimate pseudo-outcome sequence ------------------------------------------
  ord <- order(A)
  Y <- Y[ord]; W <- W[ord,]; A <- A[ord]; n <- length(A)
  pseudo.res <- .pseudo.outcomes(Y, A, W, mu, g)
  pseudo.out <- pseudo.res$pseudo.outcome
  muhat.mat <- pseudo.res$muhat.mat
  mhat.obs <- pseudo.res$mhat.obs

  # Compute bandwidth ---------------------------------------------------------
  bw.eval <- unique(c(eval.pts.1, eval.pts.2))
  ord <- order(bw.eval); bw.eval <- bw.eval[ord]
  if(control$bandwidth.method == "LOOCV"){
    if (is.null(control$bw.seq)){
      bw.seq <- seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 10)
    }
    else{
      bw.seq <- control$bw.seq
    }
    bw.seq.h <- rep(bw.seq, length(bw.seq))
    bw.seq.b <- rep(bw.seq, each=length(bw.seq))
    risk <- mapply(function(h, b){
      .robust.loocv(A, pseudo.out, h, b, eval.pt=bw.eval,
                    kernel.type=kernel.type)}, bw.seq.h, bw.seq.b)
    h.opt <- bw.seq.h[which.min(risk)]; b.opt <- bw.seq.b[which.min(risk)]
  }
  else if(control$bandwidth.method == "LOOCV(h=b)"){
    if (is.null(control$bw.seq)){
      bw.seq <- seq(0.1*n^{-1/5}, 1*n^{-1/5}, length.out = 50)
    }
    else{
      bw.seq <- control$bw.seq
    }
    risk <- mapply(function(h, b){
      .robust.loocv(A, pseudo.out, h, b, eval.pt=bw.eval,
                    kernel.type=kernel.type)}, bw.seq, bw.seq)
    h.opt <- bw.seq[which.min(risk)]; b.opt <- bw.seq[which.min(risk)]
  }
  else{
    h.opt <- lpbwselect(pseudo.out, A, eval=bw.eval,
                        bwselect="imse-dpi")$bws[,2]
    b.opt <- h.opt/tau
  }

  # Fit debiased local linear regression --------------------------------------
  est.res.1 <- .lprobust(A, pseudo.out, h.opt, b.opt,
                         eval=eval.pts.1, kernel.type=kernel.type)
  est.res.2 <- .lprobust(A, pseudo.out, h.opt, b.opt,
                         eval=eval.pts.2, kernel.type=kernel.type)

  # Estimate influence function sequence --------------------------------------
  kern <- function(t){.kern(t, kernel=kernel.type)}
  inf.cov <- mapply(function(a1.val, a2.val, h.val, b.val){
    .compute.rinfl.func(pseudo.out, A, a1.val, h.val, b.val,
                       kern, muhat.mat, mhat.obs)-
      .compute.rinfl.func(pseudo.out, A, a2.val, h.val, b.val,
                         kern, muhat.mat, mhat.obs)
  }, eval.pts.1, eval.pts.2, h.opt, b.opt)
  se.cov <- apply(inf.cov, 2, sd)/sqrt(n)
  alpha <- control$alpha.pts
  ci.ll.p <- (est.res.1[,"theta.hat"]-est.res.2[,"theta.hat"])-
    qnorm(1-alpha/2)*se.cov
  ci.ul.p <- (est.res.1[,"theta.hat"]-est.res.2[,"theta.hat"])+
    qnorm(1-alpha/2)*se.cov

  # Output data.frame ---------------------------------------------------------
  res <- data.frame(
    eval.pts.1=eval.pts.1,
    eval.pts.2=eval.pts.2,
    theta.1=est.res.1[,"theta.hat"],
    theta.2=est.res.2[,"theta.hat"],
    theta.eff=est.res.1[,"theta.hat"]-est.res.2[,"theta.hat"],
    ci.ul.pts=ci.ul.p,
    ci.ll.pts=ci.ll.p,
    bias.1=est.res.1[,"b.hat"],
    bias.2=est.res.2[,"b.hat"],
    if.val=se.cov)
  res$h <- h.opt; res$b <- b.opt
  return(res)
}
