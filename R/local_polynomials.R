# File: local_polynomials.R
# Author: Kenta Takatsu
# Description:
#   function in this file corresponds to the debiased local linear regression
#   the implementation is similar to nprobust package who studied debiasing
#   in nonparametric inference
# References:
#   Takatsu K., and Westling T., (2023).
#   "Debiased inference for a covariate-adjusted regression function"
#   Calonico, Cattaneo and Farrell., (2019).
#   nprobust: Nonparametric Kernel-Based Estimation and Robust Bias-Corrected Inference.

#' Title: .lprobust.quad
#'
#' @param x \code{n x 1} numeric vector of observed covariate values.
#' @param y \code{n x 1} numeric vector of observed outcome values.
#' @param h the bandwidth parameter for the local linear estimator
#' @param b the bandwidth parameter for the bias estimation
#' @param eval.pt the sequence of points on which the regression values are
#' estimated
#' @param kernel.type The choice of a kernel function. It must be one
#'   of the following: \code{epa} (epanechnikov), \code{tri} (triangular),
#'   \code{uni} (uniform) and \code{gau} (Gaussian). The default is \code{epa}.
#'
#' @return
.lprobust.quad <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]; n<-length(y)
  if (is.null(eval.pt)){
    x.max <- max(x); x.min <- min(x)
    eval.pt <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval.pt)
  Estimate <- matrix(NA, neval, 4)
  colnames(Estimate) <- c("eval", "theta.hat", "mu.hat", "b.hat")
  for (i in 1:neval){
    # this part of code is simply avoiding a singular matrix
    bw.min   <- sort(abs(x-eval.pt[i]))[4]
    h0 <- max(h, bw.min); b0 <- max(b, bw.min)
    k.h   <- .kern((x-eval.pt[i])/h0, kernel.type)/h0
    k.b   <- .kern((x-eval.pt[i])/b0, kernel.type)/b0
    ind.h <- k.h>0;  ind.b <- k.b>0; ind <- ind.b
    if (h>b) ind <- ind.h
    eN  <- sum(ind); eY  <- y[ind]; eX  <- x[ind]
    K.h <- k.h[ind]; K.b <- k.b[ind]
    W.b <- matrix(NA, eN, 3) # (1, u, u^2)
    for (j in 1:3)  W.b[, j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h <- matrix(NA, eN, 2) # (1, u)
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    invD.b  <- .qrXXinv((sqrt(K.b)*W.b))
    invD.h  <- .qrXXinv((sqrt(K.h)*W.h))
    # estimate fixed-h c.2
    u <- (eX-eval.pt[i])/h0
    c2.est <- (invD.h %*% crossprod(W.h*K.h,u^2))[1]
    beta.ll  <- invD.h%*%crossprod(W.h*K.h, eY) #R.p^T W.h Y
    beta.bias <- (h0/b0)^2*invD.b%*%crossprod(W.b*K.b, eY)*c2.est
    tau <- beta.ll[1,1] # local linear term
    tau.bias <- beta.bias[3,1] # estimated bias
    Estimate[i,] <- c(eval.pt[i], tau-tau.bias, tau, tau.bias)
  }
  Estimate
}
