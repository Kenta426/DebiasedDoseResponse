# File: influence_function.R
# Author: Kenta Takatsu
# Description:
#   function in this file corresponds to the influence functions
#   for debiased local linear estimators and local linear estimators.
# References:
#   Takatsu K., and Westling T., (2023).
#   "Debiased inference for a covariate-adjusted regression function"

#' Title: .compute.rinfl.func.quad
#'
#' @param Y \code{n x 1} numeric vector of estimated pseudo-outcomes.
#' @param A \code{n x 1} numeric vector of observed exposure values.
#' @param a a scalar corresponding to a point of exposure to run inference.
#' @param h the bandwidth parameter for the local linear estimator.
#' @param b the bandwidth parameter for the bias estimation.
#' @param kern function R to R corresponding to the kernel
#' @param mu.hat estimated mu(A, W) as matrix. output from .pseudo.outcomes in
#' helper.R.
#' @param theta.hat estimated int mu(A, W) dW. output from .pseudo.outcomes in
#' helper.R.
#'
#' @return
#'
.compute.rinfl.func.quad <- function(Y, A, a, h, b, kern, mu.hat, theta.hat){
  n <- length(A)
  a.std.h <- (A - a)/h; kern.std.h <- kern(a.std.h)/h
  a.std.b <- (A - a)/b; kern.std.b <- kern(a.std.b)/b
  # Compute 2x2 inverse matrix ------------------------------------------------
  c0.h <- mean(kern.std.h)
  c1.h <- mean(kern.std.h * a.std.h)
  c2.h <- mean(kern.std.h * a.std.h^2)
  Dh <- matrix(c(c0.h, c1.h,
                 c1.h, c2.h), nrow = 2)
  # Compute 3x3 inverse matrix ------------------------------------------------
  c0.b <- mean(kern.std.b)
  c1.b <- mean(kern.std.b * a.std.b)
  c2.b <- mean(kern.std.b * a.std.b^2)
  c3.b <- mean(kern.std.b * a.std.b^3)
  c4.b <- mean(kern.std.b * a.std.b^4)
  Db <- matrix(c(c0.b, c1.b, c2.b,
                 c1.b, c2.b, c3.b,
                 c2.b, c3.b, c4.b), nrow = 3)
  # Estimate fixed-h c2 -------------------------------------------------------
  w.vec <- cbind(kern.std.h, kern.std.h* a.std.h)
  c2.vec <- (solve(Dh) %*% crossprod(w.vec, a.std.h^2))/n
  c2.est <- c2.vec[1]
  # the EIF for the fixed-h c2 ------------------------------------------------
  w.1 <- cbind(1, a.std.h)
  w.1.tilde <- cbind(kern.std.h*a.std.h^2, kern.std.h*a.std.h^3)
  term1 <- (w.1 %*% solve(Dh))[,1] * c(w.1 %*% c2.vec) * kern.std.h
  term2 <- (w.1.tilde %*% solve(Dh))[,1]
  deriv2 <- (solve(Db) %*%
               crossprod(cbind(1, a.std.b, a.std.b^2)*kern.std.b, Y)/n)[3,]
  if.c2 <- h^2/2 * deriv2 *(term2-term1)
  # Integral terms in EIF -----------------------------------------------------
  int1.h <- colMeans(kern.std.h * (mu.hat - theta.hat))
  int2.h <- colMeans(a.std.h * kern.std.h * (mu.hat - theta.hat))
  int1.b <- colMeans(kern.std.b * (mu.hat - theta.hat))
  int2.b <- colMeans(a.std.b * kern.std.b * (mu.hat - theta.hat))
  int3.b <- colMeans(a.std.b^2 * kern.std.b * (mu.hat - theta.hat))
  # Terms for lower case gamma. This corresponds to empirical residuals -------
  model.h <- lm(Y ~ a.std.h, weights = kern.std.h)
  model.b <- lm(Y ~ poly(a.std.b, 2), weights = kern.std.b)
  res.h <- residuals(model.h)
  res.b <- residuals(model.b)
  # EIF for linear estimator --------------------------------------------------
  inf.fn <- t(solve(Dh) %*% rbind(res.h * kern.std.h + int1.h,
                                  a.std.h * res.h * kern.std.h + int2.h))
  # EIF for bias estimator ----------------------------------------------------
  inf.fn.robust <- t(solve(Db) %*%
                       rbind(res.b * kern.std.b + int1.b,
                             a.std.b * res.b * kern.std.b + int2.b,
                             a.std.b^2 * res.b * kern.std.b + int3.b))
  # Putting together ----------------------------------------------------------
  return (inf.fn[,1]-(h/b)^2*c2.est*inf.fn.robust[,3]-if.c2)
}
