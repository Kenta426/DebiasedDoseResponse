# File: influence_function.R
# Author: Kenta Takatsu
# Description:
#   Two functions in this file corresponds to the influence functions
#   for debiased local linear estimators and local linear estimators.
# References:
#   Takatsu K., and Westling T., (2023).
#   "Debiased inference for a covariate-adjusted regression function"
#   Kennedy EH., Ma Z., McHugh MD, Small DS (2017).
#   "Nonparametric methods for doubly robust estimation of
#   continuous treatment effects"

# Function 1: .compute.rinfl.func
# Purpose:
#   Compute the estimate of influence function sequence for debiased local
#   linear estimator of pseudo-outcome sequences.
# Arguments:
#   Y:
#   A:
#   a:
#   h:
#   b:
#   kern:
#   muhat.mat:
#   mhat.obs:
# Output: [Brief description of output]
.compute.rinfl.func <- function(Y, A, a, h, b, kern, muhat.mat, mhat.obs){

  n <- length(A); bw.min <- sort(abs(A - a))[21]
  h <- max(h, bw.min); b <- max(b, bw.min)
  a.std.h <- (A - a)/h; kern.std.h <- kern(a.std.h)/h
  a.std.b <- (A - a)/b; kern.std.b <- kern(a.std.b)/b

  # Compute 2x2 inverse matrix ------------------------------------------------
  c0.h <- mean(kern.std.h)
  c1.h <- mean(kern.std.h * a.std.h)
  c2.h <- mean(kern.std.h * a.std.h^2)
  Dh <- matrix(c(c0.h, c1.h,
                 c1.h, c2.h), nrow = 2)
  Dh.inv <- solve(Dh)

  # Compute 4x4 inverse matrix ------------------------------------------------
  c0.b <- mean(kern.std.b)
  c1.b <- mean(kern.std.b * a.std.b)
  c2.b <- mean(kern.std.b * a.std.b^2)
  c3.b <- mean(kern.std.b * a.std.b^3)
  c4.b <- mean(kern.std.b * a.std.b^4)
  c5.b <- mean(kern.std.b * a.std.b^5)
  c6.b <- mean(kern.std.b * a.std.b^6)
  Db <- matrix(c(c0.b, c1.b, c2.b, c3.b,
                 c1.b, c2.b, c3.b, c4.b,
                 c2.b, c3.b, c4.b, c5.b,
                 c3.b, c4.b, c5.b, c6.b), nrow = 4)
  Db.inv <- solve(Db)

  # Estimate local linear components ------------------------------------------
  g2.h <- (A - a)/h
  int1.h <- colMeans(kern.std.h * (muhat.mat - mhat.obs))
  int2.h <- colMeans(g2.h * kern.std.h * (muhat.mat - mhat.obs))
  gamma.h <- coef(lm(Y ~ a.std.h, weights = kern.std.h))
  res.h <- Y - (gamma.h[1] + gamma.h[2] * a.std.h)
  inf.fn <- t(Dh.inv %*% rbind(res.h * kern.std.h + int1.h,
                               g2.h * res.h * kern.std.h + int2.h))

  # Estimate local polynomial components --------------------------------------
  g2.b <- (A - a)/b
  g3.b <- ((A - a)/b)^2
  g4.b <- ((A - a)/b)^3
  int1.b <- colMeans(kern.std.b * (muhat.mat - mhat.obs))
  int2.b <- colMeans(g2.b * kern.std.b * (muhat.mat - mhat.obs))
  int3.b <- colMeans(g3.b * kern.std.b * (muhat.mat - mhat.obs))
  int4.b <- colMeans(g4.b * kern.std.b * (muhat.mat - mhat.obs))
  model.b <- lm(Y ~ poly(a.std.b, 3), weights = kern.std.b)
  res.b <- Y - predict(model.b)
  inf.fn.robust <- t(Db.inv %*% rbind(res.b * kern.std.b + int1.b,
                                         g2.b * res.b * kern.std.b + int2.b,
                                         g3.b * res.b * kern.std.b + int3.b,
                                         g4.b * res.b * kern.std.b + int4.b))

  # Build the influence function ----------------------------------------------
  c2 <- integrate(function(u){u^2 * kern(u)}, -Inf, Inf)$value
  return (inf.fn[,1] - (h/b)^2 * c2 * inf.fn.robust[,3])
}

# Compute influence function values for a local linear regression
.compute.infl.func <- function(Y, A, a, h, kern, muhat.mat, mhat.obs){
  n <- length(A)
  bw.min <- sort(abs(A - a))[21]
  h <- max(h, bw.min)
  a.std <- (A - a)/h; kern.std <- kern(a.std)/h

  # Compute 2x2 inverse matrix ------------------------------------------------
  c0 <- mean(kern.std)
  c1 <- mean(kern.std * a.std)
  c2 <- mean(kern.std * a.std^2)
  Dh <- matrix(c(c0, c1,
                 c1, c2), nrow = 2)
  Dh.inv <- solve(Dh)

  # Estimate integrals and other components -----------------------------------
  g2 <- (A - a)/h
  int1 <- colMeans(kern.std * (muhat.mat - mhat.obs))
  int2 <- colMeans(g2 * kern.std * (muhat.mat - mhat.obs))
  gamma.h <- coef(lm(Y ~ a.std, weights = kern.std))
  res.h <- Y - gamma.h[1] - gamma.h[2] * a.std

  # Build the influence function ----------------------------------------------
  inf.fn <- t(Dh.inv %*% rbind(res.h * kern.std + int1,
                                  g2 * res.h * kern.std + int2))
  return (inf.fn[,1])
}
