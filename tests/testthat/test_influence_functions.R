test_that("test pseudo-outcomes", {
  # sample problem

  # data-generaing process
  lambda <- function(w, beta, kappa=0.1){
    c(kappa + 2 * (1 - kappa) * expit(w %*% beta))
  }
  expit <- function(x) 1 / (1 + exp(-x))
  G0inv <- function(u, lambda) {
    ret <- rep(NA, length(u))
    ones <- lambda == 1
    ret[ones] <- u[ones]
    ret[!ones] <- (-lambda[!ones] + sqrt(lambda[!ones]^2 + 4 * u[!ones] * (1 - lambda[!ones])))/
      (2 * (1 -lambda[!ones]))
    ret
  }
  g0 <- function(u, lambda) c(lambda + 2 * (1 - lambda) * u)

  # generate data
  W <- matrix(rnorm(200, 0, 1), ncol=2)
  A <- G0inv(runif(100), lambda(W, c(1, 0)))
  Y <- rnorm(100, A+W[,1],1)

  # oracle nuisances
  mu.hat <- function(a, w){
    a+w[,1]
  }
  g.hat <- function(a, w) g0(a, lambda(as.matrix(w), c(1,0)))
  expect_error(.pseudo.outcomes(Y, A, W, mu.hat, g.hat), NA)
})

test_that("test influence function (local linear)", {
  # sample problem

  # data-generaing process
  lambda <- function(w, beta, kappa=0.1){
    c(kappa + 2 * (1 - kappa) * expit(w %*% beta))
  }
  expit <- function(x) 1 / (1 + exp(-x))
  G0inv <- function(u, lambda) {
    ret <- rep(NA, length(u))
    ones <- lambda == 1
    ret[ones] <- u[ones]
    ret[!ones] <- (-lambda[!ones] + sqrt(lambda[!ones]^2 + 4 * u[!ones] * (1 - lambda[!ones])))/
      (2 * (1 -lambda[!ones]))
    ret
  }
  g0 <- function(u, lambda) c(lambda + 2 * (1 - lambda) * u)

  # generate data
  W <- matrix(rnorm(200, 0, 1), ncol=2)
  A <- G0inv(runif(100), lambda(W, c(1, 0)))
  Y <- rnorm(100, A+W[,1],1)

  # oracle nuisances
  mu.hat <- function(a, w){
    a+w[,1]
  }
  g.hat <- function(a, w) g0(a, lambda(as.matrix(w), c(1,0)))
  pseudo.res <- .pseudo.outcomes(Y, A, W, mu.hat, g.hat)
  pseudo.out <- pseudo.res$pseudo.outcome
  muhat.mat <- pseudo.res$muhat.mat
  mhat.obs <- pseudo.res$mhat.obs
  expect_error(.compute.infl.func(pseudo.out, A, 0.5, 0.1,
                                  .kern, muhat.mat, mhat.obs), NA)
})


test_that("test influence function (debiased local linear)", {
  # sample problem

  # data-generaing process
  lambda <- function(w, beta, kappa=0.1){
    c(kappa + 2 * (1 - kappa) * expit(w %*% beta))
  }
  expit <- function(x) 1 / (1 + exp(-x))
  G0inv <- function(u, lambda) {
    ret <- rep(NA, length(u))
    ones <- lambda == 1
    ret[ones] <- u[ones]
    ret[!ones] <- (-lambda[!ones] + sqrt(lambda[!ones]^2 + 4 * u[!ones] * (1 - lambda[!ones])))/
      (2 * (1 -lambda[!ones]))
    ret
  }
  g0 <- function(u, lambda) c(lambda + 2 * (1 - lambda) * u)

  # generate data
  W <- matrix(rnorm(200, 0, 1), ncol=2)
  A <- G0inv(runif(100), lambda(W, c(1, 0)))
  Y <- rnorm(100, A+W[,1],1)

  # oracle nuisances
  mu.hat <- function(a, w){
    a+w[,1]
  }
  g.hat <- function(a, w) g0(a, lambda(as.matrix(w), c(1,0)))
  pseudo.res <- .pseudo.outcomes(Y, A, W, mu.hat, g.hat)
  pseudo.out <- pseudo.res$pseudo.outcome
  muhat.mat <- pseudo.res$muhat.mat
  mhat.obs <- pseudo.res$mhat.obs
  expect_error(.compute.rinfl.func(pseudo.out, A, 0.5, 0.1, 0.1,
                                  .kern, muhat.mat, mhat.obs), NA)
})

