test_that("testing plot", {
  # sample problem
  n <- 200; cols <- 3
  W <- matrix(runif(n*cols), ncol = cols) # a 200 * 3 matrix of covariates
  A <- rnorm(n, mean=W%*%rnorm(cols)) # a 200 * 1 vector of treatment variable
  Y <- rnorm(n, mean = sin(A)) # a 200 * 1 vector of response variable
  res <- debiased_inference(Y, A, W)
  expect_error(plot_curve(res), NA)
})
