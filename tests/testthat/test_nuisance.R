test_that("testing mu 3D", {
  X <- rnorm(100, 0, 1)
  W <- data.frame(matrix(rnorm(300, 0, 1), ncol=3))
  Y <- rnorm(100, X+W[,1], 1)
  expect_error(.fit.regression(Y, X, W), NA)
})

test_that("testing mu 1D", {
  X <- rnorm(100, 0, 1)
  W <- data.frame(matrix(rnorm(100, 0, 1), ncol=1))
  Y <- rnorm(100, X+W[,], 1)
  expect_error(.fit.regression(Y, X, W), NA)
})


