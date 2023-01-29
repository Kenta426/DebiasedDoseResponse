test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("testing kernel", {
  expect_equal(.kern(0), 0.75)
})

test_that("testing local linear", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.loclinear(x, y, 0.1), NA)
})

test_that("testing local polinomial degree 1", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.locpoly(x, y, 0.1), NA)
})

test_that("testing local polinomial degree 2", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.locpoly(x, y, 0.1, degree=2), NA)
})

test_that("testing local polinomial degree 3", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.locpoly(x, y, 0.1, degree=3), NA)
})

test_that("testing debiased local linear regression", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.lprobust(x, y, 0.1, 0.1), NA)
})

test_that("testing hat matrix computation", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.hatmatrix(x, y, 0.1, 0.1), NA)
})

test_that("testing leave-one-out cross-validation", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.robust.loocv(x, y, 0.1, 0.1), NA)
})

