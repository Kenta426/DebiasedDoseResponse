test_that("compute hatmatrix (epa)", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.hatmatrix(x, y, 0.1, 0.1, eval.pt=NULL, kernel.type="epa"), NA)
})

test_that("compute hatmatrix (tri)", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.hatmatrix(x, y, 0.1, 0.1, eval.pt=NULL, kernel.type="tri"), NA)
})

test_that("compute hatmatrix (uni)", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.hatmatrix(x, y, 0.1, 0.1, eval.pt=NULL, kernel.type="uni"), NA)
})

test_that("compute hatmatrix (gau)", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.hatmatrix(x, y, 0.1, 0.1, eval.pt=NULL, kernel.type="gau"), NA)
})

test_that("loocv with hatmatrix trick (epa)", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.robust.loocv(x, y, 0.1, 0.1, eval.pt=NULL, kernel.type="epa"),
               NA)
})

test_that("loocv with hatmatrix trick (tri)", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.robust.loocv(x, y, 0.1, 0.1, eval.pt=NULL, kernel.type="tri"),
               NA)
})

test_that("loocv with hatmatrix trick (uni)", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.robust.loocv(x, y, 0.1, 0.1, eval.pt=NULL, kernel.type="uni"),
               NA)
})

test_that("loocv with hatmatrix trick (gau)", {
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  expect_error(.robust.loocv(x, y, 0.1, 0.1, eval.pt=NULL, kernel.type="gau"),
               NA)
})
