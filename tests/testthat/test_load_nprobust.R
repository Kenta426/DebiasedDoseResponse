# make sure nprobust is loaded correctly

test_that("testing np robust", {
  require(nprobust)
  x <- runif(100, -1, 1)
  y <- rnorm(100, sin(pi*x), 0.1)
  x.eval <- seq(-0.8, 0.8, by=0.1)
  expect_error(
    nprobust::lpbwselect(x, y, eval=x.eval, bwselect="imse-dpi")$bws[,2], NA)
})
