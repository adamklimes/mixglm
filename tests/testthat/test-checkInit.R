test_that("checkInit stops when it should", {
  n <- 200
  x <- rnorm(n)
  y <- rnorm(n)
  dat <- data.frame(x, y)
  initVals <- c(0, 0)
  expect_error(checkInit(y ~ x, y ~ x, dat, initVals, "gaussian", "identity", 2))
})
test_that("checkInit with legitimate initial values", {
  n <- 200
  x <- rnorm(n)
  y <- rnorm(n)
  dat <- data.frame(x, y)
  initVals <- list(
    intercept_stateVal = c(0,0),
    intercept_statePrec = c(0,0),
    intercept_stateProb = c(0,0),
    x_stateVal = c(0,0),
    x_statePrec = c(0,0),
    x_stateProb = c(0,0)
  )
  expect_equal(checkInit(y ~ x, y ~ x, dat, initVals, "gaussian", "identity", 2), initVals)
})
