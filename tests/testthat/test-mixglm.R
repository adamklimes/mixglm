test_that("mixglm stops when it should", {
  n <- 200
  x <- rnorm(n)
  y <- rnorm(n)
  dat <- data.frame(x, y)
  expect_error(mixglm( ~ x, ~ 1, ~ 1, dat, 2))
})
