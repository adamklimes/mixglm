test_that("dbetabin gives expected result for specified values", {
  n <- 10
  A <- 5
  B <- 4
  x <- 0.5
  dens <- lgamma(n + 1.0) + lgamma(x + A) + lgamma(n - x + B) +
    lgamma(A + B) - lgamma(x + 1.0) - lgamma(n - x + 1.0) -
    lgamma(n + A + B) - lgamma(A) - lgamma(B)
  expect_equal(dbetabin(x, A, B, n), exp(dens))
  expect_equal(dbetabin(x, A, B, n, log = TRUE), dens)
})
