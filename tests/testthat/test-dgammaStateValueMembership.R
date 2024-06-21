test_that("dgammaStateValueMembership gives expected result for specified values", {
  x <- 3
  val <- c(1, 3)
  prec <- c(2, 2)
  prob <- c(0.5, 0.5)
  dens <- prob * dgamma(x, val * val * prec, val * prec, log = FALSE)
  expect_equal(dgammaStateValueMembership(x, val, prec, prob), sum(dens))
  expect_equal(dgammaStateValueMembership(x, val, prec, prob, TRUE), log(sum(dens)))
})
