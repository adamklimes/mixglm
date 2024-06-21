test_that("dbetaStateValueMembership gives expected result for specified values", {
  x <- 0.5
  val <- c(0.2, 0.6)
  prec <- c(8, 8)
  prob <- c(0.5, 0.5)
  dens <- prob * dbeta(x, val * val * (1.0 - val) * prec - val,
    val * (1.0 - val)^2 * prec + val - 1.0, FALSE)
  expect_equal(dbetaStateValueMembership(x, val, prec, prob), sum(dens))
  expect_equal(dbetaStateValueMembership(x, val, prec, prob, TRUE), log(sum(dens)))
})
