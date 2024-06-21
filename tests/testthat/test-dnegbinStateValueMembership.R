test_that("dnegbinStateValueMembership gives expected result for specified values", {
  x <- 1
  val <- c(1, 1.5)
  prec <- c(0.5, 0.5)
  prob <- c(0.5, 0.5)
  dens <- prob * dnbinom(x, val * val * prec / (1 - val * prec), val * prec, log = FALSE)
  expect_equal(dnegbinStateValueMembership(x, val, prec, prob), sum(dens))
  expect_equal(dnegbinStateValueMembership(x, val, prec, prob, TRUE), log(sum(dens)))
})
