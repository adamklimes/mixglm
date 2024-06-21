test_that("dnormStateValueMembership gives expected result for specified values", {
  x <- 1.5
  val <- c(1, 3)
  prec <- c(2, 2)
  prob <- c(0.5, 0.5)
  dens <- prob * dnorm(x, val, 1/sqrt(prec), log = FALSE)
  expect_equal(dnormStateValueMembership(x, val, prec, prob), sum(dens))
  expect_equal(dnormStateValueMembership(x, val, prec, prob, TRUE), log(sum(dens)))
})
