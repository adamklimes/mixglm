test_that("findMin stops when it should", {
  expect_error(findMin(LETTERS))
})
test_that("findMin with inversed order of x", {
  x <- rnorm(10)
  x2 <- rev(x)
  expect_equal(sort(x[findMin(x)]), sort(x2[findMin(x2)]))
})
test_that("findMin with extremes = FALSE", {
  x <- sort(rnorm(10))
  expect_equal(findMin(x, extremes = FALSE), numeric(0))
})
test_that("findMin with monotonic sequences", {
  x <- sort(rnorm(10))
  expect_equal(findMin(x), 1)
  x2 <- sort(rnorm(10), decreasing = TRUE)
  expect_equal(findMin(x2), length(x2))
})
