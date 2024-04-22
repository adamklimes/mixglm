test_that("setBUGSVariableName stops when it should", {
  expect_error(setBUGSVariableName(function(x) x))
})
test_that("setBUGSVariableName with a legitimate name", {
  x <- "varName"
  expect_equal(setBUGSVariableName(x), x)
})
test_that("setBUGSVariableName with illegitimate", {
  expect_equal(setBUGSVariableName("varName1"), "varNameOne")
  expect_equal(setBUGSVariableName("var$name"), "var_name")
})

