test_that("mixglmSimulation stops when it should", {
  dat <- data.frame(x = rnorm(200), y = rnorm(200))
  cfValues <- list(
    "intercept_stateVal" = c(-1.0, 0.5, 0.5),
    "x_stateVal" = c(-0.05, 0.05, 0.1),
    "intercept_stateProb" = c(NA, -1.0, 1.0),
    "intercept_statePrec" = c(0.5, 1.0, 1.5)
  )
  expect_error(mixglmSimulation(1, y ~ x, ~ z, ~ 1, dat, 2, gaussian, cfValues))
})
test_that("mixglmSimulation for different predictors", {
  numDataPoints <- 500
  padZero <- function(inVec, size) {
    outVec <- rep(0.0, size)
    outVec[1:min(size, length(inVec))] <- inVec[1:min(size, length(inVec))]
    outVec
  }
  testFrame <- data.frame(
    respVariable = rep(NA, numDataPoints),
    covA = seq(0.0, 1.0, length.out = numDataPoints),
    covB = padZero(c(seq(0.0, 1.0, length.out = floor(numDataPoints / 2)), seq(1.0, 0.0, length.out = floor(numDataPoints / 2))), numDataPoints),
    covC = as.factor(c(rep("levelA", floor(numDataPoints / 3)), rep("levelB", floor(numDataPoints / 3)), rep("levelC", numDataPoints - 2 * floor(numDataPoints / 3))))
  )
  valFormula <- respVariable ~ covB
  precFormula <- ~ covC
  probFormula <- ~ covA
  coefficientValues <- list(
    "intercept_stateVal" = c(-1.0, 0.5, 0.5),
    "covB_stateVal" = c(-0.05, 0.05, 0.1),
    "intercept_stateProb" = c(NA, -1.0, 1.0),
    "covA_stateProb" = c(NA, 1.0, -1.0),
    "intercept_statePrec" = c(0.5, 1.0, 1.5),
    "covClevelB_statePrec" = c(0.1, -0.1, 0.5),
    "covClevelC_statePrec" = c(0.4, -0.4, 0.5)
  )
  simOutput <- mixglmSimulation(1, valFormula, probFormula, precFormula, testFrame, 3, gaussian, coefficientValues)
  expect_equal(length(simOutput$respVariable[, 1]), numDataPoints)
})
