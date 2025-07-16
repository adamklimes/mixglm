# Registers new distributions when the package is loaded
.onAttach <- function(libname, pkgname) {
distributionList <- list(
  ## 1.1.2. Define the dnormStateValueMembership distribution ----
  dnormStateValueMembership = list(
    # Define the BUGS code to call the distribution
    BUGSdist = "dnormStateValueMembership(stateVal, statePrec, stateProb)",
    # Set the input and output types and dimension structure
    types = c(
      "value = double(0)", "stateVal = double(1)", "statePrec = double(1)",
      "stateProb = double(1)"),
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )
  ## 1.1.3. Define the dgammaStateValueMembership distribution ----
  , dgammaStateValueMembership = list(
  # Define the BUGS code to call the distribution
  BUGSdist = "dgammaStateValueMembership(stateVal, statePrec, stateProb)",
    # Set the input and output types and dimension structure
    types = c(
      "value = double(0)", "stateVal = double(1)", "statePrec = double(1)",
      "stateProb = double(1)"),
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    # Specify the range of values
    range = c(0, Inf),
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )
  ## 1.1.4. Define the dbetaStateValueMembership distribution ----
  , dbetaStateValueMembership = list(
    # Define the BUGS code to call the distribution
    BUGSdist = "dbetaStateValueMembership(stateVal, statePrec, stateProb)",
    # Set the input and output types and dimension structure
    types = c(
      "value = double(0)", "stateVal = double(1)", "statePrec = double(1)",
      "stateProb = double(1)"),
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    # Specify the range of values
    range = c(0, 1),
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )
  ## 1.1.5. Define the dnegbinStateValueMembership distribution ----
  , dnegbinStateValueMembership = list(
    # Define the BUGS code to call the distribution
    BUGSdist = "dnegbinStateValueMembership(stateVal, statePrec, stateProb)",
    # Set the input and output types and dimension structure
    types = c(
      "value = double(0)", "stateVal = double(1)", "statePrec = double(1)",
      "stateProb = double(1)"),
    # Define the cumulative probability and quantile function availability
    pqAvail = FALSE,
    # It is a discrete-valued distribution
    discrete = TRUE,
    # Specify the range of values
    range = c(0, Inf),
    mixedSizes = TRUE   # Turn off warnings about possible dimension mismatch
  )
)
suppressMessages({
  if(exists("distributions", nimbleUserNamespace)) {
    # Find those distributions that are defined in this source file and see if they are already registered
    distributionNames <- names(distributionList)
    isDefinedDist <- distributionNames %in% nimbleUserNamespace$distributions$namesVector
    if(any(isDefinedDist)) {
      # If the distributions are registered then deregister them
      deregisterDistributions(distributionNames[isDefinedDist])
    }
  }
  registerDistributions(distributionList)
})

}
