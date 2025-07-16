## 3. ------ DEFINE CUSTOM DISTRIBUTION FUNCTIONS ------

### 3.2. ==== Create NIMBLE functions for the mixture model (normal error) ====
## 3.2.1. dnormStateValueMembership density function ----
#' @title Probability Density of a response variable from Gaussian mixture (Normal Error)
#'
#' @description Calculate the probability density of a response variable from Gaussian mixture.
#'
#' @param x A scalar value containing the value of response variable
#' @param stateVal A numeric vector containing the mean values of the mixture components
#' @param statePrec A numeric vector containing the precision (reciprocal of the variance) of the
#' mixture components
#' @param stateProb A numeric vector containing the probabilities of belonging
#' to each mixture component (internally normalised to one)
#' @param log If \code{TRUE} then return the log density
#'
#' @return A scalar containing the probability density
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
dnormStateValueMembership <- nimbleFunction(
  run = function(
    x = double(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1),
    log = integer(0, default = 0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(stateProb)[1]))
    if(numStates < 2) {
      stop("invalid number of mixture components (there must be at least two)")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- stateVal#parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(-maxDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- statePrec#parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- stateProb#parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Calculate the total of the state probabilities
    fullProb <- sum(inStateProb[1:numStates])
    # Intialise a vector of conditional probabilities for each state
    condStateProb <- numeric(length = numStates)
    for(stateIter in 1:numStates) {
      condStateProb[stateIter] <- inStateProb[stateIter] * dnorm(x, inStateVal[stateIter], pow(inStatePrec[stateIter], -0.5), FALSE)
    }
    # Intialise the output probability
    outProb <- sum(condStateProb[1:numStates]) / fullProb
    if(log) {
      outProb <- log(outProb)
    }
    return(outProb)
  }
)
## 3.2.2. dnormStateValueMembership simulation function ----
#' @title Simulate a Variable from a Gaussian mixture (Normal Error)
#'
#' @description Simulate a variable from a Gaussian mixture.
#'
#' @param n An integer scalar containing the number of variables to simulate (NIMBLE currently
#' only allows n = 1)
#' @param stateVal A numeric vector containing the mean values of mixture components
#' @param statePrec A numeric vector containing the precision (reciprocal of the
#' variance) of the mixture components
#' @param stateProb A numeric vector containing the probabilities of belonging
#' to each mixture component (internally normalised to one)
#'
#' @return A scalar containing the simulated response variable
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
rnormStateValueMembership <- nimbleFunction(
  run = function(
    n = integer(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(stateProb)[1]))
    if(numStates < 2) {
      stop("invalid number of mixture components (there must be at least two)")
    }
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("this function only allows n = 1; using n = 1")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- stateVal#parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(-maxDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- statePrec#parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- stateProb#parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Draw a state to be in
    outState <- rcat(1, inStateProb)
    # Using this state: draw the state value
    return(rnorm(1, inStateVal[outState], pow(inStatePrec[outState], -0.5)))
  }
)

### 3.3. ==== Create NIMBLE functions for the mixture model (gamma error) ====
## 3.3.1. dgammaStateValueMembership density function ----
#' @title Probability Density of a response variable from gamma mixture (Gamma Error)
#'
#' @description Calculate the probability density of a response variable from gamma mixture.
#'
#' @param x A scalar value containing the value of response variable
#' @param stateVal A numeric vector containing the mean values of mixture components
#' @param statePrec A numeric vector containing the precision (reciprocal of the
#' variance) of the mixture components
#' @param stateProb A numeric vector containing the probabilities of belonging
#' to each mixture component (internally normalised to one)
#' @param log If \code{TRUE} then return the log density
#'
#' @return A scalar containing the probability density
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
dgammaStateValueMembership <- nimbleFunction(
  run = function(
    x = double(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1),
    log = integer(0, default = 0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(stateProb)[1]))
    if(numStates < 2) {
      stop("invalid number of mixture components (there must be at least two)")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- stateVal#parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- statePrec#parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- stateProb#parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Calculate the total of the state probabilities
    fullProb <- sum(inStateProb[1:numStates])
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inShape <- inStatePrec[1:numStates] * inStateVal[1:numStates] * inStateVal[1:numStates]
#    inScale <- 1.0 / (inStatePrec[1:numStates] * inStateVal[1:numStates])
    inRate <- inStatePrec[1:numStates] * inStateVal[1:numStates]
    # Intialise a vector of conditional probabilities for each state
    condStateProb <- numeric(length = numStates)
    for(stateIter in 1:numStates) {
      condStateProb[stateIter] <- inStateProb[stateIter] * dgamma(x, inShape[stateIter], inRate[stateIter], log = FALSE)
    }
    # Intialise the output probability
    outProb <- sum(condStateProb[1:numStates]) / fullProb
    if(log) {
      outProb <- log(outProb)
    }
    return(outProb)
  }
)
## 3.3.2. dgammaStateValueMembership simulation function ----
#' @title Simulate a Variable from a gamma mixture (Gamma Error)
#'
#' @description Simulate a variable from a gamma mixture.
#'
#' @param n An integer scalar containing the number of variables to simulate (NIMBLE currently
#' only allows n = 1)
#' @param stateVal A numeric vector containing the mean values of the mixture components
#' @param statePrec A numeric vector containing the precision (reciprocal of the
#' variance) of the mixture components
#' @param stateProb A numeric vector containing the probabilities of belonging
#' to each mixture component (internally normalised to one)
#'
#' @return A scalar containing the simulated response variable
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
rgammaStateValueMembership <- nimbleFunction(
  run = function(
    n = integer(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(stateProb)[1]))
    if(numStates < 2) {
      stop("invalid number of mixture components (there must be at least two)")
    }
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("this function only allows n = 1; using n = 1")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- stateVal#parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- statePrec#parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- stateProb#parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inShape <- inStatePrec[1:numStates] * inStateVal[1:numStates] * inStateVal[1:numStates]
    # inScale <- 1.0 / (inStatePrec[1:numStates] * inStateVal[1:numStates])
    inRate <- inStatePrec[1:numStates] * inStateVal[1:numStates]
    # Draw a state to be in
    outState <- rcat(1, inStateProb)
    # Using this state: draw the state value
    return(rgamma(1, inShape[outState], rate = inRate[outState]))
  }
)

### 3.4. ==== Create NIMBLE functions for the beta mixture (beta error) ====
## 3.4.1. dbetaStateValueMembership density function ----
#' @title Probability Density of a response variable from beta mixture (Beta Error)
#'
#' @description Calculate the probability density of a response variable from beta mixture.
#'
#' @param x A scalar value containing the value of the response variable
#' @param stateVal A numeric vector containing the mean values of the mixture components
#' @param statePrec A numeric vector containing the precision (reciprocal of the
#' variance) of the mixture components
#' @param stateProb A numeric vector containing the probabilities of belonging
#' to each mixture component (internally normalised to one)
#' @param log If \code{TRUE} then return the log density
#'
#' @return A scalar containing the probability density
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
dbetaStateValueMembership <- nimbleFunction(
  run = function(
    x = double(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1),
    log = integer(0, default = 0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(stateProb)[1]))
    if(numStates < 2) {
      stop("invalid number of mixture components (there must be at least two)")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- stateVal#parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), 1.0 - minNonZeroDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- statePrec#parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- stateProb#parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Calculate the total of the state probabilities
    fullProb <- sum(inStateProb[1:numStates])
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inAlpha <- inStateVal[1:numStates] * inStateVal[1:numStates] * (1.0 - inStateVal[1:numStates]) * inStatePrec[1:numStates] - inStateVal[1:numStates]
    inBeta <- inStateVal[1:numStates] * (1.0 - inStateVal[1:numStates]) * (1.0 - inStateVal[1:numStates]) * inStatePrec[1:numStates] + inStateVal[1:numStates] - 1.0
    # Intialise a vector of conditional probabilities for each state
    condStateProb <- numeric(length = numStates)
    for(stateIter in 1:numStates) {
      condStateProb[stateIter] <- inStateProb[stateIter] * dbeta(x, inAlpha[stateIter], inBeta[stateIter], FALSE)
    }
    # Intialise the output probability
    outProb <- sum(condStateProb[1:numStates]) / fullProb
    if(log) {
      outProb <- log(outProb)
    }
    return(outProb)
  }
)
## 3.4.2. dbetaStateValueMembership simulation function ----
#' @title Simulate a Variable from a beta mixture (Beta Error)
#'
#' @description Simulate a variable from a beta mixture.
#'
#' @param n An integer scalar containing the number of variables to simulate (NIMBLE currently
#' only allows n = 1)
#' @param stateVal A numeric vector containing the mean values of the mixture components
#' @param statePrec A numeric vector containing the precision (reciprocal of the
#' variance) of the mixture components
#' @param stateProb A numeric vector containing the probabilities of belonging
#' to each mixture component (internally normalised to one)
#'
#' @return A scalar containing the simulated response variable
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
rbetaStateValueMembership <- nimbleFunction(
  run = function(
    n = integer(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(stateProb)[1]))
    if(numStates < 2) {
      stop("invalid number of mixture components (there must be at least two)")
    }
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("this function only allows n = 1; using n = 1")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- stateVal#parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), 1.0 - minNonZeroDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- statePrec#parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- stateProb#parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inAlpha <- inStateVal[1:numStates] * inStateVal[1:numStates] * (1.0 - inStateVal[1:numStates]) * inStatePrec[1:numStates] - inStateVal[1:numStates]
    inBeta <- inStateVal[1:numStates] * (1.0 - inStateVal[1:numStates]) * (1.0 - inStateVal[1:numStates]) * inStatePrec[1:numStates] + inStateVal[1:numStates] - 1.0
    # Draw a state to be in
    outState <- rcat(1, inStateProb)
    # Using this state: draw the state value
    return(rbeta(1, inAlpha[outState], inBeta[outState]))
  }
)

### 3.5. ==== Create NIMBLE functions for the negative binomial mixture (negative binomial error) ====
## 3.5.1. dnegbinStateValueMembership density function ----
#' @title Probability Density of a response variable from negative binomial mixture (Negative Binomial Error)
#'
#' @description Calculate the probability density of a response variable from negative binomial mixture.
#'
#' @param x A scalar value containing the value of the response variable
#' @param stateVal A numeric vector containing the mean values of the mixture components
#' @param statePrec A numeric vector containing the precision (reciprocal of the
#' variance) of the mixture components
#' @param stateProb A numeric vector containing the probabilities of belonging
#' to each mixture component (internally normalised to one)
#' @param log If \code{TRUE} then return the log density
#'
#' @return A scalar containing the probability density
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
dnegbinStateValueMembership <- nimbleFunction(
  run = function(
    x = double(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1),
    log = integer(0, default = 0)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(stateProb)[1]))
    if(numStates < 2) {
      stop("invalid number of mixture components (there must be at least two)")
    }
    # Recycle the input vector if necessary and check to ensure the values are valid
    inStateVal <- stateVal#parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- statePrec#parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- stateProb#parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Calculate the total of the state probabilities
    fullProb <- sum(inStateProb[1:numStates])
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inProb <- inStateVal[1:numStates] * inStatePrec[1:numStates]
    inSize <- inStateVal[1:numStates] * inStateVal[1:numStates] * inStatePrec[1:numStates] / (1.0 - inProb[1:numStates])
    # Intialise a vector of conditional probabilities for each state
    condStateProb <- numeric(length = numStates)
    for(stateIter in 1:numStates) {
      # Test to ensure that the parameterisation is valid
      if(inProb[stateIter] < 0.0 | inProb[stateIter] > 1.0 | inSize[stateIter] < 0.0) {
        if(log) {
          return(-Inf)
        } else {
          return(0.0)
        }
      }
      # Calculate the conditional probability
      condStateProb[stateIter] <- inStateProb[stateIter] * dnbinom(x, inSize[stateIter], inProb[stateIter], log = FALSE)
    }
    # Intialise the output probability
    outProb <- sum(condStateProb[1:numStates]) / fullProb
    if(log) {
      outProb <- log(outProb)
    }
    return(outProb)
  }
)
## 3.5.2. dnegbinStateValueMembership simulation function ----
#' @title Simulate a Variable from a negative binomial mixture (Negative Binomial Error)
#'
#' @description Simulate a variable from negative binomial mixture.
#'
#' @param n An integer scalar containing the number of variables to simulate (NIMBLE currently
#' only allows n = 1)
#' @param stateVal A numeric vector containing the mean values of the mixture components
#' @param statePrec A numeric vector containing the precision (reciprocal of the
#' variance) of the mixture components
#' @param stateProb A numeric vector containing the probabilities of belonging
#' to each mixture component (internally normalised to one)
#'
#' @return A scalar containing the simulated response variable
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @export
#'
rnegbinStateValueMembership <- nimbleFunction(
  run = function(
    n = integer(0),
    stateVal = double(1),
    statePrec = double(1),
    stateProb = double(1)
  ) {
    # Specify the return type dimensionality
    returnType(double(0))
    # Retrieve the number of states in the model
    numStates <- max(dim(stateVal)[1], max(dim(statePrec)[1], dim(stateProb)[1]))
    if(numStates < 2) {
      stop("invalid number of mixture components (there must be at least two)")
    }
    # Ensure that only one sample is requested
    if(n <= 0) {
      stop("the number of requested samples must be above zero")
    } else if(n > 1) {
      print("this function only allows n = 1; using n = 1")
    }
    # Recycle the input vector if neccessary and check to ensure the values are valid
    inStateVal <- stateVal#parameterTest_vectorDouble(stateVal, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStatePrec <- statePrec#parameterTest_vectorDouble(statePrec, numeric(length = 2, value = c(minNonZeroDouble(), maxDouble())), integer(length = 1, value = numStates), 1)
    inStateProb <- stateProb#parameterTest_vectorDouble(stateProb, numeric(length = 2, value = c(0.0, maxDouble())), integer(length = 1, value = numStates), 1)
    # Reparameterise in terms of the canoconical NIMBLE parameterisation
    inProb <- inStateVal[1:numStates] * inStatePrec[1:numStates]
    inSize <- inStateVal[1:numStates] * inStateVal[1:numStates] * inStatePrec[1:numStates] / (1.0 - inProb[1:numStates])
    # Draw a state to be in
    outState <- rcat(1, inStateProb)
    # Test to ensure that the parameterisation is valid
    if(inProb[outState] < 0.0 | inProb[outState] > 1.0 | inSize[outState] < 0.0) {
      stop("invalid parameter values given for the negative binomial distribution")
    }
    # Using this state: draw the state value
    return(rnbinom(1, inSize[outState], inProb[outState]))
  }
)
