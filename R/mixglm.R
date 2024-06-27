## 1. ------ DEFINE INTERNAL FUNCTIONS ------
#' @importFrom grDevices col2rgb rgb
#' @importFrom graphics abline axTicks axis box image lines par points segments text
#' @importFrom stats as.formula coef dbeta dgamma dnbinom dnorm formula gaussian model.frame model.matrix model.response pnorm quantile rbeta rbinom rgamma rnbinom rnorm sd setNames terms
#' @importFrom utils head str tail
#' @import nimble
#'
### 1.1. ==== Change Variable Names to BUGS-Friendly Versions ====
#' @title Change Variable Names to BUGS-Friendly Versions
#'
#' @description Not all potential variable names used in R can be used in BUGS code without
#' producing a syntax error.  This function changes a vector of input names and ensures that
#' they are valid BUGS variable names.
#'
#' @param inName A character vector of variable names
#'
#' @return A character vector of BUGS-compliant variable names
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @keywords internal
#'
setBUGSVariableName <- function(inName) {
  outName <- tryCatch(as.character(inName), error = function(err) {
    stop(paste("invalid parameter name:", err, sep = " "))
  })
  # Remove any non-word characters in the name
  outName <- gsub("\\W+", "_", outName, perl = TRUE)
  # Replace any digit values with a text equivalent (this is to ensure that numbers in variable names aren't parsed incorrectly)
  outName <- gsub("0", "Zero", outName, fixed = TRUE)
  outName <- gsub("1", "One", outName, fixed = TRUE)
  outName <- gsub("2", "Two", outName, fixed = TRUE)
  outName <- gsub("3", "Three", outName, fixed = TRUE)
  outName <- gsub("4", "Four", outName, fixed = TRUE)
  outName <- gsub("5", "Five", outName, fixed = TRUE)
  outName <- gsub("6", "Six", outName, fixed = TRUE)
  outName <- gsub("7", "Seven", outName, fixed = TRUE)
  outName <- gsub("8", "Eight", outName, fixed = TRUE)
  outName <- gsub("9", "Nine", outName, fixed = TRUE)
  outName
}

### 1.2. ==== Create NIMBLE Model Structures for a mixture model====
#' @title Create NIMBLE Model Structures for a mixture model
#'
#' @description This function takes a set of model specifications for the three sub-models
#' of a mixture model and generates a set of structures
#' for initialisation of a NIMBLE model
#'
#' @param stateValModels A formula describing the regression relationship between the mean
#' of the response variable and the covariates for all mixture components, or a list of formulas
#' with each element giving the regression relationship between the mean
#' of the response variable and the covariates for each mixture component (ordered according to intercept of the
#' mixture component on the y-axis).  The response variable must be given on the
#' left-hand side of the formula.
#' @param stateProbModels A formula describing the regression relationship between the probability
#' of the existence of mixture components and the covariates, or a list
#' of formulas with each element giving the regression relationship between the probability
#' of the existence of a mixture component and its covariates (ordered
#' according to intercept of the mixture component on the y-axis).
#' @param statePrecModels A formula describing the regression relationship between the variances
#' in the mixture components and the covariates, or a list of
#' formulas with each mixture component giving the regression relationship between the variance in the
#' mixture component and the covariates for each mixture component (ordered according to
#' intercept of the mixture component on the y-axis).
#' @param inputData A data frame containing the covariate information and the response
#' variable.
#' @param numStates An integer scalar containing the number of mixture components.
#' If any of the \code{stateValModels}, \code{stateProbModels}, or \code{statePrecModels} parameters
#' is a list then \code{numStates} can be omitted and is therefore set to the length of the
#' list.
#' @param stateValError A description of the error distribution and link function to be used
#' in the model describing the response variable.  This can be from the \link[stats]{family}
#' specification or \code{character} scalar with the following possible values: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"}, or \code{"negbinomial"}.
#' @param setPriors A named list of prior distributions. Distribution are specified using character
#' strings. If sublists are not provided, values from the list are distributed to all sublist items
#' allowing to specify several priors at once. Sublist items are \code{"int"}, for specification of
#' priors on intercept parameters, and \code{"pred"}, from specification of priors on predictor/covariate
#' parameters. \code{"int"} is followed by \code{"1"} or \code{"2"} marking priors for the first
#' intercept and all the other intercepts respectively. For full structure of the list see default
#' values. Prior \code{"stateVal$Int2"} should allow only positive values to ensure distinctness of
#' states.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item{\code{modelText}}{ A character scalar containing the text of the model specification in
#' the NIMBLE BUGS dialect}
#' \item{\code{modelCode}}{ A \code{nimbleCode} object containing the model specification}
#' \item{\code{constants}}{ A list containing the constants to be provided to NIMBLE}
#' \item{\code{data}}{ A list containing the data to be provided to NIMBLE}
#' \item{\code{errorModel}}{ A factor containing the error model used for the specification of
#' the error distribution for the response variable}
#' \item{\code{linkFunction}}{ A factor containing the link function used in the specification
#' of the mixture model}
#' \item{\code{initialValues}}{ A list containing potential initial values for each of the stochastic
#' nodes in the NIMBLE model specification}
#' }
#'
#' @author Joseph D. Chipperfield \email{joechip90@@googlemail.com}, Adam Klimes
#' @keywords internal
#'
mixglmSpecification <- function(
  stateValModels,
  stateProbModels,
  statePrecModels,
  inputData,
  numStates = NULL,
  stateValError = gaussian,
  setPriors = list(
    stateVal = list(
      int1 = "dnorm(0.0, 0.001)",
      int2 = "dgamma(0.001, 0.001)",
      pred = "dnorm(0.0, 0.001)"),
    stateProb = list(
      int2 = "dnorm(0.0, 0.001)",
      pred = "dnorm(0.0, 0.001)"),
    statePrec = list(
      int = "dnorm(0.0, 0.001)",
      pred = "dnorm(0.0, 0.001)"))
) {
  # Small helper function to test whether a variable is a formula
  is.formula <- function(inVal){
    inherits(inVal, "formula")
  }
  # The supported error distributions
  supportedError <- c("gaussian", "gamma", "beta", "negbinomial") #, "betabinomial"
  # The suported link functions
  supportedLink <- c("identity", "log", "logit", "probit", "cloglog")
  # Sanity test the error distribution for the state error
  inStateValError <- stateValError
  inLinkFunction <- NA
  if(is.function(inStateValError)) {
    # If it is a function then call it to retrieve the family object
    inStateValError <- stateValError()
  }
  if(is.factor(inStateValError)) {
    # If it is a factor then convert it to a string
    inStateValError <- as.character(inStateValError)
  }
  if(is.character(inStateValError)) {
    # If it is a string then use the default link function associated with the specified family
    if(tolower(inStateValError[1]) == "gaussian") {
      inStateValError <- factor("gaussian", levels = supportedError)
      inLinkFunction <- factor("identity", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "gamma") {
      inStateValError <- factor("gamma", levels = supportedError)
      inLinkFunction <- factor("log", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "beta") {
      inStateValError <- factor("beta", levels = supportedError)
      inLinkFunction <- factor("logit", levels = supportedLink)
    } else if(tolower(inStateValError[1]) == "negbinomial") {
      inStateValError <- factor("negbinomial", levels = supportedError)
      inLinkFunction <- factor("log", levels = supportedLink)
#    } else if(tolower(inStateValError[1]) == "betabinomial") {
#      inStateValError <- factor("betabinomial", levels = supportedError)
#      inLinkFunction <- factor("logit", levels = supportedLink)
    } else {
      stop("selected error family is not supported")
    }
  } else if(inherits(inStateValError, "family")) {
    # If it is a "family" object then use the set error distribution and link function
    inLinkFunction <- factor(tolower(inStateValError$link), levels = supportedLink)
    inStateValError <- factor(tolower(inStateValError$family), levels = supportedError)
    if(is.na(inStateValError)) {
      stop("selected error family is not supported")
    }
    if(is.na(inLinkFunction)) {
      stop("selected link function is not supported")
    }
  } else {
    stop("error family specification is of invalid type")
  }
  inStateValModels <- stateValModels
  inStateProbModels <- stateProbModels
  inStatePrecModels <- statePrecModels
  inNumStates <- 1
  if(is.null(numStates)) {
    # If the maximum number of states is not set then find the largest model component to set the number
    # of states from
    inNumStates <- as.integer(max(c(
      ifelse(is.list(inStateValModels), length(inStateValModels), 1),
      ifelse(is.list(inStateProbModels), length(inStateProbModels), 1),
      ifelse(is.list(inStatePrecModels), length(inStatePrecModels), 1)
    )))
  } else {
    inNumStates <- tryCatch(as.integer(numStates), error = function(err) {
      stop("invalid entry for the number of states: ", err)
    })
  }
  if(inNumStates <= 0) {
    # Ensure that the number of states is greater than zero
    stop("invalid entry for the number of states: values equal or less than zero given")
  }
  # Define a function to facilitate the recycling and type specification of the model specification lists
  recycleModelFormulae <- function(modelInput, numStates, allowNull = FALSE) {
    inModelInput <- modelInput
    if((allowNull && is.null(modelInput)) || is.formula(modelInput)) {
      # If the model input is a formula or NULL then recycle that value to form a list
      inModelInput <- replicate(numStates, modelInput, simplify = FALSE)
    } else {
      # Otherwise force the input to be a list
      inModelInput <- tryCatch(as.list(modelInput), error = function(err) {
        stop("invalid entry for the model specification list: ", err)
      })
    }
    # Ensure that the list is of at least length one
    if(length(inModelInput)) {
      inModelInput <- lapply(X = inModelInput[((1:numStates - 1) %% length(inModelInput)) + 1], FUN = function(curEntry, allowNull) {
        outVal <- NULL
        if(is.null(curEntry)) {
          if(!allowNull) {
            stop("error encountered during the processing of the model specification list: NULL entries encountered")
          }
        } else {
          outVal <- tryCatch(as.formula(curEntry), error = function(err) {
            stop("error encountered during the processing of the model specification list: ", err)
          })
        }
        outVal
      }, allowNull = allowNull)
    } else {
      stop("invalid entry for the model specification list: list is empty")
    }
    inModelInput
  }
  # Recycle the model formulae so that they are of the correct length and type
  inStateValModels <- recycleModelFormulae(inStateValModels, inNumStates, FALSE)
  inStateProbModels <- recycleModelFormulae(inStateProbModels, inNumStates, FALSE)
  inStatePrecModels <- recycleModelFormulae(inStatePrecModels, inNumStates, TRUE)
  # Retrieve the list of formula strings for each of the models
  formulaStrings <- matrix(unlist(lapply(X = c(inStateValModels, inStateProbModels, inStatePrecModels), FUN = function(curForm) {
    outText <- NA
    if(!is.null(curForm)) {
      outText <- Reduce(paste, deparse(curForm))
    }
    outText
  })), ncol = 3, dimnames = list(NULL, c("stateVal", "stateProb", "statePrec")))
  # Retrieve the names of any response variables mentioned in any of the models
  respVariables <- unique(as.vector(gsub("\\s*~.*$", "", formulaStrings, perl = TRUE)))
  respVariables <- respVariables[!is.na(respVariables) & respVariables != ""]
  if(length(respVariables) != 1) {
    stop("invalid entry for the response variable: only one variable name must be present on the left-hand side of the formulae")
  }
  # Retrieve the names of the predictor variables mentioned in any of the models
  predVariables <- unique(gsub("^.*~\\s*", "", formulaStrings, perl = TRUE))
  predVariables <- predVariables[!is.na(predVariables) & predVariables != ""]
  predVariables <- unlist(lapply(X = predVariables, FUN = function(curVars) {
    strsplit(curVars, "\\s*\\+\\s*", perl = TRUE)[[1]]
  }))
  # Create a formula with the entrie set of predictor variables
  fullFormula <- as.formula(paste(respVariables, "~", paste(predVariables, collapse = " + "), sep = " "))
  # Retrieve the model matrix
  modelMatrix <- tryCatch(model.matrix(fullFormula, model.frame(fullFormula, inputData, na.action = NULL)), error = function(err) {
    stop("error thrown during construction of the model matrix: ", err)
  })
  # Remove the intercept term in the model matrix
  modelMatrix <- modelMatrix[, colnames(modelMatrix) != "(Intercept)", drop = FALSE]
  # Retrieve the model response variable
  respValues <- model.response(model.frame(fullFormula, inputData, na.action = NULL))
  numTrials <- NULL
  if(!is.null(dim(respValues))) {
    if(length(dim(respValues)) == 2) {
      if(ncol(respValues) == 1) {
        # Treat a single column in a same way as a dimensionless
        respValues <- tryCatch(as.numeric(respValues), error = function(err) {
          stop("error thrown during processing of response variable: ", err)
        })
        numTrials <- rep(NA, length(respValues))
      } else if(ncol(respValues) == 2) {
        # Check to ensure that the values in these columns take appropriate values
        if(any(respValues < 0.0)) {
          stop("a two-column response variable must have only positive values")
        }
        numTrials <- tryCatch(as.numeric(apply(X = as.matrix(respValues), FUN = sum, na.rm = TRUE, MARGIN = 1)), error = function(err) {
          stop("error thrown during processing of the number of trials: ", err)
        })
        respValues <- tryCatch(as.numeric(respValues[, 1]), error = function(err) {
          stop("error thrown during processing of response variable: ", err)
        })
      } else {
        stop("response variable can only have one or two columns")
      }
    } else {
      stop("response variable has invalid dimension structure")
    }
  } else {
    # Process the response variable
    respValues <- tryCatch(as.numeric(respValues), error = function(err) {
      stop("error thrown during processing of response variable: ", err)
    })
    numTrials <- rep(NA, length(respValues))
  }
  if(length(respValues) != nrow(modelMatrix)) {
    stop("response variable does not have the same length as the predictor variables")
  }
  # Convert the name of the response variable to something BUGS-friendly
  respVariablesBUGS <- setBUGSVariableName(respVariables)
  # Rename the covariates from the model matrix to something BUGS-friendly
  covariatesBUGS <- setBUGSVariableName(colnames(modelMatrix))
  # Very occasionally the conversion of the covariate names results in duplicates: this line protects from the possibility
  covariatesBUGS <- setNames(
    ifelse(duplicated(covariatesBUGS), paste(covariatesBUGS, setBUGSVariableName(as.character(1:length(covariatesBUGS))), sep = "_"), covariatesBUGS),
    colnames(modelMatrix))
  colnames(modelMatrix) <- covariatesBUGS
  # Set the link prefix and suffix for the state value model
  linkPrefix <- ""
  linkSuffix <- ""
  if(inLinkFunction != "identity") {
    linkPrefix <- paste(inLinkFunction, "(", sep = "")
    linkSuffix <- ")"
  }
  # Function to retrieve the covariate names from each sub-model formula
  getCovNames <- function(curFormula, inputData, covariatesBUGS) {
    outNames <- NA
    if(!is.na(curFormula)) {
      # Retrieve the model covariates
      modelCovars <- tryCatch(colnames(model.matrix(as.formula(curFormula), model.frame(as.formula(curFormula), inputData, na.action = NULL))), error = function(err) {
        stop("error thrown during construction of sub-model matrix: ", err)
      })
      modelCovars <- modelCovars[modelCovars != "(Intercept)"]
      outNames <- c("intercept", covariatesBUGS[modelCovars])
      if(anyNA(outNames)) {
        stop("error thrown during construction of sub-model matrix: covariates present in a sub-model that are not found in the full model")
      }
    }
    outNames
  }
  # Check priors
  if (!all(grepl("^d.+\\(.*)$", unlist(setPriors))))
    stop("unexpected prior specification")
  # Helper function which takes call of the list and modify/amends its items on all levels
  callModify <- function(oldCall, mods){
    callObj <- eval(oldCall)
    # Recursive function which modifies items of a list (and add new ones)
    recMod <- function(target, mod){
      for (i in names(mod)){
        target[i] <- if (is.list(mod[[i]]))
          list(recMod(target[[i]], mod[[i]])) else mod[[i]]
      }
      target
    }
    # Recursive function which copies structure of the list propagating items to lower levels
    recStr <- function(structure, content){
      for (i in names(structure)){
        contentItem <- if (is.null(names(content))) content else content[[i]]
        structure[i] <- if (is.list(structure[[i]]))
          list(recStr(structure[[i]], contentItem)) else contentItem
      }
      structure
    }
    recStr(callObj, recMod(callObj, mods))
  }
  # Prepare full object specifying them priors
  inSetPriors <- callModify(formals(mixglmSpecification)$setPriors, setPriors)
  # Initialise a vector to store potential initial values for the model
  initialValues <- as.character(c())
  # If the model has more than one state (99% of the times this function will be called) then create the relevant model strings
  if(inNumStates > 1) {
    # Create a matrix of strings containing the model specification
    modelStrings <- t(sapply(X = 1:inNumStates, FUN = function(curState, formulaStrings, inputData, covariatesBUGS, stateValError, linkFunction) {
      # Create a string version of the current state number
      stateString <- tolower(setBUGSVariableName(curState))
      # Retrieve the covariate names for each of the sub-models
      stateValCovs <- getCovNames(formulaStrings[curState, 1], inputData, covariatesBUGS)
      stateValCovs_nonIntercept <- stateValCovs[stateValCovs != "intercept"]
      stateProbCovs <- getCovNames(formulaStrings[curState, 2], inputData, covariatesBUGS)
      statePrecCovs <- getCovNames(formulaStrings[curState, 3], inputData, covariatesBUGS)
      # Prepare auxiliary vectors of priors
      auxPriorsPrec <- c(inSetPriors$statePrec$int, rep(inSetPriors$statePrec$pred, length(statePrecCovs) - 1))
      auxPriorsProb <- c(inSetPriors$stateProb$int2, rep(inSetPriors$stateProb$pred, length(stateProbCovs) - 1))
      # Set the prior text for the state value model
      priorStateVal <- paste(
        paste("\t# Set priors for the state variable value model for state ", stateString, sep = ""),
        if (length(stateValCovs_nonIntercept) > 0) paste("\t", stateValCovs_nonIntercept, "_stateVal[", curState, "] ~ ", inSetPriors$stateVal$pred, sep = "", collapse = "\n"),
        # The intercept of the first state has a normal prior.  All other states are forced to have positive priors in order to ensure that the
        # state labels are ordered and that MCMC doesn't just do state relabelling.
        paste("\tintercept_stateVal[", curState, "] ~ ", ifelse(curState == 1, inSetPriors$stateVal$int1, inSetPriors$stateVal$int2), sep = ""),
        sep = "\n")
      # Set the prior text for the state probability model
      priorStateProb <- "\t# The first state probability model is a baseline model so has no parameters"
      if(curState > 1) {
        priorStateProb <- paste(
          paste("\t# Set priors for the state probability model for state ", stateString, sep = ""),
          paste("\t", stateProbCovs, "_stateProb[", curState, "] ~", auxPriorsProb, sep = "", collapse = "\n"),
          sep = "\n")
      }
      # Set the prior text for the state precision model
      # Intially assume a simple multiplier model
      priorStatePrec <- paste(
        paste("\t# Set priors for the state precision model for state ", stateString, " (simple multiplier model)", sep = ""),
        paste("\tlinStateProb_statePrec[", curState, "] ~", inSetPriors$statePrec$pred, sep = ""),
        paste("\tintercept_statePrec[", curState, "] ~", inSetPriors$statePrec$int, sep = ""),
        sep = "\n")
      if(!is.na(formulaStrings[curState, 3])) {
        # If a formula has been specified for the precision model then use a linear sub-model instead
        priorStatePrec <- paste(
          paste("\t# Set priors for the state precision model for state ", stateString, sep = ""),
          paste("\t", statePrecCovs, "_statePrec[", curState, "] ~", auxPriorsPrec, sep = "", collapse = "\n"),
          sep = "\n")
      }
      # Set the model specification text for the state value model
      likelihoodStateVal <- paste(
        paste("\t\t# Set the model specification for the state value for state ", stateString, sep = ""),
        paste("\t\t", linkPrefix, "linStateVal[dataIter, ", curState, "]", linkSuffix, " <- ",
          ifelse(curState > 1, paste("sum(intercept_stateVal[1:", curState, "])", sep = ""), "intercept_stateVal[1]"), " * intercept[dataIter]", if (length(stateValCovs_nonIntercept > 0)) "+",
          if (length(stateValCovs_nonIntercept > 0)) paste(stateValCovs_nonIntercept, "_stateVal[", curState, "] * ", stateValCovs_nonIntercept, "[dataIter]", sep = "", collapse = " + "), sep = ""),
        sep = "\n")
      # Set the model specification text for the state probability model
      likelihoodStateProb <- paste("\t\t# Set the model specification for the state probability model for state ", stateString, sep = "")
      if(curState > 1) {
        likelihoodStateProb <- paste(
          likelihoodStateProb,
          paste("\t\tlog(linStateProb[dataIter, ", curState, "]) <- ", paste(stateProbCovs, "_stateProb[", curState, "] * ", stateProbCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
          sep = "\n")
      } else {
        likelihoodStateProb <- paste(
          likelihoodStateProb,
          paste("\t\tlinStateProb[dataIter, ", curState, "] <- 1.0", sep = ""),
          sep = "\n")
      }
      # Set the model specification text for the state precision model
      # Initially assume a simple multiplier model
      likelihoodStatePrec <- paste(
        paste("\t\t# Set the model specification for the state precision model for state ", stateString, sep = ""),
        paste("\t\tlog(linStatePrec[dataIter, ", curState, "]) <- intercept_statePrec[", curState, "] + linStateProb_statePrec[", curState, "] * linStateProb[dataIter, ", curState, "] / sum(linStateProb[dataIter, 1:numStates])", sep = ""),
        sep = "\n")
      if(!is.na(formulaStrings[curState, 3])) {
        # If a formula has been specified for the precision model then use a linear sub-model instead
        likelihoodStatePrec <- paste(
          paste("\t\t# Set the model specification for the state precision model for state ", stateString, sep = ""),
          paste("\t\tlog(linStatePrec[dataIter, ", curState, "]) <- ", paste(statePrecCovs, "_statePrec[", curState, "] * ", statePrecCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
          sep = "\n")
      }
      # Set an output vector with the appropriate names
      setNames(
        c(priorStateVal, priorStateProb, priorStatePrec, likelihoodStateVal, likelihoodStateProb, likelihoodStatePrec),
        c("priorValModel", "priorProbModel", "priorPrecModel", "likelihoodValModel", "likelihoodProbModel", "likelihoodPrecModel"))
    }, formulaStrings = formulaStrings, inputData = inputData, covariatesBUGS = covariatesBUGS, stateValError = as.character(inStateValError), linkFunction = as.character(inLinkFunction)))
    # Assign the error distribution for the ecosystem state value model
    errorStrings <- paste("\t\t# Set the error specification model for the state value sub-model", switch(as.character(inStateValError),
      "gaussian" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dnormStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
      sep = ""),
      "gamma" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dgammaStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
        sep = ""),
      "beta" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dbetaStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
        sep = ""),
      "negbinomial" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dnegbinStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates])",
        sep = ""),
      "betabinomial" = paste(
        "\t\t", respVariablesBUGS, "[dataIter] ~ dbetabinStateValueMembership(linStateVal[dataIter, 1:numStates], linStatePrec[dataIter, 1:numStates], linStateProb[dataIter, 1:numStates], numTrials[dataIter])",
        sep = "")
    ), sep = "\n")
    # Create a vector of potential initial values for the model parameters
    initialValues <- unlist(lapply(X = 1:inNumStates, FUN = function(curState, formulaStrings, inputData, covariatesBUGS) {
      # Retrieve the covariate names for each of the sub-models
      stateValCovs <- getCovNames(formulaStrings[curState, 1], inputData, covariatesBUGS)
      stateValCovs_nonIntercept <- stateValCovs[stateValCovs != "intercept"]
      stateProbCovs <- getCovNames(formulaStrings[curState, 2], inputData, covariatesBUGS)
      statePrecCovs <- getCovNames(formulaStrings[curState, 3], inputData, covariatesBUGS)
      # Initialise a vecotr of output values
      outValuesNames <- c("intercept", stateValCovs_nonIntercept)
      outValues <- setNames(c(
        rnorm(length(stateValCovs_nonIntercept), 0.0, 1.0),
        ifelse(curState > 1, abs(rnorm(1, 0.0, 1.0)), rnorm(1, 0.0, 1.0))
      ), paste(outValuesNames, "_stateVal[", curState, "]", sep = ""))
      if(curState > 1) {
        # Add the the probability sub-model parameters if the current state is greater than 1
        outValues <- c(outValues, setNames(
          rnorm(length(stateProbCovs), 0.0, 1.0),
          paste(stateProbCovs, "_stateProb[", curState, "]", sep = "")
        ))
      }
      # Add the precision sub-model parameters
      if(is.na(formulaStrings[curState, 3])) {
        outValues <- c(outValues, setNames(
          rnorm(2, 0.0, 1.0),
          paste(c("intercept_statePrec", "linStateProb_statePrec"), "[", curState, "]", sep = "")
        ))
      } else {
        outValues <- c(outValues, setNames(
          rnorm(length(statePrecCovs), 0.0, 1.0),
          paste(statePrecCovs, "_statePrec[", curState, "]", sep = "")
        ))
      }
    }, formulaStrings = formulaStrings, inputData = inputData, covariatesBUGS = covariatesBUGS))
  } else {
    # If the model has only one state then completely remove the multinomial state latent state components
    # Retrieve the covariate names for each of the sub-models
    stateValCovs <- getCovNames(formulaStrings[1, 1], inputData, covariatesBUGS)
    statePrecCovs <- getCovNames(formulaStrings[1, 3], inputData, covariatesBUGS)
    # Create a matrix of model text
    auxPriorsVal <- c(inSetPriors$stateVal$int1, rep(inSetPriors$stateVal$pred, length(stateValCovs) - 1))
    auxPriorsPrec <- c(inSetPriors$statePrec$int, rep(inSetPriors$statePrec$pred, length(statePrecCovs) - 1))
    modelStrings <- matrix(nrow = 1, dimnames = list(NULL, c("priorValModel", "priorProbModel", "priorPrecModel", "likelihoodValModel", "likelihoodProbModel", "likelihoodPrecModel")), data = c(
      # Set the prior for the state value model
      paste(
        "\t# Set priors for the state variable value model",
        paste("\t", stateValCovs, "_stateVal ~", auxPriorsVal, sep = "", collapse = "\n"),
      sep = "\n"),
      # Set the prior for the state probability model: there is no state probability model because there is only one state
      "\t# There is no state probability model because there is only one state",
      # Set the prior for the state precision model
      ifelse(is.na(formulaStrings[1, 3]),
        paste("\t# Set priors for the state precision model\n\tintercept_statePrec ~", inSetPriors$statePrec$int),
        paste(
          "\t# Set priors for the state precision model",
          paste("\t", statePrecCovs, "_statePrec ~", auxPriorsPrec, sep = "", collapse = "\n"),
        sep = "\n")
      ),
      # Set the model specification text for the state value model
      paste(
        "\t\t# Set the model specification for the state value",
        paste("\t\t", linkPrefix, "linStateVal[dataIter]", linkSuffix, " <- ", paste(stateValCovs, "_stateVal * ", stateValCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
        sep = "\n"),
      # Set the model specification text for the state probability model: there is no state probability model because there is only one state
      "\t\t# There is no state probability model because there is only one state",
      # Set teh model specification text for the state precision model
      ifelse(is.na(formulaStrings[1, 3]),
        "\t\t# Set the model specification for the state precision model\n\t\tlog(linStatePrec[dataIter]) <- intercept_statePrec",
        paste(
          "\t\t# Set the model specification for the state precision model",
          paste("\t\tlog(linStatePrec[dataIter]) <- ", paste(statePrecCovs, "_statePrec * ", statePrecCovs, "[dataIter]", sep = "", collapse = " + "), sep = ""),
          sep = "\n")
      )
    ))
    # Assign the error distribution for the ecosystem state precision model
    errorStrings <- paste("\t\t# Set the error specification model for the state value sub-model", switch(as.character(inStateValError),
      "gaussian" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dnorm(linStateVal[dataIter], linStatePrec[dataIter])", sep = ""),
      "gamma" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dgamma(mean = linStateVal[dataIter], sd = pow(linStatePrec[dataIter], -0.5))", sep = ""),
      "beta" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dbeta(mean = linStateVal[dataIter], sd = pow(linStatePrec[dataIter], -0.5))", sep = ""),
#      "negbinomial" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dnegbin(\n\t\t\t1.0 - linStateVal[dataIter] * linStatePrec[dataIter], \n\t\t\tlinStateVal[dataIter] * linStateVal[dataIter] * linStatePrec[dataIter] / (1.0 - linStateVal[dataIter] * linStatePrec[dataIter]))", sep = ""),
      "negbinomial" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dnegbin(\n\t\t\tlinStateVal[dataIter] * linStatePrec[dataIter], \n\t\t\tlinStateVal[dataIter] * linStateVal[dataIter] * linStatePrec[dataIter] / (1.0 - linStateVal[dataIter] * linStatePrec[dataIter]))", sep = ""),
      "betabinomial" = paste("\t\t", respVariablesBUGS, "[dataIter] ~ dbetabin(mean = linStateVal[dataIter], prec = linStatePrec[dataIter], size = numTrials[dataIter])")
    ), sep = "\n")
    # Create a vector of potential initial values for the model parameters
    initialValues <- setNames(rnorm(length(stateValCovs), 0.0, 1.0), paste(stateValCovs, "_stateVal", sep = ""))
    if(is.na(formulaStrings[1, 3])) {
      initialValues <- c(initialValues, setNames(rnorm(1, 0.0, 1.0), "intercept_statePrec"))
    } else {
      initialValues <- c(initialValues, setNames(rnorm(length(statePrecCovs), 0.0, 1.0), paste(statePrecCovs, "_statePrec", sep = "")))
    }
  }
  # Create NIMBLE model code
  modelText <- paste(
    "nimbleCode({",
    "\n\t#### PRIOR SPECIFICATION ####",
    "\n\t# Set priors for the ecosystem state value sub-model ----",
    paste(modelStrings[, "priorValModel"], collapse = "\n"),
    "\n\t# Set priors for the ecosystem state probability sub-model ----",
    paste(modelStrings[, "priorProbModel"], collapse = "\n"),
    "\n\t# Set priors for the ecosystem state precision sub-model ----",
    paste(modelStrings[, "priorPrecModel"], collapse = "\n"),
    "\n\t#### MODEL STRUCTURE ####",
    "\n\t# Iterate over each data point",
    "\tfor(dataIter in 1:numData) {",
    "\n\t\t# Set model structure for the ecosystem state value sub-model ----",
    paste(modelStrings[, "likelihoodValModel"], collapse = "\n"),
    "\n\t\t# Set model structure for the ecosystem state probability sub-model ----",
    paste(modelStrings[, "likelihoodProbModel"], collapse = "\n"),
    "\n\t\t# Set model structure for the ecosystem state precision sub-model ----",
    paste(modelStrings[, "likelihoodPrecModel"], collapse = "\n"),
    "\n\t\t# Set the error model for the ecosystem state value sub-model ----",
    errorStrings,
    "\t}",
    "})",
    sep = "\n")
  # Parse the NIMBLE model code to create a code object
  modelCode <- eval(parse(text = modelText))
  # Create a set of constants to use in the model
  modelConstants <- append(list(
    numData = nrow(modelMatrix),
    numStates = inNumStates,
    intercept = rep(1.0, nrow(modelMatrix))
  ), as.list(as.data.frame(modelMatrix)))
  if(as.character(inStateValError) == "betabinomial") {
    # Add the number of trials if the betabinomial error distribution is being used
    modelConstants <- append(modelConstants, list(numTrials = numTrials))
  }
  # Restructrue the initial values as a list
  vectorNames <- unique(gsub("\\[.*$", "", names(initialValues), perl = TRUE))
  initialValuesList <- setNames(lapply(X = vectorNames, FUN = function(curVecName, initialValues, inNumStates) {
    # Initialise an output vector
    outVec <- rep(0.0, inNumStates)
    if(inNumStates > 1) {
      # Retrieve the covairates with the current vector
      curCovs <- initialValues[grepl(paste("^", curVecName, "\\[\\d+\\]$", sep = ""), names(initialValues), perl = TRUE)]
      # Fill in the covariate values in the respective places
      outVec[as.integer(gsub("\\]$", "", gsub("^.*\\[", "", names(curCovs), perl = TRUE), perl = TRUE))] <- as.double(curCovs)
    } else {
      outVec <- initialValues[curVecName]
    }
    outVec
  }, initialValues = initialValues, inNumStates = inNumStates), vectorNames)
  # Create a set of data to use in the model
  modelData <- list(response = respValues)
  names(modelData) <- respVariablesBUGS
  list(
    modelText = modelText,
    modelCode = modelCode,
    constants = modelConstants,
    data = modelData,
    errorModel = inStateValError,
    linkFunction = inLinkFunction,
    initialValues = initialValuesList
  )
}

### 1.2. ==== Check initial values and generate new ones if necessary ====
#' @title Check initial values for mixglm and generate new ones if necessary
#'
#' @description This function checks initial values for mixglm and if they are illegitimate, it tries
#' to generate new ones.
#'
#' @param stateValModels A formula describing the regression relationship between the mean
#' of the response variable and the covariates for all mixture components, or a list of formulas
#' with each element giving the regression relationship between the mean
#' of the response variable and the covariates for each mixture component (ordered according to intercept of the
#' mixture component on the y-axis).  The response variable must be given on the
#' left-hand side of the formula.
#' @param statePrecModels A formula describing the regression relationship between the variances
#' in the mixture components and the covariates, or a list of
#' formulas with each mixture component giving the regression relationship between the variance in the
#' mixture component and the covariates for each mixture component (ordered according to
#' intercept of the mixture component on the y-axis).
#' @param inputData A data frame containing the covariate information.
#' @param setInit List of initial values. It is assumed that for each parameter, there is a vector
#' of initial values with length corresponding to \code{Nstates}
#' @param errorModel A string specifying the error distribution to be used
#' in the mixglm. This can be: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"} or \code{"negbinomial"}.
#' @param linkFunction A string specifying the link function. This can be: \code{"identity"},
#' \code{"log"}, \code{"logit"}, \code{"probit"}, or \code{"cloglog"}.
#' @param Nstates An integer specifying number of components in the mixture.
#' @param genNew A logical value specifying if new initial values should be generated even if
#' \code{setInit} contains legitimate values.
#' @param genFn A function used to generate new initial values. It has one argument which is
#' the component rank and it should return one value. Note that for intercept of precision of beta
#' distribution, \code{genFn(1) + 6} is used.
#'
#' @return A list of initial values
#'
#' @author Adam Klimes
#' @keywords internal
#'
checkInit <- function(stateValModels, statePrecModels, inputData, setInit,
  errorModel, linkFunction, Nstates, genNew = FALSE, genFn = function(x) rnorm(1, 0, 2)){
  evalModMat <- function(modForm, formType, inputData, initVal, state){
    if (is.list(modForm)) modForm <- modForm[[state]]
    if (formType == "stateVal") initVal$intercept_stateVal <- cumsum(initVal$intercept_stateVal)
    modMat <- model.matrix(modForm, inputData)
    initNames <- paste0("_", formType, "$")
    IDaux <- grep(initNames, names(initVal))
    initAux <- initVal[IDaux]
    names(initAux) <- sub(initNames, "", names(initAux))
    names(initAux)[names(initAux) == "intercept"] <- "(Intercept)"
    matchID <- match(colnames(modMat), names(initAux))
    initAux <- initAux[matchID]
    init <- unlist(lapply(initAux, function(x, state) x[state], state))
    list(val = as.vector(modMat %*% init), ID = IDaux[matchID])
  }
  auxID <- grep("^intercept", names(setInit))
  initVal <- c(setInit[auxID], setInit[-auxID])
  checkFun <- switch(as.character(errorModel),
                     gaussian = function(M, P) TRUE,
                     gamma = function(M, P) M > 0 & P > 0,
                     beta = function(M, P) {A <- M * M * (1 - M) * P - M; B <- M * (1 - M)^2 * P + M - 1; A > 0 & B > 0},
                     negbinomial = function(M, P) {r <- M * M * P / (1 - P * M); r >= 0 & r <= 1})
  invlink <- switch(as.character(linkFunction), identity = function(x) x, log = exp,
                    logit = function(x) exp(x)/(1+exp(x)), probit = pnorm, cloglog = function(x) 1 - exp(-exp(x)))
  genP <- if (errorModel == "beta") function(x) genFn(x) + c(0, 6)[(x == 1) + 1] else function(x) genFn(x)
  updateItem <- function(x, state, fn) Map(function(x, ID, state, fn) {x[state] <- fn(ID); x}, x, 1:length(x), list(state), list(fn))
  for (state in 1:Nstates){
    i <- 1
    if (genNew) {
      M <- evalModMat(stateValModels, "stateVal", inputData, initVal, state)
      P <- evalModMat(statePrecModels, "statePrec", inputData, initVal, state)
      initVal[M$ID] <- updateItem(initVal[M$ID], state, genFn)
      initVal[P$ID] <- updateItem(initVal[P$ID], state, genP)
      if (state > 1) initVal[grep("_stateProb$", names(initVal))] <- updateItem(initVal[grep("_stateProb$", names(initVal))], state, genFn)
    }
    repeat{
      M <- evalModMat(stateValModels, "stateVal", inputData, initVal, state)
      P <- evalModMat(statePrecModels, "statePrec", inputData, initVal, state)
      if (all(checkFun(invlink(M$val), exp(P$val)))) break
      i <- i + 1
      if (i > 999) stop("Legitimate initial values fail to be fould automatically. You can provide them as an argument. Standardization of predictors might also help.")
      initVal[M$ID] <- updateItem(initVal[M$ID], state, genFn)
      initVal[P$ID] <- updateItem(initVal[P$ID], state, genP)
    }
  }
  initVal
}

## 2. ------ DEFINE MODELLING FUNCTIONS ------

### 2.1. ==== Simulate Instances of a mixture model ====
#' @title Simulate Instances of a mixture model
#'
#' @description This function takes a set of model specifications for the three sub-models
#' of a mixture model and generates a simulation from this
#' model specification
#'
#' @param numSims The number of simulations to draw from the mixture model.
#' @param stateValModels A formula describing the regression relationship between the mean
#' of the response variable and the covariates for all mixture components, or a list of formulas
#' with each element giving the regression relationship between the mean
#' of the response variable and the covariates for each mixture component (ordered according to intercept of the
#' mixture component on the y-axis).  The response variable must be given on the
#' left-hand side of the formula.
#' @param stateProbModels A formula describing the regression relationship between the probability
#' of the existence of mixture components and the covariates, or a list
#' of formulas with each element giving the regression relationship between the probability
#' of the existence of a mixture component and its covariates (ordered
#' according to intercept of the mixture component on the y-axis).
#' @param statePrecModels A formula describing the regression relationship between the variances
#' in the mixture components and the covariates, or a list of
#' formulas with each mixture component giving the regression relationship between the variance in the
#' mixture component and the covariates for each mixture component (ordered according to
#' intercept of the mixture component on the y-axis).
#' @param inputData A data frame containing the covariate information and the response
#' variable.
#' @param numStates An integer scalar containing the number of components in the mixture.
#' If any of the \code{stateValModels}, \code{stateProbModels}, or \code{statePrecModels} parameters
#' is a list then \code{numStates} can be omitted and is therefore set to the length of the
#' list.
#' @param stateValError A description of the error distribution and link function to be used
#' in the model describing the ecosystem state value.  This can be from the \link[stats]{family}
#' specification or \code{character} scalar with the following possible values: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"}, \code{"negbinomial"}, or \code{"betabinomial"}.
#' @param coefficientValues A list containing the values of the coefficients to use in the
#' simulation.
#'
#' @return A list containg a vector of simulated values for each stochastic node.  In addition the
#' following elements are appended to the list:
#' \itemize{
#' \item{\code{linStateVal}}{ A matrix containing the predicted response variable values at each row
#' of the input covariate data.frame.  Each column represents the predicted response variable value
#' for each mixture component}
#' \item{\code{linStatePrec}}{ A matrix containing the predicted precision of the response variable
#' at each row of the input covariate data.frame.  Each column represents the predicted
#' precision of the response variable value for each mixture component}
#' \item{\code{linStateProb}}{ A matrix containing the probability of each mixture component existing
#' at each row of the input covariate data.frame.  Each column represents the probability of each
#' mixture component existing}
#' }
#'
#' @author Joseph D. Chipperfield, \email{joechip90@@googlemail.com}
#' @examples {
#' dat <- data.frame(x = rnorm(200), y = rnorm(200))
#' cfValues <- list(
#'   "intercept_stateVal" = c(-1.0, 0.5, 0.5),
#'   "x_stateVal" = c(-0.05, 0.05, 0.1),
#'   "intercept_stateProb" = c(NA, -1.0, 1.0),
#'   "intercept_statePrec" = c(0.5, 1.0, 1.5)
#'  )
#'  mixglmSimulation(
#'    numSims = 1,
#'    stateValModels = y ~ x,
#'    stateProbModels = ~ 1,
#'    statePrecModels = ~ 1,
#'    inputData = dat,
#'    numStates = 2,
#'    stateValError = gaussian,
#'    coefficientValues = cfValues
#'  )
#' }
#' @export
#'
mixglmSimulation <- function(
  numSims,
  stateValModels,
  stateProbModels,
  statePrecModels,
  inputData,
  numStates = NULL,
  stateValError = gaussian,
  coefficientValues = NULL
) {
  # Ensure that the coefficient values are correctly specified
  inCoefficientValues <- list()
  if(!is.null(coefficientValues)) {
    inCoefficientValues <- tryCatch(as.list(coefficientValues), error = function(err) {
      stop("error encountered during processing of the coefficient value list: ", err)
    })
  }
  # Ensure that the number of simulations input is correctly specified
  inNumSims <- tryCatch(as.integer(numSims)[1], error = function(err) {
    stop("error encountered during processing of the number of simulations: ", err)
  })
  if(inNumSims <= 0) {
    stop("error encountered during processing of the number of simulations: number of simulations requested is less than 1")
  }
  # Create a NIMBLE model specification
  modelSpecification <- mixglmSpecification(stateValModels, stateProbModels, statePrecModels, inputData, numStates, stateValError)
  # Initialise the data with some dummy (but plausable) values for the relevant error model.  This is so that the NIMBLE
  # model is fully initialised (avoids some error messages and ensures NIMBLE starts from a valid state)
  modelSpecification$data[[1]] <- rep(switch(as.character(modelSpecification$errorModel),
    "gaussian" = 0.0,
    "gamma" = 0.1,
    "beta" = 0.5,
    "negbinomial" = 0,
    "betabinomial" = 0
    ), length(modelSpecification$data[[1]]))
  # If coefficient values have been provided then use those in the model initialisation
  if(length(inCoefficientValues) > 0 && !is.null(names(inCoefficientValues))) {
    # Overwrite the list of initial values with the input coefficients (where appropriate)
    modelSpecification$initialValues <- setNames(lapply(X = names(modelSpecification$initialValues), FUN = function(curName, curInits, inputInits) {
      # Initialise the current output
      outValues <- curInits[[curName]]
      if(curName %in% names(inputInits) && length(outValues) > 0) {
        outValues <- ifelse(is.na(inputInits[[curName]][(1:length(outValues) - 1) %% length(inputInits) + 1]), outValues, inputInits[[curName]][(1:length(outValues) - 1) %% length(inputInits) + 1])
      }
      outValues
    }, curInits = modelSpecification$initialValues, inputInits = inCoefficientValues), names(modelSpecification$initialValues))
    # Check the variable names to ensure that they are present in the model
    isInModel <- names(inCoefficientValues) %in% names(modelSpecification$initialValues)
    if(any(!isInModel)) {
      warning("some coefficient names are not present in the model: ", print(names(inCoefficientValues)[!isInModel], collapse = ", "))
    }
  }
  # Initialise the NIMBLE model
  modelObject <- nimbleModel(modelSpecification$modelCode, constants = modelSpecification$constants, data = modelSpecification$data, inits = modelSpecification$initialValues)
  # Specify a function to simulate data (of a particular data node)
  singleSimulation <- function(curDataName, modelObject, modelSpecification) {
    # Simulate the data for the model (overwrites the previously set dummy data)
    simulate(modelObject, names(modelSpecification$data), includeData = TRUE)
    # Retrieve the simulated data
    modelObject[[curDataName]]
  }
  # Initialise an output object with the simulations of the data
  outObject <- setNames(lapply(X = names(modelSpecification$data), FUN = function(curDataName, modelObject, modelSpecification, numSims) {
    # Call simulation function for each requested simulation
    outVec <- do.call(cbind, replicate(numSims, singleSimulation(curDataName, modelObject, modelSpecification), simplify = FALSE))
    # Ensure that the output object has two dimensions
    if(is.null(dim(outVec)) || length(dim(outVec)) <= 1) {
      dim(outVec) <- c(length(outVec), 1)
    }
    outVec
  }, modelObject = modelObject, modelSpecification = modelSpecification, numSims = inNumSims), names(modelSpecification$data))
  # Retrieve the linear predictors for the value and precision sub-models also
  outObject <- append(outObject, list(
    linStateVal = modelObject[["linStateVal"]],
    linStatePrec = modelObject[["linStatePrec"]]
  ))
  # Ensure consistent dimensions of the linear predictor objects
  if(is.null(dim(outObject$linStateVal)) || length(dim(outObject$linStateVal)) == 1) {
    dim(outObject$linStateVal) <- c(length(outObject$linStateVal), 1)
  }
  if(is.null(dim(outObject$linStatePrec)) || length(dim(outObject$linStatePrec)) == 1) {
    dim(outObject$linStatePrec) <- c(length(outObject$linStatePrec), 1)
  }
  # Add the probability model outputs if they are missing
  if(is.null(modelObject$linStateProb)) {
    outObject <- append(outObject, list(linStateProb = rep(1.0, length(outObject$linStateVal))))
    dim(outObject$linStateProb) <- c(length(outObject$linStateVal), 1)
  # Otherwise normalise the probability values and add it to the output list
  } else {
    outObject <- append(outObject, list(linStateProb = t(apply(X = modelObject$linStateProb, FUN = function(curRow) {
      curRow / sum(curRow, na.rm = TRUE)
    }, MARGIN = 1))))
  }
  outObject
}

### 2.2. ==== Fit the mixture model ====
#' @title Fit the mixture model
#'
#' @description This function takes a set of model specifications for the three sub-models
#' of the mixture model and fits them to a data set.
#'
#' @param stateValModels A formula describing the regression relationship between the mean
#' of the response variable and the covariates for all mixture components, or a list of formulas
#' with each element giving the regression relationship between the mean
#' of the response variable and the covariates for each mixture component (ordered according to intercept of the
#' mixture component on the y-axis).  The response variable must be given on the
#' left-hand side of the formula.
#' @param stateProbModels A formula describing the regression relationship between the probability
#' of the existence of mixture components and the covariates, or a list
#' of formulas with each element giving the regression relationship between the probability
#' of the existence of a mixture component and its covariates (ordered
#' according to intercept of the mixture component on the y-axis).
#' @param statePrecModels A formula describing the regression relationship between the variances
#' in the mixture components and the covariates, or a list of
#' formulas with each mixture component giving the regression relationship between the variance in the
#' mixture component and the covariates for each mixture component (ordered according to
#' intercept of the mixture component on the y-axis).
#' @param inputData A data frame containing the covariate information and the response
#' variable.
#' @param numStates An integer scalar containing the number of components in the mixture.
#' If any of the \code{stateValModels}, \code{stateProbModels}, or \code{statePrecModels} parameters
#' is a list then \code{numStates} can be omitted and is therefore set to the length of the
#' list.
#' @param stateValError A description of the error distribution and link function to be used
#' in the mixture model.  This can be from the \link[stats]{family}
#' specification or \code{character} scalar with the following possible values: \code{"gaussian"},
#' \code{"gamma"}, \code{"beta"}, or \code{"negbinomial"}.
#' @param mcmcIters An integer scalar providing the number of post-burn-in samples to draw from the
#' MCMC sampler (per chain).
#' @param mcmcBurnin An integer scalar providing the number of MCMC samples to use for the
#' adaption or burn-in portion of the process (per chain).
#' @param mcmcChains An integer scalar giving the number of MCMC chains to use.
#' @param mcmcThin An integer scalar giving the thinning frequency in the MCMC chains.  For example,
#' a value of \code{4} results in every fourth values being retained.
#' @param setPriors A named list of prior distributions. Distribution are specified using character
#' strings. If sublists are not provided, values from the list are distributed to all sublist items
#' allowing to specify several priors at once. Sublist items are \code{"int"}, for specification of
#' priors on intercept parameters, and \code{"pred"}, from specification of priors on predictor
#' parameters. \code{"int"} is followed by \code{"1"} or \code{"2"} marking priors for the first
#' intercept and all the other intercepts respectively. For full structure of the list see default
#' values. Prior \code{"stateVal$Int2"} should allow only positive values to ensure distinctness of
#' mixture components.
#' @param setInit list of initial values which overwrites generated ones.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item{\code{mcmcSamples}}{ An \link[coda]{mcmc} object if \code{mcmcChains == 1} or \link[coda]{mcmc.list}
#' object if \code{mcmcChains > 1} containing the sampled values}
#' \item{\code{compiledModel}}{ An R interface object containing an interface for the compiled NIMBLE model.
#' This is the same output as the \link[nimble]{compileNimble} function applied to the model object}
#' \item{\code{modelText}}{ A character scalar containing the text of the model specification in
#' the NIMBLE BUGS dialect}
#' \item{\code{modelCode}}{ A \code{nimbleCode} object containing the model specification}
#' \item{\code{constants}}{ A list containing the constants provided to NIMBLE}
#' \item{\code{data}}{ A list containing the data provided to NIMBLE}
#' \item{\code{errorModel}}{ A factor containing the error model used for the specification of
#' the error distribution for the response variable}
#' \item{\code{linkFunction}}{ A factor containing the link function used in the specification
#' of the mixture model}
#' \item{\code{initialValues}}{ A list containing potential initial values used for each of the stochastic
#' nodes in the NIMBLE model specification}
#' }
#'
#' @author Joseph D. Chipperfield \email{joechip90@@googlemail.com}, Adam Klimes
#' @examples \dontrun{
#' set.seed(10)
#' n <- 200
#' x <- rnorm(n)
#' group <- rbinom(n, 1, 0.5)
#' y <- rnorm(n, 1 + 0.5 * x * c(-1, 1)[group + 1], 0.1)
#' plot(y ~ x)
#' dat <- data.frame(x, y)
#'
#' mod <- mixglm(
#'   stateValModels = y ~ x,
#'   stateProbModels = ~ x,
#'   statePrecModels = ~ x,
#'   inputData = dat,
#'   numStates = 2)
#' plot(mod)}
#' @export
#'
mixglm <- function(
  stateValModels,
  stateProbModels,
  statePrecModels,
  inputData,
  numStates = NULL,
  stateValError = gaussian,
  mcmcIters = 10000,
  mcmcBurnin = 5000,
  mcmcChains = 4,
  mcmcThin = 1,
  setPriors = list(
    stateVal = list(
      int1 = "dnorm(0.0, 0.001)",
      int2 = "dgamma(0.001, 0.001)",
      pred = "dnorm(0.0, 0.001)"),
    stateProb = list(
      int2 = "dnorm(0.0, 0.001)",
      pred = "dnorm(0.0, 0.001)"),
    statePrec = list(
      int = "dnorm(0.0, 0.001)",
      pred = "dnorm(0.0, 0.001)")),
  setInit = NULL
) {
  # Create a NIMBLE model specification
  modelSpecification <- mixglmSpecification(stateValModels, stateProbModels, statePrecModels, inputData, numStates, stateValError, setPriors)
  # Change initial values if provided
  modelSpecification$initialValues <- if (!is.null(setInit)) setInit else checkInit(stateValModels, statePrecModels, inputData, modelSpecification$initialValues, modelSpecification$errorModel, modelSpecification$linkFunction, modelSpecification$constants$numStates)
  modelObject <- nimbleModel(modelSpecification$modelCode, constants = modelSpecification$constants, data = modelSpecification$data, inits = modelSpecification$initialValues)
  # Build the MCMC object and compile it
  # varsToMonitor <- c(modelObject$getVarNames(), "linStateVal", "linStatePrec")
  # if (grepl("linStateProb", modelSpecification$modelText)) varsToMonitor <- c(varsToMonitor, "linStateProb")
  # Monitor only parameters
  # varsToMonitor <- modelObject$getVarNames()
  # varsToMonitor <- varsToMonitor[!varsToMonitor %in% c(names(modelSpecification$data), "linStateVal", "linStatePrec", "linStateProb")]
  varsToMonitor <- names(modelSpecification$initialValues)
  varsToMonitor <- varsToMonitor[varsToMonitor %in% modelObject$getVarNames()]
  mcmcObject <- buildMCMC(modelObject, enableWAIC = TRUE, monitors = varsToMonitor)
  mcmcObjectCompiled <- compileNimble(mcmcObject, modelObject)
  # Run the MCMC
  mcmcOutput <- runMCMC(mcmcObjectCompiled$mcmcObject, niter = mcmcIters, nburnin = mcmcBurnin, thin = mcmcThin, nchains = mcmcChains, WAIC = TRUE, samplesAsCodaMCMC = TRUE)
  # Structure the compiled model, the MCMC samples, and the model specification into a list
  out <- append(list(mcmcSamples = mcmcOutput, compiledModel = mcmcObjectCompiled), modelSpecification)
  class(out) <- "mixglm"
  out
}

## 3. ------ DEFINE GENERAL MODEL METHODS ------

### 3.1. ==== Plot results of a mixture model ====
#' @title Plot results of a mixture model
#'
#' @description This function plots results of a mixture model on the
#' current graphical device.
#'
#' @param x an object of class "mixglm"
#' @param form formula, such as y ~ pred, specifying variables to be plotted. By
#' default the first predictor is plotted on x axis.
#' @param byChains logical value indicating whether to plot mixture components for each chain
#' @param transCol logical value indicating usage of transparent colours
#' @param addWAIC logical value indication display of WAIC in upper right corner of the plot
#' @param setCol vector of colours to be used for mixture components
#' @param drawAxes logical value indicating whether values should be marked on axes
#' @param SDmult scalar multiplying visualized standard deviation (to make
#' lines for small standard deviation visible)
#' @param xlab a label for the x axis, defaults to predictor name.
#' @param ylab a label for the y axis, defaults to response name.
#' @param ... additional arguments passed to \link[graphics]{plot}
#'
#' @author Adam Klimes
#' @examples \dontrun{
#' set.seed(10)
#' n <- 200
#' x <- rnorm(n)
#' group <- rbinom(n, 1, 0.5)
#' y <- rnorm(n, 1 + 0.5 * x * c(-1, 1)[group + 1], 0.1)
#' plot(y ~ x)
#' dat <- data.frame(x, y)
#'
#' mod <- mixglm(
#'   stateValModels = y ~ x,
#'   stateProbModels = ~ x,
#'   statePrecModels = ~ x,
#'   inputData = dat,
#'   numStates = 2)
#' plot(mod)}
#' @export
#'
plot.mixglm <- function(x, form = NULL, byChains = TRUE, transCol = TRUE,
  addWAIC = FALSE, setCol = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666"),
  drawAxes = TRUE, SDmult = 1, xlab = NULL, ylab = NULL, ...) {
  # input check
  form <- if (is.null(form))
    formula(paste(names(x$data), "~", names(x$constants)[4])) else formula(form)
  svar <- labels(terms(form))
  if (length(svar) > 1) stop("Specify only one predictor in 'form'")
  if (!all.vars(form)[1] %in% names(x$data))
    stop("Response specified in 'form' does not match 'mod'")
  if (!svar %in% names(x$constants))
    stop("Predictor specified in 'form' does not match 'mod'")
  if (!is.logical(byChains)) stop("'byChains' has to be logical")
  if (!is.logical(transCol)) stop("'transCol' has to be logical")
  if (!is.logical(addWAIC)) stop("'addWAIC' has to be logical")
  if (!is.logical(drawAxes)) stop("'drawAxes' has to be logical")
  if (!is.numeric(SDmult)) stop("'SDmult' has to be numeric")
  #_
  resp <- x$data[[1]]
  dat <- data.frame(x$data,
    x$constants[vapply(x$constants, length, FUN.VALUE = 1) == length(resp)])
  if (is.null(xlab)) xlab <- svar
  if (is.null(ylab)) ylab <- names(x$data)
  auxRange <- max(resp) - min(resp)
  invlink <- switch(as.character(x$linkFunction), identity = function(x) x, log = exp,
    logit = function(x) exp(x)/(1+exp(x)), probit = pnorm, cloglog = function(x) 1 - exp(-exp(x)))
  plot(form, data = dat, ylim = c(min(resp) - 0.05 * auxRange, max(resp) + 0.3 * auxRange),
       yaxs = "i", axes = FALSE, ann = FALSE, ...)
  usr <- par("usr")
  maxY <- max(resp) + 0.05 * auxRange
  axis(1, labels = c("", ""), at = c(2*usr[1]-usr[2], 2*usr[2]-usr[1]))
  axis(2, labels = c("", ""), at = c(2*usr[3]-usr[4], maxY), lwd.ticks = 0)
  axis(1, labels = xlab, at = mean(usr[1:2]), line = 2, tick = FALSE)
  if (drawAxes) {
    axis(1)
    yaxis <- head(axTicks(2), -1)
    axis(2, labels = yaxis, at = yaxis, las = 2)
  }
  axis(2, labels = 0:1, at = max(resp) + c(0.1, 0.25) * auxRange, las = 2)
  abline(h = maxY, lwd = 3)
  abline(h = max(resp) + 0.1 * auxRange, lty = 2)
  abline(h = max(resp) + 0.25 * auxRange, lty = 2)
  axis(2, labels = "Probability", at = max(resp) + 0.175 * auxRange, line = 2, tick = FALSE)
  axis(2, labels = ylab, at = mean(range(resp)), line = 2, tick = FALSE)
  if (addWAIC) text(par("usr")[2] - (par("usr")[2] - par("usr")[1]) * 0.2,
    max(resp) + 0.175 * auxRange, paste("WAIC:", round(x$mcmcSamples$WAIC$WAIC, 1)))
  parsTab <- summary(x, byChains = byChains, absInt = TRUE, digit = NULL)
  parsTab <- lapply(parsTab, function(x) {x[is.na(x)] <- 0; x})
  auxLines <- function(parsChain, dat, mod){
    nstates <- mod$constants$numStates
    xx <- seq(min(dat[, svar]), max(dat[, svar]), length.out = 100)
    ind <- NULL
    cNames <- rownames(parsChain)
    if (nstates > 1) {
      ind <- paste0("[", 1:nstates, "]")
      probInt <- parsChain[paste0("intercept_stateProb", ind), "mean"]
    }
    valInt <- parsChain[paste0("intercept_stateVal", ind), "mean"]
    precInt <- parsChain[paste0("intercept_statePrec", ind), "mean"]
    valCov <- if (paste0(svar, "_stateVal", ind[1]) %in% cNames)
      parsChain[paste0(svar, "_stateVal", ind), "mean"] else rep(0, nstates)
    precCov <- if (paste0(svar, "_statePrec", ind[1]) %in% cNames)
      parsChain[paste0(svar, "_statePrec", ind), "mean"] else rep(0, nstates)
    probCov <- if (paste0(svar, "_stateProb", ind[1]) %in% cNames)
      parsChain[paste0(svar, "_stateProb", ind), "mean"] else rep(0, nstates)
    if (nstates > 1) {
      probVals <- as.matrix(data.frame(Map(function(int, cov) exp(int + cov * xx), probInt, probCov)))
      probVals[is.na(probVals)] <- 1
      probVals <- probVals / rowSums(probVals)
      probVals[is.nan(probVals)] <- 1
      }
    for (i in 1:nstates){
      cols <- setCol[i]
      if (nstates > 1) {
        lines(xx, max(resp) + 0.1 * auxRange + probVals[, i] * 0.15 * auxRange,
          col = setCol[i], lwd = 3)
        if (transCol) {
          rgbVec <- col2rgb(cols)[, 1]
          cols <- rgb(rgbVec[1], rgbVec[2], rgbVec[3],
            alpha = 40 + probVals[, i] * 215, maxColorValue = 255)
        }
      }
      sdVals <- 1 / sqrt(exp(precInt[i] + precCov[i] * xx))
      yEst <- do.call(invlink, list(valInt[i] + valCov[i] * xx))
      uci <- do.call(invlink, list(valInt[i] + valCov[i] * xx + sdVals * SDmult))
      lci <- do.call(invlink, list(valInt[i] + valCov[i] * xx - sdVals * SDmult))
      yEstSel <- yEst < maxY
      uciSel <- uci < maxY
      lciSel <- lci < maxY
      segments(head(xx[yEstSel], -1), head(yEst[yEstSel], -1), x1 = tail(xx[yEstSel], -1),
        y1 = tail(yEst[yEstSel], -1), col = if (length(cols) > 1) cols[yEstSel] else cols, lwd = 3)
      lines(xx[uciSel], uci[uciSel], col = setCol[i], lty = 2, lwd = 1)
      lines(xx[lciSel], lci[lciSel], col = setCol[i], lty = 2, lwd = 1)
    }
  }
  lapply(parsTab, auxLines, dat, x)
  invisible()
}

### 3.2. ==== Summary of a mixture model ====
#' @title Summarize a mixture model
#'
#' @description This function calculates posterior quantiles of parameters of
#' a mixture model across all chains or for each chain separately
#'
#' @param object an object of class "mixglm"
#' @param byChains logical value indicating if the summary should be calculated
#' for each chain separately
#' @param digit integer specifying the number of decimal places to be used.
#' Use \code{"NULL"} for no rounding.
#' @param absInt logical value indicating if intercepts for state values should
#' be absolute (by default, they represent differences)
#' @param randomSample integer specifying how many random samples from posterior
#' distribution to take instead of summary. Use \code{"NULL"} for summary.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Returns a list with one item being data.frame of quantiles of
#' posterior distributions of parameters. If \code{byChains} is \code{TRUE},
#' a data.frame is returned for each chain.
#'
#' @author Adam Klimes
#' @examples \dontrun{
#' set.seed(10)
#' n <- 200
#' x <- rnorm(n)
#' group <- rbinom(n, 1, 0.5)
#' y <- rnorm(n, 1 + 0.5 * x * c(-1, 1)[group + 1], 0.1)
#' plot(y ~ x)
#' dat <- data.frame(x, y)
#'
#' mod <- mixglm(
#'   stateValModels = y ~ x,
#'   stateProbModels = ~ x,
#'   statePrecModels = ~ x,
#'   inputData = dat,
#'   numStates = 2)
#' summary(mod)
#'
#' # random samples of parameters from posterior distribution
#' summary(mod, randomSample = 5)}
#' @export
#'
summary.mixglm <- function(object, byChains = FALSE, digit = 4L,
  absInt = FALSE, randomSample = NULL, ...){
  # input check
  if (!is.logical(byChains)) stop("'byChains' has to be logical")
  if (!(is.numeric(digit) | is.null(digit))) stop("'digit' has to be NULL or numeric")
  if (!is.logical(absInt)) stop("'absInt' has to be logical")
  if (!(is.numeric(randomSample) | is.null(randomSample)))
    stop("'randomSample' has to be NULL or numeric")
  #_
  mcmcList <- if (is.list(object$mcmcSamples$samples))
    object$mcmcSamples$samples else list(object$mcmcSamples$samples)
  varsSamples <- lapply(mcmcList,
    function(x) x[, !grepl(paste0("^lifted|^linState|^", names(object$data)), colnames(x))])
  if (!byChains) varsSamples <- list(do.call(rbind, varsSamples))
  sepInt <- function(samp){
    scol <- grepl("intercept_stateVal", colnames(samp))
    samp[, scol] <- t(apply(samp, 1, function(x, scol) cumsum(x[scol]), scol))
    samp
  }
  if (absInt) varsSamples <- lapply(varsSamples, sepInt)
  if (is.null(randomSample)){
    auxSummary <- function(x)
      c(mean = mean(x), sd = sd(x), quantile(x, c(0.025,0.25,0.75,0.975), na.rm = TRUE))
    out <- lapply(varsSamples, function(x) t(apply(x, 2, auxSummary)))
  } else {
    out <- lapply(varsSamples,
      function(x) t(x[sample(1:nrow(varsSamples[[1]]), randomSample), , drop = FALSE]))
  }
  if (!is.null(digit)) out <- lapply(out, round, digit)
  out
}

### 3.3. ==== Predict resilience metrics based on a mixture model ====
#' @title Predict from mixglm
#'
#' @description This function calculates stability curves and resilience metrics based
#' on a mixture model
#'
#' @param object an object of class "mixglm"
#' @param newdata dataframe of predictor values.
#'   If it contains column with response variable, \code{obsDat} is returned.
#'   If not provided, prediction is done for modelled data.
#' @param samples number of samples to take along the response variable
#' @param threshold number from 0 to 1 denoting how pronounced stable states
#' should be marked as considerable
#' @param ... further arguments passed to or from other methods.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item{\code{sampledResp}}{ A numeric vector of samples along response variable}
#' \item{\code{probCurves}}{ A data frame of stability curves for each observation}
#' \item{\code{tipStable}}{ A list of data.frames for each observation with all
#' stable states and tipping points, their scaled probability density and response
#' variable value, whether they are stable states, and whether their satisfy the threshold}
#' \item{\code{obsDat}}{ A data frame containing values of response variable, potential energy,
#' distance to closest tipping point and stable state for each observation,
#' potential depth defined as a difference between potential energy of the
#' observation and the closest tipping point}
#' }
#'
#' @author Adam Klimes
#' @examples \dontrun{
#' set.seed(10)
#' n <- 200
#' x <- rnorm(n)
#' group <- rbinom(n, 1, 0.5)
#' y <- rnorm(n, 1 + 0.5 * x * c(-1, 1)[group + 1], 0.1)
#' plot(y ~ x)
#' dat <- data.frame(x, y)
#'
#' mod <- mixglm(
#'   stateValModels = y ~ x,
#'   stateProbModels = ~ x,
#'   statePrecModels = ~ x,
#'   inputData = dat,
#'   numStates = 2)
#' predict(mod)
#'
#' # predict for newdata
#' pre <- predict(mod, newdata = data.frame(x = c(0.5, 1), y = c(0.3, 0.2)))
#' plot(pre$sampledResp, pre$probCurves$X1, type = "l")
#' plot(pre$sampledResp, pre$probCurves$X2, type = "l")}
#' @export
#'
predict.mixglm <- function(object, newdata = NULL, samples = 1000, threshold = 0.0, ...){
  # input check
  if (!(is.null(newdata) | is.data.frame(newdata)))
    stop("'newdata' has to be NULL or data.frame")
  #_
  respVal <- NULL
  if (names(object$data) %in% colnames(newdata)) {
    respVal <- newdata[, names(object$data)]
    newdata <- newdata[, -which(colnames(newdata) == names(object$data)), drop = FALSE]
  }
  if (is.null(newdata)) {
    respVal <- object$data[[1]]
    newdata <- as.data.frame(object$constants[-(1:3)])
  }
  if (is.null(respVal)) respVal <- object$data[[1]]
  form <- formula(paste(names(object$data), "~", colnames(newdata)[1]))
  slices <- sliceMixglm(object, form, value = newdata, byChains = FALSE, doPlot = FALSE, samples = samples)
  tipStableAll <- getMinMax(slices, threshold)
  tipStable <- tipStableAll$tipStable[[1]]
  probCurve <- data.frame(tipStableAll$matsSt[[1]])
  resp <- slices$resp
  auxFn <- function(x){
    dist <- abs(resp - x)
    which(dist == min(dist))[1]
  }
  potentEn <- unlist(Map(function(x, y) -x[y]+1, probCurve, vapply(respVal, auxFn, FUN.VALUE = 1)))
  getClosest <- function(x, targetResp, state = 1, getVar = "respDist"){
    xSel <- x[x$state == state & x$catSt == 1, ]
    distT <- abs(xSel$resp - targetResp)
    if (length(distT) == 0) distT <- NA
    cbind(xSel[distT == min(distT), ], respDist = min(distT))[, getVar]
  }
  obsDat <- NULL
  if (!is.null(respVal)){
    obsDat <- data.frame(
      respVal = respVal,
      potentEn = potentEn,
      distToTip = unlist(Map(getClosest, tipStable, respVal, state = 0)),
      distToState = unlist(Map(getClosest, tipStable, respVal, state = 1)),
      potentialDepth = 1 - potentEn - unlist(Map(getClosest, tipStable, respVal, state = 0, getVar = "probDens")))
  }
  list(sampledResp = resp,
              probCurves = probCurve,
              tipStable = tipStable,
              obsDat = obsDat)
}

### 3.4. ==== Extract model coefficients from a mixture model ====
#' @title Extract model coefficients from a mixture model
#'
#' @description Posterior mean values from mixglm
#'
#' @param object an object of class "mixglm"
#' @param digit integer specifying the number of decimal places to be used.
#' Use \code{"NULL"} for no rounding.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A List of matrices of posterior mean values for each parameter
#'
#' @author Adam Klimes
#' @export
#'
coef.mixglm <- function(object, digit = NULL, ...){
  s <- summary(object, digit = digit)[[1]]
  getPars <- function(state, s){
    auxGetPars <- function(type, state, s){
      stateIndex <- if (Nstates == 1) "$" else paste0("\\[", state, "\\]$")
      aux <- s[grep(paste0("_state", type, stateIndex), rownames(s)), "mean", drop = FALSE]
      rownames(aux) <- sub(paste0("_state", type, stateIndex), "", rownames(aux))
      data.frame(ID = rownames(aux), aux)
    }
    types <- if (Nstates == 1) c("Val", "Prec") else c("Val", "Prec", "Prob")
    parsPerType <- lapply(types, auxGetPars, state, s)
    mergeAll <- function(x, y) merge(x, y, by = "ID", all = TRUE)
    out <- Reduce(mergeAll, parsPerType)
    rownames(out) <- out$ID
    colnames(out) <- c("ID", types)
    intID <- which(rownames(out) == "intercept")
    out[c(intID, (1:nrow(out))[-intID]), -1, drop = FALSE]
  }
  Nstates <- object$constants$numStates
  res <- lapply(1:Nstates, getPars, s)
  names(res) <- paste0("State", 1:object$constants$numStates)
  res
}

### 3.5. ==== Print a mixture model ====
#' @title Print a mixture model
#'
#' @description This function prints summary information about mixglm
#'
#' @param x an object of class "mixglm"
#' @param ... further arguments passed to or from other methods.
#'
#' @return Invisibly returns \code{x}
#'
#' @author Adam Klimes
#' @export
#'
print.mixglm <- function(x, ...){
  WAIC <- x$mcmcSamples$WAIC
  cat("Multinomial Ecosystem State Model\n")
  cat("WAIC:", WAIC$WAIC, "\n")
  cat("pWAIC:", WAIC$pWAIC, "\n")
  cat("Posterior mean values:\n")
  print(coef(x, 3))
  invisible(x)
}

### 3.6. ==== Compactly Display the Structure of a mixglm object ====
#' @title Compactly Display the Structure of a mixglm object
#'
#' @description Compactly display the internal structure of a mixglm object
#'
#' @param object an object of class "mixglm"
#' @param max.level integer giving maximal level of nesting of displayed structures
#' @param give.attr logical indicating if attributes should be shown
#' @param ... arguments passed to \link[base]{NextMethod}
#'
#' @author Adam Klimes
#' @export
#'
str.mixglm <- function(object, max.level = 2, give.attr = FALSE, ...){
  cat("List of", length(object), "\n")
  NextMethod(str, object, max.level = max.level, give.attr = give.attr, ...)
}

## 4. ------ DEFINE HELPER FUNCTIONS ------

### 4.1. ==== Find positions of local minima in a vector ====
#' @title Find positions of local minima in a vector
#'
#' @description Finds position of all local minima in a vector including start
#'   and end point. For flat minima (identical subsequent values), it denotes middle point
#'
#' @param x numeric vector
#' @param extremes logical indicating if first and last values should be considered
#'
#' @return A vector of positions of minima in x
#'
#' @author Adam Klimes
#' @keywords internal
#'
findMin <- function(x, extremes = TRUE){
  dfXin <- diff(x)
  seqCount <- diff(c(0, which(dfXin != 0), length(x)))
  Nflat <- rep(seqCount, seqCount) - 1
  xClear <- x[c(TRUE,  dfXin != 0)]
  dfX <- diff(xClear)
  loc <- which(diff(sign(dfX)) == 2) + 1
  if (extremes){
    if (dfX[1] > 0) loc <- c(1, loc)
    if (tail(dfX, 1) < 0) loc <- c(loc, length(xClear))
  }
  inLoc <- seq_along(x)[c(TRUE, dfXin != 0)][loc]
  inLoc[inLoc %in% which(dfXin == 0)] <-
    0.5 * Nflat[inLoc[inLoc %in% which(dfXin == 0)]] + inLoc[inLoc %in% which(dfXin == 0)]
  inLoc
}

### 4.2. ==== Find local minima and maxima in mixglm ====
#' @title Find local minima and maxima in mixglm
#'
#' @description Finds local minima and maxima in slices of stability landscape
#'   from mixglm which represent multimodality over
#'   specified threshold.
#'
#' @param slices slices of stability landscape from sliceMixglm() function
#' @param threshold numerical value specifying threshold for marking minima
#'   and maxima. Marked maxima/minima differ in scaled probability density by
#'   threshold value.
#'
#' @return A list containing the following components:
#' \itemize{
#' \item{\code{tipStable}}{ A list of lists of dataframes with
#'   tipping points (state == 0) and stable states (state == 1), categorized
#'   based on satifying the threshold (catSt == 1), with their scaled [0-1]
#'   probability density and response variable value}
#' \item{\code{mats}}{ An array of stability landscapes}
#' \item{\code{matsSt}}{ A list of scaled stability landscapes}}
#'
#' @author Adam Klimes
#' @keywords internal
#'
getMinMax <- function(slices, threshold = 0.0){
  if (!is.numeric(threshold)) stop("'threshold' has to be numeric")
  if (threshold > 1) warning("maximum value for 'threshold' is 1")
  getMin <- function(x, resp, extremes = TRUE, inv = FALSE) {
    id <- findMin(x, extremes = extremes)
    respOut <- resp[id] * (1 - id %% 1) + resp[min(id + 1, length(resp))] * id %% 1
    probDens <- x[id]
    if (length(respOut) == 0) {
      respOut <- NA
      probDens <- NA
    }
    if (inv) probDens <- -probDens
    data.frame(probDens = probDens, resp = respOut)
  }
  combStates <- function(tip, state, threshold){
    comb <- rbind(cbind(tip, state = 0), cbind(state, state = 1))
    combSort <- comb[order(comb$resp), ]
    if (nrow(combSort) < 3 | threshold == 0) catSt <- 1 else {
      dif <- (abs(diff(-combSort$probDens)) > threshold) + 0
      groups <- lapply(0:sum(dif), function(x) which(x == c(0,cumsum(dif))))
      sdf <- sign(diff(-combSort$probDens))
      catDf <- -diff(c(-1, sdf[dif == 1], 1))/2
      mins <- vapply(groups[catDf == -1],
        function(x) x[which(-combSort$probDens[x] == min(-combSort$probDens[x]))], FUN.VALUE = 1)
      maxs <- vapply(groups[catDf == 1],
        function(x) x[which(-combSort$probDens[x] == max(-combSort$probDens[x]))], FUN.VALUE = 1)
      auxChange <- function(x, signC) {
        comp <- diff(-combSort$probDens[c(min(x)-1, max(x))])
        if (length(comp) == 0) comp <- -1
        if ((signC * comp) > 0) NULL else range(x)
      }
      addInc <- lapply(groups[catDf == 0 & c(-1, sdf[dif == 1]) == 1], auxChange, signC = 1)
      addDec <- lapply(groups[catDf == 0 & c(-1, sdf[dif == 1]) == -1], auxChange, signC = -1)
      minsAll <- sort(c(mins, unlist(lapply(addInc, "[", 2)), unlist(lapply(addDec, "[", 1))))
      maxsAll <- sort(c(maxs, unlist(lapply(addInc, "[", 1)), unlist(lapply(addDec, "[", 2))))
      catSt <- rep(0, nrow(combSort))
      catSt[c(minsAll, maxsAll)] <- 1
    }
    cbind(combSort, catSt = catSt)
  }
  mat <-  slices[[1]]
  mats <- array(do.call(c, mat), dim = c(dim(mat[[1]]), length(mat)))
  scaleCurve <- function(x) (x - min(x)) / max((x - min(x)))
  matsSt <- apply(mats, 3, function(x) apply(x, 2, scaleCurve), simplify = FALSE)
  maxs <- lapply(matsSt, function(y) apply(y, 2, findMin, extremes = FALSE, simplify = FALSE))
  mins <- lapply(matsSt, function(y) apply(y, 2, function(x) findMin(-x), simplify = FALSE))
  tipPoints <- lapply(matsSt, function(y) apply(y, 2, getMin, slices$resp, extremes = FALSE))
  stableStates <- lapply(matsSt, function(y) apply(-y, 2, getMin, slices$resp, inv = TRUE))
  tipStable <- Map(Map, list(combStates), tipPoints, stableStates, threshold)
  list(tipStable = tipStable, mats = mats, matsSt = matsSt)
}

### 4.3. ==== Plot slice from mixglm ====
#' @title Plot slice from mixglm
#'
#' @description This function plots probability density for given predictor value
#'
#' @param mod an object of class "mixglm"
#' @param form formula with one predictor specifying which variables to plot
#' @param value numeric vector of values of the preditor specified by
#'   \code{"form"} where the slice is done or data.frame with values of predictors
#'   in named columns
#' @param byChains logical value indicating if slice should be done for each chain separately
#' @param xlab string used as label for x-axis
#' @param doPlot logical value indicating if plotting should be done
#' @param setCol vector of colours to be used for visualization of estimated states
#' @param plotEst logical value indicating if estimated mixture components should be visualized
#' @param xaxis logical value indicating if values should be marked on x-axis
#' @param addEcos logical value indicating if observations within \code{"ecosTol"}
#' from \code{"value"} should be visualized on the stability curve
#' @param ecosTol scalar specifying range of predictor from the \code{"value"}
#' to select observations to be visualized
#' @param samples scalar specifying number of samples to be taken along predictor
#' @param randomSample integer specifying how many random samples from posterior
#' distribution to take instead of summary. Use \code{"NULL"} for summary.
#' @param getParsD logical value indicating if the output should be parameters
#' of distributions instead of the slice. Requires \code{doPlot} to be \code{FALSE}
#' @param respIn mainly for internal use. Numeric vector specifying values along
#' response variable for which to compute probability density.
#'
#' @return Returns invisibly a list containing the following components:
#' \itemize{
#' \item{\code{chainX}}{ One or more items of a matrix of probability density
#' values along the response with a column for each \code{value}. One item is
#' returned if \code{byChains} is not \code{TRUE}, otherwise one item for each
#' chain. If \code{getParsD} is \code{TRUE}, parameters for probability density
#' distributions are returned instead.}
#' \item{\code{densFun}}{ If \code{getParsD} is \code{TRUE}, function for
#' calculation of probability density is returned.}
#' \item{\code{resp}}{ A vector of response values for which stability curves are
#' evaluated}}
#'
#' @author Adam Klimes
#' @examples \dontrun{
#' set.seed(10)
#' n <- 200
#' x <- rnorm(n)
#' group <- rbinom(n, 1, 0.5)
#' y <- rnorm(n, 1 + 0.5 * x * c(-1, 1)[group + 1], 0.1)
#' plot(y ~ x)
#' dat <- data.frame(x, y)
#'
#' mod <- mixglm(
#'   stateValModels = y ~ x,
#'   stateProbModels = ~ x,
#'   statePrecModels = ~ x,
#'   inputData = dat,
#'   numStates = 2)
#' sliceMixglm(mod, value = 1)
#'
#' # plot including observations around x == 1
#' sliceMixglm(mod, value = 1, addEcos = TRUE)
#'
#' # plot 5 random samples from posterior distribution
#' sliceMixglm(mod, value = 1, randomSample = 5, byChain = FALSE)}
#' @export
#'
sliceMixglm <- function(mod, form = NULL, value = 0, byChains = TRUE,
  xlab = names(mod$data), doPlot = TRUE,
  setCol = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"),
  plotEst = TRUE, xaxis = TRUE, addEcos = FALSE, ecosTol = 0.1, samples = 1000,
  randomSample = NULL, getParsD = FALSE, respIn = NULL){
  # input check
  if (!inherits(mod, "mixglm")) stop("'mod' must be an object of class 'mixglm'")
  form <- if (is.null(form))
    formula(paste(names(mod$data), "~", names(mod$constants)[4])) else formula(form)
  svar <- labels(terms(form))
  if (length(svar) > 1) stop("Specify only one predictor in 'form'")
  if (!all.vars(form)[1] %in% names(mod$data))
    stop("Response specified in 'form' does not match 'mod'")
  if (!svar %in% names(mod$constants))
    stop("Predictor specified in 'form' does not match 'mod'")
  if (!(is.numeric(value) | is.data.frame(value)))
    stop("'value' has to be a numeric or data.frame")
  if (is.data.frame(value)) if (any(!colnames(value) %in% names(mod$constants)))
    stop("colnames of 'value' do not match 'mod'")
  if (!is.logical(doPlot)) stop("'doPlot' has to be logical")
  if (!is.logical(plotEst)) stop("'plotEst' has to be logical")
  if (!is.logical(xaxis)) stop("'xaxis' has to be logical")
  if (!is.logical(addEcos)) stop("'addEcos' has to be logical")
  if (!is.numeric(ecosTol)) stop("'ecosTol' has to be numeric")
  if (!is.numeric(samples)) stop("'samples' has to be numeric")
  if (!is.logical(getParsD)) stop("'getParsD' has to be logical")
  if (getParsD & doPlot) stop("please set 'doPlot' to FALSE to use 'getParsD'")
  if (!is.numeric(respIn) & !is.null(respIn)) stop("'respIn' has to be numeric")
  #_
  resp <- mod$data[[1]]
  parsTab <- summary(mod, byChains = byChains, absInt = TRUE, digit = NULL, randomSample = randomSample)
  parsTab <- lapply(parsTab, function(x) {x[is.na(x)] <- 0; x})
  if (is.null(randomSample)) parsTab <- lapply(parsTab, function(x) x[, "mean", drop = FALSE])
  Nstates <- mod$constants$numStates
  invlink <- switch(as.character(mod$linkFunction), identity = function(x) x, log = exp,
    logit = function(x) exp(x)/(1+exp(x)), probit = pnorm, cloglog = function(x) 1 - exp(-exp(x)))
  xx <- seq(min(resp), max(resp), length.out = samples)
  if (mod$errorModel == "negbinomial") {
    xx <- seq(min(resp), max(resp), by = trunc((max(resp)-min(resp)+1) / samples) + 1)
    samples <- length(xx)
  }
  if (addEcos) {
    pred <- mod$constants[[svar]]
    xx <- c(xx, resp[abs(pred - value) < ecosTol])
  }
  if (!is.null(respIn)) xx <- respIn
  if (is.null(dim(value))) value <- matrix(value, length(value), dimnames = list(NULL, svar))
  value <- as.matrix(value)
  calcParsVal <- function(pars, value, mod){
    getPars <- function(curState, pars, value){
      pars <- as.matrix(pars)
      auxExtract <- function(toGet, curState, pars, value){
        stateIndex <- if (Nstates == 1) NULL else paste0("[", curState, "]")
        pos <- match(paste0(c("intercept", colnames(value)), "_", toGet, stateIndex), rownames(pars))
        sel <- !is.na(pos[-1])
        pars[pos[1], 1] + as.vector(value[, sel, drop = FALSE] %*% pars[pos[-1][sel], 1])
      }
      est <- auxExtract("stateVal", curState, pars, value)
      prec <- auxExtract("statePrec", curState, pars, value)
      prob <- auxExtract("stateProb", curState, pars, value)
      cbind(est = do.call(invlink, list(est)), sd = 1 / sqrt(exp(prec)), prob = prob)
    }
    parsVal <- vapply(1:Nstates, getPars, FUN.VALUE = array(0, dim = c(nrow(value), 3)), pars, value)
    parsVal[, "prob", 1] <- rep(0, nrow(value))
    parsVal[, "prob", ] <- exp(parsVal[, "prob", ]) / rowSums(exp(parsVal[, "prob", , drop = FALSE]))
    rownames(parsVal) <- paste0("value", 1:nrow(value))
    parsVal
  }
  makeDensFun <- function(){
    dfun <- switch(as.character(mod$errorModel),
                   gaussian = "dnorm",
                   gamma = "dgamma",
                   beta = "dbeta",
                   negbinomial = "dnbinom")
    i <- 1:Nstates
    strFun <- paste0("(", paste(paste0("probs[", i, "] * ", dfun,
      "(x, ", "parOne[", i, "], ", "parTwo[", i, "])"), collapse = " + "), ")/", Nstates)
    argsList <- alist(x =, parOne =, parTwo =, probs =)
    as.function(c(argsList, str2lang(strFun)))
  }
  densFun <- makeDensFun()
  calcDens <- function(parsVal, getParsD = FALSE){
    aux <- parsVal[, "est", ]
    auxDim <- c(nrow(value), 3, Nstates)
    parsD <- switch(as.character(mod$errorModel),
                    gaussian = parsVal,
                    gamma = array(rbind(aux^2/parsVal[, "sd", ]^2, aux/parsVal[, "sd", ]^2, parsVal[, "prob", ]), dim = auxDim),
                    beta = array(rbind(aux*(aux*(1-aux)/parsVal[, "sd", ]^2-1), (aux*(1-aux)/parsVal[, "sd", ]^2-1)*(1-aux), parsVal[, "prob", ]), dim = auxDim),
                    negbinomial = array(rbind(aux * aux / (parsVal[, "sd", ]^2 - aux), aux/parsVal[, "sd", ]^2, parsVal[, "prob", ]), dim = auxDim))
    if (getParsD) parsD else apply(parsD, 1, function(y) densFun(xx, y[1, ], y[2, ], y[3, ]))
  }
  parsValList <- lapply(parsTab, apply, 2, calcParsVal, value, mod, simplify = FALSE)
  densOut <- lapply(parsValList, lapply, calcDens, getParsD = getParsD)
  if (getParsD) densOut <- list(densOut, densFun = densFun)
  plotSlice <- function(sliceDens){
    lines(xx[1:samples], sliceDens[1:samples])
    if (addEcos) points(tail(xx, -samples), tail(sliceDens, -samples), pch = 16)
  }
  if (doPlot) {
    if (nrow(value) > 1) warning("Only curve(s) for the first value is/are plotted.")
    yrange <- c(0, 1.05 * max(unlist(densOut)))
    plot(range(resp), yrange, type = "n", ylab = "Probability density",
      xlab = xlab, ylim = yrange, axes = FALSE, yaxs = "i")
    if (xaxis) axis(1)
    axis(2, las = 2)
    box(bty = "l")
    lapply(densOut, lapply, plotSlice)
    if (plotEst) {
      plotEstFn <- function(parsVal){
        rgbVec <- col2rgb(setCol)
        cols <- rgb(rgbVec[1, ], rgbVec[2, ], rgbVec[3, ],
          alpha = 40 + parsVal[1, "prob", ] * 215, maxColorValue = 255)
        abline(v = parsVal[1, "est", ], lty = 2, lwd = 3, col = cols)
      }
      lapply(parsValList, lapply, plotEstFn)
    }
  }
  invisible(c(densOut, list(resp = xx)))
}

### 4.4. ==== Stability landscape from a mixture model ====
#' @title Plot stability landscape from a mixture model
#'
#' @description This function plots stability landscape for given predictor
#'
#' @param mod an object of class "mixglm"
#' @param form formula with one predictor specifying which variables to plot
#' @param threshold numerical value denoting minimum relative
#'   importance of visualized stable states and tipping points
#' @param addPoints logical value indicating if observations should be visualized
#' @param addMinMax logical value indicating if stable states and tipping points
#' should be visualized
#' @param randomSample integer specifying how many random samples from posterior
#' distribution to take instead of mean. Use \code{"NULL"} for mean. Plots instead
#' standard deviation of probability density of these samples.
#' @param otherPreds named vector of values of predictors not specified by form.
#' Default are zeros
#' @param eqiCol vector of colors of length 2 specifying colors for stable states
#' and tipping points respectively
#' @param ... parameters passed to \link[graphics]{image}
#'
#' @return Returns invisibly a list with scaled probability density matrix
#' (for each \code{randomSample}).
#'
#' @author Adam Klimes
#' @examples \dontrun{
#' set.seed(10)
#' n <- 200
#' x <- rnorm(n)
#' group <- rbinom(n, 1, 0.5)
#' y <- rnorm(n, 1 + 0.5 * x * c(-1, 1)[group + 1], 0.1)
#' plot(y ~ x)
#' dat <- data.frame(x, y)
#'
#' mod <- mixglm(
#'   stateValModels = y ~ x,
#'   stateProbModels = ~ x,
#'   statePrecModels = ~ x,
#'   inputData = dat,
#'   numStates = 2)
#' landscapeMixglm(mod)
#'
#' # uncertainty and stable states and tipping points for random samples
#' landscapeMixglm(mod, randomSample = 10)}
#' @export
#'
landscapeMixglm <- function(mod, form = NULL, threshold = 0, addPoints = TRUE,
  addMinMax = TRUE, randomSample = NULL, otherPreds = NULL, eqiCol = c("blue", "red"), ...){
  # input check
  if (!inherits(mod, "mixglm")) stop("'mod' must be an object of class 'mixglm'")
  form <- if (is.null(form))
    formula(paste(names(mod$data), "~", names(mod$constants)[4])) else formula(form)
  svar <- labels(terms(form))
  if (length(svar) > 1) stop("Specify only one predictor in 'form'")
  if (!svar %in% names(mod$constants))
    stop("Predictor specified in 'form' does not match 'mod'")
  if (!is.logical(addPoints)) stop("'addPoints' has to be logical")
  if (!is.logical(addMinMax)) stop("'addMinMax' has to be logical")
  # _
  resp <- mod$data[[1]]
  pred <- mod$constants[[svar]]
  grad <- seq(min(pred), max(pred), length.out = 500)
  if (!is.null(otherPreds)){
    valueDF <- data.frame(grad, matrix(otherPreds, 1, length(otherPreds),
                                       dimnames = list(NULL, names(otherPreds))))
    colnames(valueDF)[1] <- svar
  }  else valueDF <- grad
  slices <- sliceMixglm(mod, form, value = valueDF, byChains = FALSE,
    doPlot = FALSE, randomSample = randomSample)
  tipStableAll <- getMinMax(slices, threshold)
  tipStable <- tipStableAll$tipStable
  mats <- tipStableAll$matsSt
  is.one <- function(x) {
    out <- !is.null(x)
    if (out) out <- x == 1
    out
  }
  matPlot <- if (!(is.null(randomSample) | is.one(randomSample)))
    apply(array(do.call(c, mats), dim = c(dim(mats[[1]]), length(mats))),
    2, apply, 1, sd) else mats[[1]]
  arg <- list(...)
  if (is.null(arg$xlab)) arg$xlab <- svar
  if (is.null(arg$ylab)) arg$ylab <- names(mod$data)
  do.call(image, c(list(grad), list(slices$resp), list(t(matPlot)), arg))
  box()
  plotMinMax <- function(tipStable, xCoors, state, col, cex = 0.5){
    selStates <- tipStable[tipStable$state == state & tipStable$catSt == 1, ]
    points(rep(xCoors, nrow(selStates)), selStates$resp, pch = 16, cex = cex, col = col)
  }
  if (addPoints) points(pred, resp, cex = 0.4, pch = 16)
  if (addMinMax) {
    cex <- if (is.null(randomSample)) 0.5 else 0.1
    lapply(tipStable, function(x) Map(plotMinMax, x, grad, 0, col = eqiCol[2], cex = cex))
    lapply(tipStable, function(x) Map(plotMinMax, x, grad, 1, col = eqiCol[1], cex = cex))
  }
  invisible(mats)
}

### 4.5. ==== Randomly generate data from a mixture model ====
#' @title Randomly generate data from a mixture model
#'
#' @description Random generation of response variable from identified mixture model
#'
#' @param mod an object of class "mixglm"
#' @param n number of generated values per predictor value(s)
#' @param newdata data.frame of predictor values to be predicted for.
#'   If not provided, prediction is done for modelled data.
#'
#' @return A matrix with a column for each \code{n}
#'
#' @author Adam Klimes
#' @examples \dontrun{
#' set.seed(10)
#' n <- 200
#' x <- rnorm(n)
#' group <- rbinom(n, 1, 0.5)
#' y <- rnorm(n, 1 + 0.5 * x * c(-1, 1)[group + 1], 0.1)
#' plot(y ~ x)
#' dat <- data.frame(x, y)
#'
#' mod <- mixglm(
#'   stateValModels = y ~ x,
#'   stateProbModels = ~ x,
#'   statePrecModels = ~ x,
#'   inputData = dat,
#'   numStates = 2)
#' rMixglm(mod, n = 3)
#'
#' # for newdata
#' rMixglm(mod, n = 3, newdata = data.frame(x = c(0, 0.5)))}
#' @export
#'
rMixglm <- function(mod, n = 1, newdata = NULL){
  # input check
  if (!is.numeric(n)) stop("'n' has to be numeric")
  if (!(is.null(newdata) | is.data.frame(newdata)))
    stop("'newdata' has to be NULL or data.frame")
  if (names(mod$data) %in% colnames(newdata))
    newdata <- newdata[, -which(colnames(newdata) == names(mod$data)), drop = FALSE]
  if (is.null(newdata)) newdata <- as.data.frame(mod$constants[-(1:3)])
  if (!all(colnames(newdata) %in% names(mod$constants)))
    stop("some colnames of 'newdata' have not been found in 'mod'")
  #_
  form <- formula(paste(names(mod$data), "~", colnames(newdata)[1]))
  slices <- sliceMixglm(mod, form, value = newdata, byChains = FALSE,
    doPlot = FALSE, getParsD = TRUE)
  Nstates <- mod$constants$numStates
  sampleDist <- function(pars){
    comp <- sample(1:Nstates, prob = pars[3, ], size = n, replace = TRUE)
    rnorm(n, pars[1, comp], pars[2, comp])
  }
  out <- t(as.data.frame(apply(slices[[1]][[1]][[1]], 1, sampleDist, simplify = FALSE)))
  rownames(out) <- paste0("obs", 1:nrow(out))
  out
}
