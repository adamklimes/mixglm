% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixglm.R
\name{landscapeMixglm}
\alias{landscapeMixglm}
\title{Plot stability landscape from a mixture model}
\usage{
landscapeMixglm(
  mod,
  form = NULL,
  threshold = 0,
  addPoints = TRUE,
  addMinMax = TRUE,
  randomSample = NULL,
  otherPreds = NULL,
  eqiCol = c("blue", "red"),
  ...
)
}
\arguments{
\item{mod}{an object of class "mixglm"}

\item{form}{formula with one predictor specifying which variables to plot}

\item{threshold}{numerical value denoting minimum relative
importance of visualized stable states and tipping points}

\item{addPoints}{logical value indicating if observations should be visualized}

\item{addMinMax}{logical value indicating if stable states and tipping points
should be visualized}

\item{randomSample}{integer specifying how many random samples from posterior
distribution to take instead of mean. Use \code{"NULL"} for mean. Plots instead
standard deviation of probability density of these samples.}

\item{otherPreds}{named vector of values of predictors not specified by form.
Default are zeros}

\item{eqiCol}{vector of colors of length 2 specifying colors for stable states
and tipping points respectively}

\item{...}{parameters passed to \link[graphics]{image}}
}
\value{
Returns invisibly a list with scaled probability density matrix
(for each \code{randomSample}).
}
\description{
This function plots stability landscape for given predictor
}
\examples{
\dontrun{
set.seed(10)
n <- 200
x <- rnorm(n)
group <- rbinom(n, 1, 0.5)
y <- rnorm(n, 1 + 0.5 * x * c(-1, 1)[group + 1], 0.1)
plot(y ~ x)
dat <- data.frame(x, y)

mod <- mixglm(
  stateValModels = y ~ x,
  stateProbModels = ~ x,
  statePrecModels = ~ x,
  inputData = dat,
  numStates = 2)
landscapeMixglm(mod)

# uncertainty and stable states and tipping points for random samples
landscapeMixglm(mod, randomSample = 10)}
}
\author{
Adam Klimes
}
