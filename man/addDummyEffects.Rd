% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/3b_wrapperFunctions.R
\name{addDummyEffects}
\alias{addDummyEffects}
\alias{addFixedEffects}
\title{addDummyEffects}
\usage{
addDummyEffects(effectsObject, state)

addFixedEffects(effectsObject, state)
}
\arguments{
\item{effectsObject}{The monan Effects object to which the dummy
location effects should be added.}

\item{state}{The monan state (data) object that the effects object applies to}
}
\value{
An object of type effectsList.monan
}
\description{
A function to add one \link{alter_covariate} effect for each location
dummy variable included in the state. Dummy variables in the state object
are named "dummyXX", where XX can take any value.
This is particularly useful when including fixed effects for each location
to explicitly model the indegree of each location.
}
\examples{
myEffects <- createEffects(myState) |>
  addEffect(loops) |>
  addDummyEffects(myState)
}
