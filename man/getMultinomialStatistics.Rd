% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_auxiliaryFunctions.R
\name{getMultinomialStatistics}
\alias{getMultinomialStatistics}
\title{getMultinomialStatistics}
\usage{
getMultinomialStatistics(state, effects, cache = NULL)
}
\arguments{
\item{state}{An object of class "processState.monan" that stores all information to be used in the model.}

\item{effects}{An object of class "effectsList.monan" for which the statistics of a multinomial
model should be calculated.}

\item{cache}{Outdated parameter, no need to specify.}
}
\value{
A data frame with N * M rows (N = mobile individuals, M = number of locations)
that specifies for each observation the statistics associated with moving to this location.
}
\description{
One updating step in simulating the mobility network model can be expressed
as a multinomial logit model. Extracting the statistics for such a model allows
a straight-forward estimation of a multinomial logit model to get initial
estimates for the full mobility model, which increases the chances of model
convergence in the first run of the estimation considerably.
}
\examples{
myStatisticsFrame <- getMultinomialStatistics(myState, myEffects)
}
\seealso{
\code{\link[=createProcessState]{createProcessState()}}, \code{\link[=createEffectsObject]{createEffectsObject()}}
}
