% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2c_effectFunctions_concentration.R
\name{concentration_prop}
\alias{concentration_prop}
\title{concentration_prop}
\usage{
concentration_prop(
  dep.var = 1,
  state,
  cache,
  i,
  j,
  edge,
  update,
  getTargetContribution = FALSE
)
}
\arguments{
\item{dep.var}{}

\item{state}{}

\item{cache}{}

\item{i}{}

\item{j}{}

\item{edge}{}

\item{update}{}

\item{getTargetContribution}{}
}
\value{
Returns the change statistic or target statistic of the effect for
internal use by the estimation algorithm.
}
\description{
Is there a bandwagon effect in mobility, i.e. do mobile individuals move to locations
that are the destination of many others from their origin? The functional form of this
statistic assumes that individuals consider the proportions of individuals (coming from
the same origin) going to a certain destination, instead of the total number.
}
\keyword{internal}
