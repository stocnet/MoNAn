########## auxiliaryFunctions


# autoCorrelationTestMoNAn
#' Title
#'
#' @param dep.var
#' @param ans
#'
#' @return
#' @export
#'
#' @examples
autoCorrelationTestMoNAn <- function(dep.var, ans) {
  # give error if no deps in ans obj
  if (is.null(ans$deps))
    stop("ans object does not have simulations stored; use returnDeps = T in estimation")

  # get number of simulated nets - 1
  nSims <- length(ans$deps) - 1

  # loop through all of them to see how many resources are re-allocated
  means <- c()
  for (i in 1:nSims) {
    means[i] <-
      mean(ans$deps[[i]]$state[[dep.var]]$data[, 2] == ans$deps[[(i + 1)]]$state[[dep.var]]$data[, 2])
  }

  mean(means)
}


# extractTracesMoNAn
#' Title
#'
#' @param dep.var
#' @param ans
#' @param effects
#'
#' @return
#' @export
#'
#' @examples
extractTracesMoNAn <- function(dep.var, ans, effects) {
  # give error if no deps in ans
  if (is.null(ans$deps))
    stop("ans object does not have simulated states stored; use returnDeps = T in estimation")

  # get statistics of the observed network using state and cache from ans object
  obsStats <-
    getNetworkStatistics(dep.var, ans$state, ans$cache, effects)

  # get a matrix of statistics from the simulated data in the ans object
  simStatsMatrix <-
    Reduce("rbind", lapply(ans$deps, function(x)
      getNetworkStatistics(dep.var, x$state, x$cache, effects)))

  results <- list(
    observedStats = obsStats,
    simulatedStats = simStatsMatrix,
    effectNames = effects$name
  )

  class(results) <- "traces.monan"

  results
}


# getInitialEstimates
#' Title
#'
#' @param state
#' @param cache
#' @param effects
#' @param dep.var
#' @param joint
#'
#' @return
#' @export
#'
#' @examples
getInitialEstimates <-
  function(state, cache, effects, dep.var, joint = T) {
    statisticsFrame <-
      getMultinomialStatistics(state, cache, effects, dep.var)

    if (joint) {
      initialParam <- getMultinomialEstimates(statisticsFrame)
    } else {
      initialParam <- getIndividualEstimates(statisticsFrame)
    }

    initialParam
}


# gofDistributionNetwork
#' Title
#'
#' @param ans
#' @param simulations
#' @param gofFunction
#' @param dep.var
#' @param lvls
#'
#' @return
#' @export
#'
#' @examples
gofDistributionNetwork <-
  function(ans,
           simulations,
           gofFunction,
           dep.var = NULL,
           lvls = NULL) {
    if (is.null(dep.var))
      dep.var <- names(ans$state)[1]

    # generate a list that contains all states for the dep.var with the observed in first place
    allStates <- list()
    allStates[[1]] <- ans$state
    for (i in 2:(length(simulations) + 1)) {
      allStates[[i]] <- simulations[[i - 1]]$state
    }
    allCaches <- list()
    allCaches[[1]] <- ans$cache
    for (i in 2:(length(simulations) + 1)) {
      allCaches[[i]] <- simulations[[i - 1]]$cache
    }

    gofStats <-
      lapply(1:length(allCaches), function(i)
        gofFunction(
          state = allStates[[i]],
          cache = allCaches[[i]],
          dep.var = dep.var,
          lvls = lvls
        ))
    simStats <- gofStats
    simStats[[1]] <- NULL
    gofRes <- list(observed = gofStats[[1]], simulated = simStats)
    class(gofRes) <- "gof.stats.monan"
    return(gofRes)
}


# plot.gof.stats.monan
#' Title
#'
#' @param gofObject
#' @param lvls
#'
#' @return
#' @export
#'
#' @examples
plot.gof.stats.monan <- function(gofObject, lvls = NULL) {
  if (is.null(lvls))
    lvls <- 1:length(gofObject$observed)
  simStats <- Reduce(rbind, gofObject$simulated)
  boxplot(simStats[, lvls])
  lines(gofObject$observed, col = "red")
}


# plot.traces.monan
#' Title
#'
#' @param xx
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.traces.monan <- function(xx, ...) {
  nParams <- length(xx[[1]])
  nSims <- length(xx[[2]][, 1])
  for (i in 1:nParams) {
    plot(xx[[2]][, i], main = xx[[3]][i])
    lines(x = 1:nSims,
          y = rep(xx[[1]][i], nSims),
          col = "red")
  }
}


# print.result.monan
#' Title
#'
#' @param x
#' @param covMat
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.result.monan <- function(x, covMat = F, ...) {
  reslt <-
    data.frame(
      Effects = names(x$estimates),
      Estimates =  as.numeric(x$estimates),
      StandardErrors = x$standardErrors,
      Convergence = as.numeric(x$convergenceStatistics)
    )

  cat("Results\n")

  print(reslt)
  cat("\n")

  if (covMat) {
    cat("Covariance Matrix\n")

    print(x$covarianceMatrix)
  }
}


# print.scoretest.monan
#' Title
#'
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
print.scoretest.monan <- function(x, ...) {
  reslt <-
    data.frame(
      Effects = x$effects$name,
      pValuesParametric = as.numeric(x$pValuesParametric),
      pValuesNonParametric = as.numeric(x$pValuesNonParametric)
    )

  cat("Results\n")

  print(reslt)
  cat("\n")
  cat(
    " Parametric p-values: small = more significant \n Non-parametric p-values: further away from 0.5 = more significant"
  )
}


# scoreTestMoNAn
#' Title
#'
#' @param dep.var
#' @param ans
#' @param effects
#'
#' @return
#' @export
#'
#' @examples
scoreTestMoNAn <- function(dep.var, ans, effects) {
  # give error if no deps in ans
  if (is.null(ans$deps))
    stop("ans object does not have simulated states stored; use returnDeps = T in estimation")

  # get statistics of the observed network using state and cache from ans object
  obsStats <-
    getNetworkStatistics(dep.var, ans$state, ans$cache, effects)

  # get a matrix of statistics from the simulated data in the ans object
  simStatsMatrix <-
    Reduce("rbind", lapply(ans$deps, function(x)
      getNetworkStatistics(dep.var, x$state, x$cache, effects)))

  # get parametric p-values using mean and sd
  p.vals.par <-
    2 * pnorm(-abs(((obsStats - colMeans(simStatsMatrix)) / apply(simStatsMatrix, 2, sd)
    )))

  # get non-parametric p-values using ecdf
  p.vals.np <-
    (apply(rbind(obsStats, simStatsMatrix), 2, function(x)
      ecdf(x[2:length(x)])(x[1])))

  # return results
  result.score <- list(
    effects = effects,
    observedStatstics = obsStats,
    simulatedStatistics = simStatsMatrix,
    pValuesParametric = p.vals.par,
    pValuesNonParametric = p.vals.np
  )

  class(result.score) <- "scoretest.monan"

  return(result.score)
}
