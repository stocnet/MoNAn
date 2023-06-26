########## auxiliaryFunctions


# autoCorrelationTest
#' Title
#'
#' @param dep.var
#' @param ans
#'
#' @return
#' @export
#'
#' @seealso [estimateMobilityNetwork()]
#'
#' @examples
#' # regression diagnostics
#' autoCorrelationTest(myDependentVariable, myResDN)
autoCorrelationTest <- function(dep.var, ans) {
  # give error if no deps in ans obj
  if (is.null(ans$deps)) {
    stop("ans object does not have simulations stored; use returnDeps = T in estimation")
  }

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


# extractTraces
#' Title
#'
#' @param dep.var
#' @param ans
#' @param effects
#'
#' @return
#' @export
#'
#' @seealso [estimateMobilityNetwork()], [createEffectsObject()]
#'
#' @examples
#' # regression diagnostics
#' traces <- extractTraces(myDependentVariable, myResDN, myEffects)
#' plot(traces)
extractTraces <- function(dep.var, ans, effects) {
  # give error if no deps in ans
  if (is.null(ans$deps)) {
    stop("ans object does not have simulated states stored; use returnDeps = T in estimation")
  }

  # get statistics of the observed network using state and cache from ans object
  obsStats <-
    getNetworkStatistics(dep.var, ans$state, ans$cache, effects)

  # get a matrix of statistics from the simulated data in the ans object
  simStatsMatrix <-
    Reduce("rbind", lapply(ans$deps, function(x) {
      getNetworkStatistics(dep.var, x$state, x$cache, effects)
    }))

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


# getMultinomialStatistics
#' Title
#'
#' @param state 
#' @param cache 
#' @param effects 
#' @param dep.var 
#'
#' @return
#' @export
#'
#' @examples
#' myStatisticsFrame <- getMultinomialStatistics(myState, myCache, myEffects, myDependentVariable)
getMultinomialStatistics <-
  function(state, cache, effects, dep.var) {
    # create a list that will store all statistics
    statsVec <- list()
    
    # get number of options where the resource can go
    actors1 <- state[[dep.var]]$nodeSet[1]
    nActors1 <- length(state[[actors1]]$ids)
    
    # loop through all resources
    for (number in 1:nrow(state[[dep.var]]$data)) {
      transferUnderAnalysis <- state[[dep.var]]$data[number, ]
      
      # get a new cache without the transfer under analysis
      
      reducedCache <- cache
      reducedCache[[dep.var]] <-
        updateWeightedCache(
          cache = reducedCache[[dep.var]],
          state = state,
          resourceID = number,
          dep.var = dep.var,
          sender = transferUnderAnalysis[1],
          receiver = transferUnderAnalysis[2],
          update = -1
        )
      
      
      # calculate the statistics for the every option
      
      stats <- matrix(NA, nActors1, length(effects$effectFormulas) + 3)
      
      for (target in 1:nActors1) {
        newCache <- reducedCache
        # newCache[[dep.var]] <- updateWeightedCache(cache = newCache[[dep.var]],
        #                                            sender = transferUnderAnalysis[1],
        #                                            receiver = target,
        #                                            update = 1)
        
        # create a state that would be with the transfer in that place
        newState <- state
        newState[[dep.var]]$data[number, ][2] <- target
        
        stats[target, 4:(length(effects$effectFormulas) + 3)] <-
          unlist(lapply(effects$effectFormulas, function(f) {
            f(
              dep.var = dep.var,
              state = newState,
              cache = newCache,
              i = transferUnderAnalysis[1],
              j = target,
              edge = number,
              update = 1
            )
          }))
      }
      
      stats[, 1] <- number
      stats[, 2] <- 1:nActors1
      stats[, 3] <- 1:nActors1 == transferUnderAnalysis[2]
      
      statsVec[[number]] <- stats
    }
    statsMat <- (Reduce(rbind, statsVec))
    
    colnames(statsMat) <-
      c("resource", "target", "choice", effects$name)
    
    return(as.data.frame(statsMat))
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
#' @seealso [estimateMobilityNetwork()], [getIndegree()],
#' [getTieWeights()]
#'
#' @examples
#' # goodness of fit
#' myGofIndegree <- gofDistributionNetwork(ans = myResDN, 
#'                                         simulations = myResDN$deps, 
#'                                         gofFunction = getIndegree, 
#'                                         lvls = 1:100)
#' 
#' myGofTieWeight <- gofDistributionNetwork(ans = myResDN, 
#'                                          simulations = myResDN$deps, 
#'                                          gofFunction = getTieWeights, 
#'                                          lvls = 1:30)
gofDistributionNetwork <-
  function(ans,
           simulations,
           gofFunction,
           dep.var = NULL,
           lvls = NULL) {
    if (is.null(dep.var)) {
      dep.var <- names(ans$state)[1]
    }

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
      lapply(1:length(allCaches), function(i) {
        gofFunction(
          state = allStates[[i]],
          cache = allCaches[[i]],
          dep.var = dep.var,
          lvls = lvls
        )
      })
    simStats <- gofStats
    simStats[[1]] <- NULL
    gofRes <- list(observed = gofStats[[1]], simulated = simStats)
    class(gofRes) <- "gof.stats.monan"
    return(gofRes)
  }


# plot.gof.stats.monan
#' Title
#'
#' @rdname gofDistributionNetwork
#' @param x a gofObject
#'
#' @return
#' @export
#'
#' @examples
plot.gof.stats.monan <- function(x, ...) {
  if (is.null(lvls)) {
    lvls <- 1:length(x$observed)
  }
  simStats <- Reduce(rbind, x$simulated)
  boxplot(simStats[, lvls])
  lines(x$observed, col = "red")
}


# plot.traces.monan
#' Title
#'
#' @rdname extractTraces
#' @param x
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.traces.monan <- function(x, ...) {
  nParams <- length(x[[1]])
  nSims <- length(x[[2]][, 1])
  for (i in 1:nParams) {
    plot(x[[2]][, i], main = x[[3]][i])
    lines(
      x = 1:nSims,
      y = rep(x[[1]][i], nSims),
      col = "red"
    )
  }
}


# print.result.monan
#' Title
#'
#' @rdname estimateMobilityNetwork
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
      Estimates = as.numeric(x$estimates),
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


# scoreTest
#' Title
#'
#' @param dep.var
#' @param ans
#' @param effects
#'
#' @return
#' @export
#'
#' @seealso [createEffectsObject()], [estimateMobilityNetwork()]
#'
#' @examples
#' # test whether other effects should be included
#' myEffects2 <- createEffectsObject(
#'   list(
#'     list("loops"),
#'     list("min_reciprocity"),
#'     list("dyadic_covariate", attribute.index = "sameRegion"),
#'     list("alter_covariate", attribute.index = "size"),
#'     list("resource_covar_to_node_covar", attribute.index = "region", 
#'           resource.attribute.index = "sex"),
#'     list("loops_resource_covar", resource.attribute.index = "sex"),
#'     list("min_transitivity")
#'   )
#' )
#' 
#' test_ME.2 <- scoreTest(myDependentVariable, myResDN, myEffects2)
#' test_ME.2
scoreTest <- function(dep.var, ans, effects) {
  # give error if no deps in ans
  if (is.null(ans$deps)) {
    stop("ans object does not have simulated states stored; use returnDeps = T in estimation")
  }

  # get statistics of the observed network using state and cache from ans object
  obsStats <-
    getNetworkStatistics(dep.var, ans$state, ans$cache, effects)

  # get a matrix of statistics from the simulated data in the ans object
  simStatsMatrix <-
    Reduce("rbind", lapply(ans$deps, function(x) {
      getNetworkStatistics(dep.var, x$state, x$cache, effects)
    }))

  # get parametric p-values using mean and sd
  p.vals.par <-
    2 * pnorm(-abs(((obsStats - colMeans(simStatsMatrix)) / apply(simStatsMatrix, 2, sd)
    )))

  # get non-parametric p-values using ecdf
  p.vals.np <-
    (apply(rbind(obsStats, simStatsMatrix), 2, function(x) {
      ecdf(x[2:length(x)])(x[1])
    }))

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


# print.scoretest.monan
#' Title
#'
#' @rdname scoreTest
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