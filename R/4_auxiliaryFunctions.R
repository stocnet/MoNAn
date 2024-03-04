########## auxiliaryFunctions


#' autoCorrelationTest
#' 
#' The autoCorrelationTest indicates the degree to which the values of the dependent 
#' variable of consecutive draws from the chain in phase 3 are correlated. Here lower 
#' values are better. Values above 0.5 are very problematic and indicate that a 
#' higher thinning is needed.
#' 
#' @param ans An object of class "result.monan" resulting from an estimation with the function [estimateMobilityNetwork()].
#'
#' @return A number indicating the auto-correlation.
#' @export
#'
#' @examples
#' # regression diagnostics
#' autoCorrelationTest(myResDN)
autoCorrelationTest <- function(ans) {
  dep.var <- ans$state$dep.var
  # give error if no deps in ans obj
  if (is.null(ans$deps)) {
    stop("ans object does not have simulations stored; use returnDeps = TRUE in estimation")
  }

  # get number of simulated nets - 1
  nSims <- length(ans$deps) - 1

  # loop through all of them to see how many resources are re-allocated
  means <- c()
  for (i in 1:nSims) {
    means[i] <-
      mean(ans$deps[[i]]$data[, 2] == ans$deps[[(i + 1)]]$data[, 2])
  }

  mean(means)
}


#' extractTraces
#' 
#' This function shows the values of simulated statistics in Phase 3 for subsequent draws from the chain.
#' Ideally, the plots show points randomly scattered around the red line, which indicates the statistics in the data.
#' 
#' @param ans An object of class "result.monan" resulting from an estimation with the function [estimateMobilityNetwork()].
#' @param effects An object of class "effectsList.monan" used in the estimation.
#'
#' @return The function `extractTraces` returns a list that includes 
#' (1) the observed statistics for all effects, 
#' (2) the distribution of statistics for all simulations and
#' (3) effect names.
#' It is recommended to use the plotting function to inspect the traces.
#' @export
#'
#' @seealso [createEffectsObject()]
#'
#' @examples
#' \donttest{
#' # regression diagnostics
#' traces <- extractTraces(myResDN, myEffects)
#' }
extractTraces <- function(ans, effects) {
  dep.var <- ans$state$dep.var
  # give error if no deps in ans
  if (is.null(ans$deps)) {
    stop("ans object does not have simulated states stored; use returnDeps = TRUE in estimation")
  }
  
  # get statistics of the observed network using state and cache from ans object
  obsStats <-
    getNetworkStatistics(dep.var, ans$state, ans$cache, effects)
  
  # get a matrix of statistics from the simulated data in the ans object
  simStatsMatrix <-
    Reduce("rbind", lapply(ans$deps, function(x) {
      getNetStatsFromDeps(dep.var, x, ans, effects)
    }))
  
  results <- list(
    observedStats = obsStats,
    simulatedStats = simStatsMatrix,
    effectNames = effects$name
  )
  
  class(results) <- "traces.monan"
  
  results
}


#' getMultinomialStatistics
#'
#' One updating step in simulating the mobility network model can be expressed
#' as a multinomial logit model. Extracting the statistics for such a model allows
#' a straight-forward estimation of a multinomial logit model to get initial 
#' estimates for the full mobility model, which increases the chances of model 
#' convergence in the first run of the estimation considerably.
#'
#' @param state An object of class "processState.monan" that stores all information to be used in the model.
#' @param cache A cache object that contains intermediate information used for estimation.
#' @param effects An object of class "effectsList.monan" for which the statistics of a multinomial
#' model should be calculated.
#'
#' @return A data frame with N * M rows (N = mobile individuals, M = number of locations)
#' that specifies for each observation the statistics associated with moving to this location.
#' @export
#'
#' @seealso [createProcessState()], [createWeightedCache()], [createEffectsObject()]
#' 
#' @examples
#' myStatisticsFrame <- getMultinomialStatistics(myState, myCache, myEffects)
getMultinomialStatistics <-
  function(state, cache, effects) {
    
    dep.var <- state$dep.var
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



#' gofMobilityNetwork
#' 
#' Akin to ERGMs, goodness of fit testing is available to see whether auxiliary 
#' statistics are well captured by the model. The logic behind gof testing for network models is outlined in 
#' Hunter et al. (2008) and Lospinoso and Snijders (2019).
#'
#' @param ans An object of class "result.monan" resulting from an estimation 
#' with the function [estimateMobilityNetwork()] using the option deps = TRUE.
#' @param gofFunction A gof function that specifies which auxiliary outcome should be used, 
#' e.g., "getIndegree" or "getTieWeights".
#' @param lvls The values for which the gofFunction should be calculated/plotted.
#'
#' @return The function `gofMobilityNetwork` returns a list containing 
#' (1) the observed values of the auxiliary statistics and
#' (2) a list of the simulated values of the auxiliary statistics.
#' @export
#'
#' @seealso [getIndegree()], [getTieWeights()]
#'
#' @references Hunter, D. R., Goodreau, S. M., & Handcock, M. S. (2008). Goodness of fit of social network models. 
#' \emph{Journal of the american statistical association}, 103(481), 248-258.
#' 
#' Lospinoso, J., & Snijders, T. A. (2019). 
#' Goodness of fit for stochastic actor-oriented models. \emph{Methodological Innovations}, 12(3).
#' 
#' 
#' @examples
#' # goodness of fit
#' myGofIndegree <- gofMobilityNetwork(ans = myResDN, 
#'                                         gofFunction = getIndegree, 
#'                                         lvls = 1:100)
#' 
#' myGofTieWeight <- gofMobilityNetwork(ans = myResDN, 
#'                                          gofFunction = getTieWeights, 
#'                                          lvls = 1:30)
gofMobilityNetwork <-
  function(ans,
           gofFunction,
           lvls = NULL) {
    dep.var <- ans$state$dep.var
    
    # generate a list that contains all states for the dep.var with the observed in first place
    allStates <- list()
    allStates[[1]] <- ans$state
    for (i in 2:(length(ans$deps) + 1)) {
      allStates[[i]] <- allStates[[1]]
      allStates[[i]][[dep.var]] <- ans$deps[[(i-1)]]
    }
    allCaches <- list()
    allCaches[[1]] <- ans$cache
    resCovsInCache <- names(ans$cache[[dep.var]]$resourceNetworks)
    for (i in 2:(length(ans$deps) + 1)) {
      allCaches[[i]] <- createWeightedCache(allStates[[i]], resourceCovariates = resCovsInCache)
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

#' gofDistributionNetwork
#'
#' @rdname gofMobilityNetwork
gofDistributionNetwork <- gofMobilityNetwork

#' plot.gof.stats.monan
#'
#' @rdname gofMobilityNetwork
#' @param x An object of class "gof.stats.monan".
#' @param lvls The values for which the gofFunction should be calculated/plotted.
#' @param ... Additional plotting parameters, use discouraged.
#'
#' @return The function `plot.gof.stats.monan` returns violin plots of the 
#' gof tests with observed values superimposed in red.
#' @export
#'
#' @examples
#' plot(myGofIndegree,  lvls = 20:70)
#' plot(myGofTieWeight, lvls = 1:15)
plot.gof.stats.monan <- function(x, lvls, ...) {
  if (is.null(lvls)) {
    lvls <- 1:length(x$observed)
  }
  simStats <- Reduce(rbind, x$simulated)
  boxplot(simStats[, lvls])
  lines(x$observed[lvls], col = "red")
}


#' plot.traces.monan
#'
#' @rdname extractTraces
#' @param x An object of class "traces.monan".
#' @param ... Additional plotting parameters, use not recommended.
#'
#' @return The function `plot.traces.monan` shows a scatter plot of the
#' statistics of simulated networks from phase three of the esimtation.
#' @export
#'
#' @examples
#' \donttest{
#' plot(traces)
#' }
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


#' print.result.monan
#'
#' @rdname estimateMobilityNetwork
#' @param x An object of class "result.monan".
#' @param covMat Logical: indicating whether the covariance matrix should be printed.
#' @param ... For internal use only.
#'
#' @return The function `print.result.monan` prints the results from a monan
#' estimation with three columns indicating the estimate, the standard error,
#' and the convergence statistic.
#' @export
#' 
#' @examples
#' myResDN
print.result.monan <- function(x, covMat = FALSE, ...) {
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



#' scoreTest
#'
#' Based on an estimated model, a score-type test is available that shows whether 
#' statistics representing non-included effects are well represented. If this is 
#' not the case, it is likely that including them will result in significant estimates.
#'
#' @param ans An object of class "result.monan" resulting from an estimation with the function [estimateMobilityNetwork()].
#' @param effects An object of class "effectsList.monan" in which the non included effects that should
#' be tested are specified.
#'
#' @return The function `scoreTest` returns basic values to calculate
#' parametric and non-parametric p-values
#' for each tested effect.
#' @export
#'
#' @seealso [createEffectsObject()]
#'
#' @examples
#' \donttest{
#' # test whether other effects should be included
#' myEffects2 <- createEffects(myState) |>
#'   addEffect(transitivity_min)
#' 
#' test_ME.2 <- scoreTest(myResDN, myEffects2)
#' }
scoreTest <- function(ans, effects) {
  dep.var <- ans$state$dep.var
  # give error if no deps in ans
  if (is.null(ans$deps)) {
    stop("ans object does not have simulated states stored; use returnDeps = TRUE in estimation")
  }

  # get statistics of the observed network using state and cache from ans object
  obsStats <-
    getNetworkStatistics(dep.var, ans$state, ans$cache, effects)

  # get a matrix of statistics from the simulated data in the ans object
  simStatsMatrix <-
    Reduce("rbind", lapply(ans$deps, function(x) {
      getNetStatsFromDeps(dep.var, x, ans, effects)
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


#' print.scoretest.monan
#'
#' @rdname scoreTest
#' @param x An object of class "scoretest.monan".
#' @param ... For internal use only.
#'
#' @export
#'
#' @return The function `print.scoretest.monan` shows parametric and non-parametric p-values
#' for each tested effect.
#' @examples
#' \donttest{
#' test_ME.2
#' }
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