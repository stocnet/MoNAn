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
#' @param effects An object of class "effectsList.monan" for which the statistics of a multinomial
#' model should be calculated.
#' @param cache Outdated parameter, no need to specify.
#'
#' @return A data frame with N * M rows (N = mobile individuals, M = number of locations)
#' that specifies for each observation the statistics associated with moving to this location.
#' @export
#'
#' @seealso [createProcessState()], [createEffectsObject()]
#' 
#' @examples
#' myStatisticsFrame <- getMultinomialStatistics(myState, myEffects)
getMultinomialStatistics <-
  function(state, effects, cache = NULL) {
    
    if(!is.null(cache)){
      warning(paste("The cache object is automatically included in the process state", 
                    "since MoNAn version 1.0.0 and does not need to be specified anymore."))
    }
    if(is.null(state$cache)){
      stop(paste("The cache object is automatically included in the process state", 
                    "since MoNAn version 1.0.0. Please recreate your process state."))
    }
    
    dep.var <- state$dep.var
    cache <- state$cache
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
#' @param simulations outdated parameter, no need to specify
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
#' \donttest{
#' # goodness of fit
#' myGofIndegree <- gofMobilityNetwork(ans = myResDN, 
#'                                     gofFunction = getIndegree, 
#'                                     lvls = 1:100)
#' 
#' myGofTieWeight <- gofMobilityNetwork(ans = myResDN, 
#'                                      gofFunction = getTieWeights, 
#'                                      lvls = 1:30)
#' }
gofMobilityNetwork <-
  function(ans,
           gofFunction,
           lvls = NULL,
           simulations = NULL) {
    dep.var <- ans$state$dep.var
    
    if(!is.null(simulations)){
      warning("parameter simulations does not need to be specified anymore 
              since MoNAn version 1.0.0")
    }
    if (is.null(ans$deps)) {
      stop("ans object does not have simulated states stored; 
           use returnDeps = TRUE in estimation")
    }
    
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
      allCaches[[i]] <- createInternalCache(allStates[[i]], resourceCovariates = resCovsInCache)
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


#' monanGOF
#'
#' @rdname gofMobilityNetwork
monanGOF <- gofMobilityNetwork


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
#' \donttest{
#' plot(myGofIndegree,  lvls = 20:70)
#' plot(myGofTieWeight, lvls = 1:15)
#' }
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


#' print.effectsList.monan
#'
#' @param x An object of class "effectsList.monan".
#' @param ... For internal use only.
#'
#' @return The function `print.effectsList.monan` gives an overview of the 
#' specified effects.
#' @export
#'
#' @examples
#' myEffects
print.effectsList.monan <- function(x, ...) {
  
  df <- as.data.frame(matrix(nrow = 0, ncol = 4))
  colnames(df) <- c("name", paste0("covariate_", x$nodeset1), paste0("covariate_", x$nodeset2), "parameter")
  
  # fill table
  for (i in 1:length(x[["name"]])) {
    
    df[nrow(df) + 1, ]$name <- strsplit(x[["name"]][i], split = " ")[[1]][1]
    
    if (!is.null(formals(x[["effectFormulas"]][[i]])[["attribute.index"]])) {
      df[nrow(df), paste0("covariate_", x$nodeset1)] <- formals(x[["effectFormulas"]][[i]])[["attribute.index"]]
    } else {
      df[nrow(df), paste0("covariate_", x$nodeset1)] <- "-"
    }
    
    if (!is.null(formals(x[["effectFormulas"]][[i]])[["resource.attribute.index"]])) {
      df[nrow(df), paste0("covariate_", x$nodeset2)] <- formals(x[["effectFormulas"]][[i]])[["resource.attribute.index"]]
    } else {
      df[nrow(df), paste0("covariate_", x$nodeset2)] <- "-"
    }
    
    if (!is.null(formals(x[["effectFormulas"]][[i]])[["alpha"]])) {
      df[nrow(df), ]$parameter <- formals(x[["effectFormulas"]][[i]])[["alpha"]]
    } 
    if (!is.null(formals(x[["effectFormulas"]][[i]])[["lambda"]])) {
      df[nrow(df), ]$parameter <- formals(x[["effectFormulas"]][[i]])[["lambda"]]
    }
    if (is.na(df[nrow(df), ]$parameter)) {
      df[nrow(df), ]$parameter <- "-"
    }
    
  }
  rownames(df) <- c(1:nrow(df))
  
  colnames(df) <- c("effect name", paste0("cov. ", x$nodeset1), paste0("cov. ", x$nodeset2), "parameter")
  
  names(df)[1] <- format(names(df)[1], 
                          width = max(nchar(names(df[1])), max(nchar(df[,1]))),
                          justify = "left")
  names(df)[2] <- format(names(df)[2], 
                          width = max(nchar(names(df[2])), max(nchar(df[,2]))),
                          justify = "centre")
  names(df)[3] <- format(names(df)[3], 
                          width = max(nchar(names(df[3])), max(nchar(df[,3]))),
                          justify = "centre")
  names(df)[4] <- format(names(df)[4], 
                         width = max(nchar(names(df[3])), max(nchar(df[,3]))),
                         justify = "centre")
  df[,2] <- format(df[,2], 
                  width = max(nchar(names(df[2])), max(nchar(df[,2]))),
                  justify = "centre")
  df[,3] <- format(df[,3], 
                  width = max(nchar(names(df[3])), max(nchar(df[,3]))),
                  justify = "centre")
  df[,4] <- format(df[,4], 
                  width = max(nchar(names(df[3])), max(nchar(df[,3]))),
                  justify = "centre")
  df[,1] <- format(df[,1], justify = "left")
  
  cat("Effects\n")
  print(df, row.names = FALSE)
}


#' print.processState.monan
#'
#' @param x An object of class "processState.monan".
#' @param ... For internal use only.
#'
#' @return The function `print.processState.monan` gives an overview of the information 
#' included in the state object.
#' @export
#'
#' @examples
#' myState
print.processState.monan <- function(x, ...) {
  
  dep.var <- x$dep.var
  nodesets <- x[[dep.var]]$nodeSet
  covars <- names(x)[!(names(x) %in% c(dep.var, nodesets, "dep.var", "cache"))]

  cat(paste("dependent variable:", dep.var, "with",
            length(x[[nodesets[3]]]$ids), nodesets[3], "mobile between",
            length(x[[nodesets[1]]]$ids), nodesets[1]), "\n")
  cat("\n")

  # covariates of nodeset 1
  df1 <- as.data.frame(matrix(nrow = 0, ncol = 3))
  colnames(df1) <- c("cov. name", "range", "mean")

  # fill table
  if(length(covars) > 0){
    for (i in 1:length(covars)) {
      
      if (nodesets[1] %in% x[[covars[i]]]$nodeSet) {
        
        df1[nrow(df1) + 1, ]$'cov. name' <- covars[i]
        df1[nrow(df1), ]$range <-
          paste0(round(min(x[[covars[i]]]$data, na.rm = TRUE), digits = 2), "-",
                 round(max(x[[covars[i]]]$data, na.rm = TRUE), digits = 2))
        df1[nrow(df1), ]$mean <- round(mean(x[[covars[i]]]$data, na.rm = TRUE),
                                       digits = 2)
      }
    }

    # if there are covariates for this nodeset, assign rownames
    if (nrow(df1) > 0) {
      rownames(df1) <- c(1:nrow(df1))
    }
  }
  
  cat(paste("covariates of", nodesets[1]), "\n")
  
  # if there are covariates for this nodeset, print table
  if (nrow(df1) > 0) {
    
    names(df1)[1] <- format(names(df1)[1], 
                            width = max(nchar(names(df1[1])), max(nchar(df1[,1]))),
                            justify = "left")
    names(df1)[2] <- format(names(df1)[2], 
                            width = max(nchar(names(df1[2])), max(nchar(df1[,2]))),
                            justify = "centre")
    names(df1)[3] <- format(names(df1)[3], 
                            width = max(nchar(names(df1[3])), max(nchar(df1[,3]))),
                            justify = "centre")
    df1[,2:3] <- format(df1[,2:3], justify = "centre")
    df1[,1] <- format(df1[,1], justify = "left")
    print(df1, row.names = FALSE)
    
  } else { # if there are no covariates for the resp. nodeset
    cat("no covariates available for nodeset", nodesets[1], "\n")
  }
  
  cat("\n")

  # covariates of nodeset 2
  df2 <- as.data.frame(matrix(nrow = 0, ncol = 3))
  colnames(df2) <- c("cov. name", "range", "mean")
  # fill table
  if(length(covars) > 0){
    for (i in 1:length(covars)) {
      
      if (nodesets[3] %in% x[[covars[i]]]$nodeSet) {
        
        df2[nrow(df2) + 1, ]$'cov. name' <- covars[i]
        df2[nrow(df2), ]$range <-
          paste0(round(min(x[[covars[i]]]$data, na.rm = TRUE), digits = 2), "-",
                 round(max(x[[covars[i]]]$data, na.rm = TRUE), digits = 2))
        df2[nrow(df2), ]$mean <- round(mean(x[[covars[i]]]$data, na.rm = TRUE),
                                       digits = 2)
      }
    }
    
    # if there are covariates for this nodeset, assign rownames
    if (nrow(df2) > 0) {
      rownames(df2) <- c(1:nrow(df2))
    }
  }

  cat(paste("covariates of", nodesets[3]), "\n")
  
  # if there are covariates for this nodeset, print table
  if (nrow(df2) > 0) {
    
    names(df2)[1] <- format(names(df2)[1], 
                            width = max(nchar(names(df2[1])), max(nchar(df2[,1]))),
                            justify = "left")
    names(df2)[2] <- format(names(df2)[2], 
                            width = max(nchar(names(df2[2])), max(nchar(df2[,2]))),
                            justify = "centre")
    names(df2)[3] <- format(names(df2)[3], 
                            width = max(nchar(names(df2[3])), max(nchar(df2[,3]))),
                            justify = "centre")
    
    df2[,2:3] <- format(df2[,2:3], justify = "centre")
    df2[,1] <- format(df2[,1], justify = "left")
    print(df2, row.names = FALSE)
    
  }  else { # if there are no covariates for the resp. nodeset
    cat("no covariates available for nodeset", nodesets[3], "\n")
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