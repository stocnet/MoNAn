########## hiddenFunctions

# binarizeNetwork
binarizeNetwork <- function(network) {
  net <- 1 * (network >= 1)
  diag(net) <- 0
  net
}


# getCovarianceMatrix
getCovarianceMatrix <- function(statistics) {
  meanStatistics <- colMeans(statistics)
  meanStatsMatrix <- meanStatistics %*% t(meanStatistics)

  nObs <- dim(statistics)[1]
  observationMatrices <-
    lapply(1:nObs, function(i) {
      statistics[i, ] %*% t(statistics[i, ])
    })
  meanObservationMatrix <- Reduce("+", observationMatrices) / nObs

  covMatrix <- meanObservationMatrix - meanStatsMatrix
  return(covMatrix)
}


# getNetworkStatistics
getNetworkStatistics <- function(dep.var, state, cache, effects) {
  actors1 <- state[[dep.var]]$nodeSet[1]
  actors2 <- state[[dep.var]]$nodeSet[2]
  nActors1 <- length(state[[actors1]]$ids)
  nActors2 <- length(state[[actors2]]$ids)

  targetStatistics <- unlist(lapply(
    effects$effectFormulas,
    function(f) {
      sum(apply(
        expand.grid(1:nActors1, 1:nActors2),
        1,
        function(v) {
          f(
            dep.var = dep.var,
            state = state,
            cache = cache,
            i = v[[1]],
            j = v[[2]],
            update = NULL,
            edge = NULL,
            getTargetContribution = T
          )
        }
      ))
    }
  ))

  names(targetStatistics) <- effects$name

  return(targetStatistics)
}


# runPhase1
runPhase1 <- function(dep.var,
                      state,
                      cache,
                      effects,
                      initialParameters,
                      burnInN1,
                      iterationsN1,
                      thinningN1,
                      gainN1,
                      multinomialProposal = F,
                      allowLoops,
                      verbose = F) {
  # simulate statistic matrix
  statisticsMatrix <-
    simulateStatisticVectors(
      dep.var,
      state,
      cache,
      effects,
      initialParameters,
      burnIn = burnInN1,
      iterations = iterationsN1,
      thinning = thinningN1,
      multinomialProposal = multinomialProposal,
      allowLoops,
      verbose = verbose
    )

  # calculate covariance matric
  covarianceMatrix <- getCovarianceMatrix(statisticsMatrix)

  # calculate sensitivity
  sensitivityVector <- 1 / diag(covarianceMatrix)

  # initial update step if initial parameters are zero
  if (all(initialParameters == 0)) {
    observedStatistics <-
      getNetworkStatistics(dep.var, state, cache, effects)
    averageStatistics <- colMeans(statisticsMatrix)
    updatedParameters <-
      initialParameters - gainN1 * sensitivityVector * (averageStatistics - observedStatistics)
  } else {
    updatedParameters <- initialParameters
  }

  # returns initial estimates and sensitivity vector
  return(list(
    estimates = updatedParameters,
    sensitivityVector = sensitivityVector
  ))
}


# runPhase2
runPhase2 <- function(dep.var,
                      state,
                      cache,
                      effects,
                      initialParameters,
                      sensitivityVector,
                      burnInN2,
                      nsubN2,
                      initGain,
                      thinningN2,
                      initialIterationsN2,
                      parallel,
                      cpus,
                      multinomialProposal = F,
                      allowLoops,
                      verbose = F) {
  # calculate observed statistics
  observedStatistics <-
    getNetworkStatistics(dep.var, state, cache, effects)

  # initialize parallel computing
  if (parallel && cpus > 1) {
    sfInit(parallel = T, cpus = cpus)
    # TODO. Replace this long command with sfLibrary("NetDist") once the package is packaged
    sfLibrary("MoNAn", character.only=TRUE)
  } else {
    parallel <- F
  }


  # TODO PARALLEL: n burn-ins with n states
  # burn in
  if (parallel) {
    res <-
      sfLapply(
        1:cpus,
        simulateNSteps,
        dep.var = dep.var,
        state = state,
        cache = cache,
        effects = effects,
        parameters = initialParameters,
        allowLoops = allowLoops,
        multinomialProposal = multinomialProposal,
        n = burnInN2
      )
  } else {
    res <-
      simulateNSteps(
        dep.var,
        state,
        cache,
        effects,
        initialParameters,
        allowLoops,
        multinomialProposal = multinomialProposal,
        n = burnInN2
      )
    state <- res$state
    cache <- res$cache
  }

  # run n sub phases
  gain <- initGain
  iterations <- initialIterationsN2
  parameters <- initialParameters
  for (i in 1:nsubN2) {
    if (verbose) {
      cat(paste("Starting sub phase", i, "\n"))
    }

    # TODO PARALLEL: sample n chains and pass averaged parameters 
    # (but real states, caches) to next sub phase simulations
    if (parallel) {
      sfExport(
        list = c(
          "iterations",
          "thinningN2",
          "dep.var",
          "allowLoops",
          "effects",
          "parameters",
          "gain",
          "sensitivityVector",
          "multinomialProposal",
          "observedStatistics"
        )
      )
      res <- sfLapply(res, function(r) {
        runSubphase2(
          dep.var,
          state = r$state,
          cache = r$cache,
          effects,
          parameters,
          sensitivityVector,
          observedStatistics,
          gain,
          thinningN2,
          iterations,
          multinomialProposal = multinomialProposal,
          allowLoops = allowLoops,
          verbose = F
        )
      })
      parameterList <- lapply(res, "[[", "parameters")
      parameters <- colMeans(Reduce(rbind, parameterList))
      names(parameters) <- effects$name
    } else {
      resSub2 <- runSubphase2(
        dep.var,
        state,
        cache,
        effects,
        parameters,
        sensitivityVector,
        observedStatistics,
        gain,
        thinningN2,
        iterations,
        multinomialProposal = multinomialProposal,
        allowLoops = allowLoops,
        verbose = verbose
      )
      state <- resSub2$state
      cache <- resSub2$cache
      parameters <- resSub2$parameters
      names(parameters) <- effects$name
    }

    if (verbose) {
      cat(paste("New parameters:\n"))
    }
    if (verbose) {
      cat(paste(names(parameters), "\n", parameters, "\n"))
    }

    # determine number of iterations
    iterations <- iterations * 1.75
    # determine gain
    gain <- gain / 2
  }

  if (parallel) {
    sfStop()
  }

  return(parameters)
}


# runPhase3
runPhase3 <- function(dep.var,
                      state,
                      cache,
                      effects,
                      parameters,
                      observedStatistics,
                      iterationsN3,
                      burnInN3,
                      thinningN3,
                      parallel,
                      cpus,
                      allowLoops,
                      verbose = F,
                      returnDeps = F,
                      multinomialProposal = F,
                      fish = fish) {
  # simulate statistic matrix
  # if parallel computing, initialize several simulation chains and rbind the results at the end
  if (parallel) {
    iterationsPerCPU <- rep(iterationsN3 %/% cpus, cpus)
    rest <- iterationsN3 %% cpus
    if (rest > 0) {
      iterationsPerCPU[1:rest] <- iterationsPerCPU[1:rest] + 1
    }

    sfInit(parallel = T, cpus = cpus)
    # TODO. Replace this long command with sfLibrary("NetDist") once the package is packaged
    sfLibrary("MoNAn", character.only=TRUE)

    statsA <-
      sfLapply(iterationsPerCPU, function(nIt) {
        simulateStatisticVectors(
          dep.var,
          state,
          cache,
          effects,
          parameters,
          burnIn = burnInN3,
          iterations = nIt,
          thinning = thinningN3,
          allowLoops,
          verbose = verbose,
          returnDeps = returnDeps,
          multinomialProposal = multinomialProposal,
          fish = F
        )
      })


    if (returnDeps) {
      stats <- list()
      stats[[1]] <- lapply(statsA, function(x) {
        x[[1]]
      })
      stats[[2]] <-
        unlist(lapply(statsA, function(x) {
          x[[2]]
        }), recursive = F)
      statisticsMatrix <- Reduce("rbind", stats[[1]])
    } else {
      stats <- statsA
      statisticsMatrix <- Reduce("rbind", stats)
    }

    sfStop()
  } else {
    stats <- simulateStatisticVectors(
      dep.var,
      state,
      cache,
      effects,
      parameters,
      burnIn = burnInN3,
      iterations = iterationsN3,
      thinning = thinningN3,
      allowLoops,
      verbose = verbose,
      returnDeps = returnDeps,
      multinomialProposal = multinomialProposal,
      fish = fish
    )
    if (returnDeps) {
      statisticsMatrix <- stats[[1]]
    } else {
      statisticsMatrix <- stats
    }
  } # no parallel computing

  # calculate convergence statistics
  simulatedMeans <- colMeans(statisticsMatrix)
  simulatedSDs <- apply(statisticsMatrix, 2, sd)
  convergenceStatistic <-
    (simulatedMeans - observedStatistics) / simulatedSDs

  # calculate covariance matrix
  covarianceMatrix <- getCovarianceMatrix(statisticsMatrix)

  # calculate covariance matrix of parameters
  covarianceMatrixParameters <- solve(covarianceMatrix)

  # create empty object for the deps
  deps <- NULL

  if (returnDeps) {
    deps <- stats[[2]]
  }

  # returns covariance matrix and convergence statistics
  return(
    list(
      covarianceMatrix = covarianceMatrixParameters,
      convergenceStatistics = convergenceStatistic,
      deps = deps
    )
  )
}


# runSubphase2
runSubphase2 <- function(dep.var,
                         state,
                         cache,
                         effects,
                         initialParameters,
                         sensitivityVector,
                         observedStatistics,
                         gain,
                         thinningN2,
                         iterations,
                         multinomialProposal = F,
                         allowLoops,
                         verbose = F) {
  parameters <- c()

  for (i in 1:iterations) {
    if (verbose) {
      cat(".")
    }
    res <-
      simulateNSteps(
        dep.var,
        state,
        cache,
        effects,
        initialParameters,
        allowLoops,
        multinomialProposal = multinomialProposal,
        n = (thinningN2 + 1)
      )
    state <- res$state
    cache <- res$cache

    # calculate statistics
    statistics <-
      getNetworkStatistics(dep.var, state, cache, effects)

    # calculate and store new parameters
    updatedParameters <-
      initialParameters - gain * sensitivityVector * (statistics - observedStatistics)
    initialParameters <- updatedParameters
    parameters <- rbind(parameters, updatedParameters)
  }
  if (verbose) {
    cat("\n")
  }

  # return mean of the parameters
  return(list(
    parameters = colMeans(parameters),
    state = state,
    cache = cache
  ))
}


# simulateNSteps
simulateNSteps <-
  function(dep.var,
           state,
           cache,
           effects,
           parameters,
           allowLoops = T,
           senderFixed = T,
           receiverFixed = F,
           multinomialProposal = F,
           debug = F,
           n = 1) {
    for (i in 1:n) {
      res <-
        simulateOneStep(
          dep.var,
          state,
          cache,
          effects,
          parameters,
          allowLoops,
          senderFixed,
          receiverFixed,
          multinomialProposal,
          debug
        )
      state <- res$state
      cache <- res$cache
    }
    return(res)
  }


# simulateOneStep
simulateOneStep <-
  function(dep.var,
           state,
           cache,
           effects,
           parameters,
           allowLoops,
           senderFixed = T,
           receiverFixed = F,
           multinomialProposal = F,
           debug = F) {
    resourceName <- state[[dep.var]]$nodeSet[3]
    nodeName <- state[[dep.var]]$nodeSet[1]
    resourceSet <-
      which(state[[resourceName]]$considerWhenSampling)

    # sample edges to swap
    randomEdge <- sample(resourceSet, size = 1)
    i <- state[[dep.var]]$data[randomEdge, 1]
    j <- state[[dep.var]]$data[randomEdge, 2]

    # nEdges <- state[[dep.var]]$size[1]
    nodeSet <-
      state[[nodeName]]$ids[state[[nodeName]]$considerWhenSampling]
    if (!allowLoops) {
      nodeSet <- setdiff(nodeSet, i)
    }

    ## MULTINOMIAL PROPOSAL ## CONSIDER REFACTOR TO FUNCTION

    if (multinomialProposal) {
      if (debug) {
        print(paste("Proposed multinomial change from", i, j))
      }

      # remove resource from process state
      effectFunctions <- effects$effectFormulas
      statisticsDrop <-
        unlist(lapply(effectFunctions, function(f) {
          f(
            dep.var = dep.var,
            state = state,
            cache = cache,
            i = i,
            j = j,
            update = -1,
            edge = randomEdge
          )
        }))

      # update process state and cache
      cacheNew <- cache
      cacheNew[[dep.var]] <-
        updateWeightedCache(
          cache[[dep.var]],
          i,
          j,
          resourceID = randomEdge,
          state = state,
          dep.var = dep.var,
          update = -1
        )
      stateNew <- state
      stateNew[[dep.var]]$data[randomEdge, ] <- rep(NA, 2)

      # for each possible receiver
      statisticsCreate <-
        sapply(nodeSet, function(k) {
          unlist(lapply(effectFunctions, function(f) {
            f(
              dep.var = dep.var,
              state = state,
              cache = cache,
              i = i,
              j = k,
              update = 1,
              edge = randomEdge
            )
          }))
        },
        simplify = T
        )
      objValues <-
        colSums((statisticsCreate + statisticsDrop) * parameters)
      pChange <- exp(objValues) / sum(exp(objValues))

      # draw receiver node
      draw <- which(cumsum(pChange) >= runif(1))[1]
      l <- nodeSet[draw]

      # update process state and cache
      cache[[dep.var]] <-
        updateWeightedCache(
          cacheNew[[dep.var]],
          sender = i,
          receiver = l,
          resourceID = randomEdge,
          state = state,
          dep.var = dep.var,
          update = 1
        )
      state[[dep.var]]$data[randomEdge, ] <- c(i, l)
    }


    ## BINOMIAL PROPOSAL ## CONSIDER REFACTOR TO FUNCTION

    if (!multinomialProposal) {
      if (senderFixed) {
        k <- i
      } # TODO else define k differently

      l <-
        sample(nodeSet, size = 1) # TODO consider case of receiver fixed

      if (debug) {
        print(paste("Proposed change from", i, j, "to", k, l))
      }


      if (j == l) {
        return(list(state = state, cache = cache))
      }

      effectFunctions <- effects$effectFormulas
      statisticsDrop <-
        unlist(lapply(effectFunctions, function(f) {
          f(
            dep.var = dep.var,
            state = state,
            cache = cache,
            i = i,
            j = j,
            update = -1,
            edge = randomEdge
          )
        }))

      # update process state and cache
      cacheNew <- cache
      cacheNew[[dep.var]] <-
        updateWeightedCache(
          cache[[dep.var]],
          i,
          j,
          resourceID = randomEdge,
          state = state,
          dep.var = dep.var,
          update = -1
        )
      stateNew <- state
      stateNew[[dep.var]]$data[randomEdge, ] <- rep(NA, 2)

      statisticsCreate <-
        unlist(lapply(effectFunctions, function(f) {
          f(
            dep.var = dep.var,
            state = stateNew,
            cache = cacheNew,
            i = k,
            j = l,
            update = +1,
            edge = randomEdge
          )
        }))

      changeContribution <-
        sum((statisticsDrop + statisticsCreate) * parameters)

      acceptanceProbability <- min(1, exp(changeContribution))

      if (acceptanceProbability >= runif(1)) {
        cache[[dep.var]] <-
          updateWeightedCache(
            cacheNew[[dep.var]],
            sender = k,
            receiver = l,
            resourceID = randomEdge,
            state = state,
            dep.var = dep.var,
            update = 1
          )
        state[[dep.var]]$data[randomEdge, ] <- c(k, l)
      }
    }

    return(list(state = state, cache = cache))
  }


# simulateStatisticVectors
simulateStatisticVectors <- function(dep.var,
                                     state,
                                     cache,
                                     effects,
                                     initialParameters,
                                     burnIn,
                                     iterations,
                                     thinning,
                                     allowLoops,
                                     verbose = F,
                                     multinomialProposal = F,
                                     fish = F,
                                     returnDeps = F) {
  statisticsMatrix <- c()

  if (returnDeps) {
    deps <- list()
  }

  # burn in
  if (verbose) {
    cat(paste("Starting burn-in with", burnIn, "steps\n"))
  }
  res <-
    simulateNSteps(
      dep.var,
      state,
      cache,
      effects,
      initialParameters,
      allowLoops,
      multinomialProposal = multinomialProposal,
      n = burnIn
    )
  cache <- res$cache
  state <- res$state

  # thinning and calculate statistics
  for (i in 1:iterations) {
    if (fish) {
      int <- i %% 11
      s <- "_,.-'``'-.,"
      cat(substr(s, int + 1, int + 1))
      if (runif(1) < 0.02) {
        cat("><(((A*>")
      }
    } else if (verbose) {
      cat(".")
    }
    res <-
      simulateNSteps(
        dep.var,
        state,
        cache,
        effects,
        initialParameters,
        allowLoops,
        multinomialProposal = multinomialProposal,
        n = thinning + 1
      )
    state <- res$state
    cache <- res$cache

    # calculate and save statistics
    stats <-
      getNetworkStatistics(dep.var,
        state = state,
        cache = cache,
        effects = effects
      )
    statisticsMatrix <- rbind(statisticsMatrix, stats)

    if (returnDeps) {
      deps[[i]] <- list(state = state, cache = cache)
    }
  }

  if (fish) {
    cat("\n")
  } else if (verbose) {
    cat("\n")
  }


  if (returnDeps) {
    return(list(statisticsMatrix, deps))
  } else {
    return(statisticsMatrix)
  }
}


# updateWeightedCache
updateWeightedCache <- function(cache,
                                sender,
                                receiver,
                                resourceID = NULL,
                                state = NULL,
                                dep.var = NULL,
                                # refers to the element of the cache to be updated
                                update,
                                debug = FALSE) {
  if (debug) {
    print(paste("update cache for s-r-u", sender, receiver, update))
  }

  # Updates of weighted resource caches
  resourceCovariates <- names(cache$resourceNetworks)
  for (ressCovar in resourceCovariates) {
    v <- cache$resourceNetworks[[ressCovar]][sender, receiver]
    cache$resourceNetworks[[ressCovar]][sender, receiver] <-
      v + update * state[[ressCovar]]$data[resourceID]
  }


  # The following uses information about the net flow network before the update
  if (update > 0 && cache$netFlowsNetwork[sender, receiver] < 0) {
    cache$minNetwork[sender, receiver] <-
      cache$minNetwork[sender, receiver] +
      # TODO: allow non-1 updates
      # min(- cache$netFlowsNetwork[sender, receiver], update)
      1
  }
  if (update < 0 &&
    cache$netFlowsNetwork[sender, receiver] <= 0 &&
    sender != receiver) {
    cache$minNetwork[sender, receiver] <-
      cache$minNetwork[sender, receiver] -
      # TODO: allow updates of more / less than one
      1
  }
  if (cache$minNetwork[sender, receiver] < 0) {
    stop(
      paste(
        "Error in cache update, negative min tie.",
        sender,
        receiver,
        update,
        cache$netFlowsNetwork[sender, receiver]
      )
    )
  }
  ## cache$minNetwork[sender, receiver] <- cache$minNetwork[sender, receiver] + max(0, min(- cache$netFlowsNetwork[sender, receiver], update) )
  cache$minNetwork[receiver, sender] <-
    cache$minNetwork[sender, receiver]

  # weighted updates
  cache$valuedNetwork[sender, receiver] <-
    cache$valuedNetwork[sender, receiver] + update
  cache$netFlowsNetwork[sender, receiver] <-
    cache$netFlowsNetwork[sender, receiver] + update
  cache$netFlowsNetwork[receiver, sender] <-
    cache$netFlowsNetwork[receiver, sender] - update

  cache
}
