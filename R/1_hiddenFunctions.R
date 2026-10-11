########## hiddenFunctions

# binarizeNetwork
binarizeNetwork <- function(network) {
  net <- 1 * (network >= 1)
  diag(net) <- 0
  net
}


# checkProcessState
checkProcessState <- function(state) {
  
  # extract dependent variable
  dep.var <- state$dep.var
  
  # make sure that dependent variable has corresponding edgelist
  if (!(dep.var %in% names(state))) {
    stop("Dependent variable has no corresponding edgelist in the state object.")
  }
  
  # is edgelist of class 'edgelist.monan'?
  if (!(inherits(state[[dep.var]], "edgelist.monan"))) {
    stop("Dependent variable is not of class 'edgelist.monan'.")
  }
  
  
  # extract nodeset names defined in edgelist
  nodesets <- state[[dep.var]]$nodeSet
  
  # do nodeset names from edgelist have a corresponding nodeset and are they of class 'nodeset.monan'?
  for (i in 1:length(nodesets)) {
    if (!(nodesets[i] %in% names(state))) {
      stop(paste(nodesets[i], "has no corresponding object.",
                 "Nodesets specified in the edgelist must have a corresponding object with the same name."))
    }
    if (!(inherits(state[[nodesets[i]]], "nodeSet.monan"))) {
      stop(paste(nodesets[i], "is not of class 'nodeSet.monan'."))
    }
  }
  
  # do dependent variable and individual nodeset have the same number of observations?
  if (!(dim(state[[dep.var]]$data)[1] == length(state[[nodesets[3]]]$ids))) {
    stop("Dependent variable and individual nodeset do not have the same number of observations.")
  }
  
  # extract covar names
  covars <- names(state)[!(names(state) %in% c(dep.var, nodesets, "dep.var"))]
  
  # if there are covars, check them
  if (length(covars) != 0) {
    
    # are covars of class 'nodeVar.monan' or 'network.monan'? do they have correctly
    # specified nodeset names and are these nodesets of the correct size?
    for (i in 1:length(covars)) {
      #class
      if (!(inherits(state[[covars[i]]], "nodeVar.monan") || 
            inherits(state[[covars[i]]], "network.monan"))) {
        stop(paste(covars[i], "must either be of class 'nodeVar.monan' or 'network.monan'."))
      }
      # nodeset names
      if (!(all(state[[covars[i]]]$nodeSet %in% names(state)))) {
        stop(paste("The nodeset of covar '", covars[i], "' has no corresponding object.",
                   "Nodesets assigned when creating covars must have a corresponding object with the same name."           ))
      }
      # nodeset sizes
      if (!(state[[covars[i]]]$size[1] == length(state[[state[[covars[i]]]$nodeSet[1]]]$ids))) {
        stop(paste("The covar '", covars[i], "' is not of the same size as its assigned nodeset."))
      }
    }
  }
}


# createInternalCache
createInternalCache <-
  function(processState,
           resourceCovariates = NULL) {
    cache <- list()
    cacheObjectNames <- processState$dep.var
    
    for (name in cacheObjectNames) {
      if (!(class(processState[[name]]) %in% c("network.monan", "edgelist.monan"))) {
        stop(paste(name, "is not a network or edgelist."))
      }
      
      nodeSet1 <- processState[[processState[[name]]$nodeSet[1]]]$ids
      nodeSet2 <- processState[[processState[[name]]$nodeSet[2]]]$ids
      nActors1 <- length(nodeSet1)
      nActors2 <- length(nodeSet2)
      
      cache[[name]] <- list()
      if (is(processState[[name]], "network.monan")) {
        cache[[name]]$valuedNetwork <- processState[[name]]$data
      }
      if (is(processState[[name]], "edgelist.monan")) {
        # create valued network from edge list
        m <- matrix(0, nActors1, nActors2)
        
        # create weighted resource networks
        m.resource <-
          lapply(resourceCovariates, function(v) {
            matrix(0, nrow = nActors1, ncol = nActors2)
          })
        names(m.resource) <- resourceCovariates
        
        for (i in 1:processState[[name]]$size[1]) {
          sender <- processState[[name]]$data[i, 1]
          receiver <- processState[[name]]$data[i, 2]
          v <- m[sender, receiver]
          m[sender, receiver] <- v + 1
          
          # cache network * resource covariate matrices
          for (ressCovar in resourceCovariates) {
            v <- m.resource[[ressCovar]][sender, receiver]
            m.resource[[ressCovar]][sender, receiver] <- v +
              processState[[ressCovar]]$data[i]
          }
        }
        cache[[name]]$valuedNetwork <- m
        cache[[name]]$resourceNetworks <- m.resource
        
        # edge ids
        edgeSet <- processState[[processState[[name]]$nodeSet[3]]]$ids
        nEdges <- length(edgeSet)
      }
      if (nActors1 == nActors2) {
        cache[[name]]$netFlowsNetwork <-
          cache[[name]]$valuedNetwork - t(cache[[name]]$valuedNetwork)
        cache[[name]]$minNetwork <-
          matrix(
            mapply(min, cache[[name]]$valuedNetwork, t(cache[[name]]$valuedNetwork)),
            nActors1,
            nActors2
          )
        diag(cache[[name]]$minNetwork) <- 0
      }
    } # end for loop
    
    cache
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

# getNetStatsFromDeps
getNetStatsFromDeps <- function(dep.var, oneDep, ans, effects){
  depState <- ans$state
  depState[[dep.var]] <- oneDep
  
  resCovsInCache <- names(ans$cache[[dep.var]]$resourceNetworks)
  if(is.null(resCovsInCache)){
    depCache <- createInternalCache(depState)
  } else {
    depCache <- createInternalCache(depState, resourceCovariates = resCovsInCache)
  }
  getNetworkStatistics(dep.var, depState, depCache, effects)
}

# getNetworkStatistics
getNetworkStatistics <- function(dep.var, state, cache, effects) {
  if (useCppSampler(effects, state)) {
    return(cppNetworkStatistics(dep.var, state, cache, effects))
  }

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
            getTargetContribution = TRUE
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
                      multinomialProposal = FALSE,
                      allowLoops,
                      verbose = FALSE,
                      parallel,
                      cpus) {
  # simulate statistic matrix
  
  # initialize parallel computing
  if (parallel && cpus > 1) {
    sfInit(parallel = TRUE, cpus = cpus)
    # TODO. Replace this long command with sfLibrary("NetDist") once the package is packaged
    sfLibrary("MoNAn", character.only=TRUE)
  } else {
    parallel <- FALSE
  }
  
  if(parallel){
    cat(paste("Phase 1: \n burn-in", burnInN1, "steps\n", 
              iterationsN1, " iterations\n thinning", 
              thinningN1, "\n", cpus, "cpus\n"))
    
    statisticsMatrix <- sfLapply(1:cpus, simulateStatisticVectors, 
                                 dep.var = dep.var, 
                                 state = state, 
                                 cache = cache, 
                                 effects = effects, 
                                 initialParameters = initialParameters, 
                                 burnIn = burnInN1, 
                                 iterations = round(iterationsN1 / cpus), 
                                 thinning = thinningN1, 
                                 multinomialProposal = multinomialProposal, 
                                 allowLoops = allowLoops, 
                                 verbose = verbose)
    statisticsMatrix <- Reduce("rbind", statisticsMatrix)
    sfStop()
  } else {
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
  }

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
                      multinomialProposal = FALSE,
                      allowLoops,
                      verbose = FALSE) {
  # calculate observed statistics
  observedStatistics <-
    getNetworkStatistics(dep.var, state, cache, effects)

  # initialize parallel computing
  if (parallel && cpus > 1) {
    sfInit(parallel = TRUE, cpus = cpus)
    # TODO. Replace this long command with sfLibrary("NetDist") once the package is packaged
    sfLibrary("MoNAn", character.only=TRUE)
  } else {
    parallel <- FALSE
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
      cat(paste("\n"))
      cat(paste0("Sub phase", i, ":\n burn-in ", burnInN2, " steps\n ", 
                iterations, " iterations\n thinning ", 
                thinningN2, "\n ", cpus, " cpus\n"))
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
          verbose = FALSE
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
      cat(paste("\n"))
      cat(paste("New parameters:\n"))
    }
    if (verbose) {
      cat(paste(names(parameters), "\n", parameters, "\n"))
    }

    # determine number of iterations
    iterations <- round(iterations * 1.75)
    # determine gain
    gain <- gain / 2
  }

  if (parallel) {
    sfStop()
  }

  return(parameters)
}


# runPhase2GT
# Alternative to runPhase2: Geyer-Thompson MCMC maximum likelihood estimation
# (Geyer & Thompson 1992; as in MCMC MLE for ERGMs). Networks are simulated at
# the current parameters theta0, and the parameters are updated with
# geyerThompsonStep. Iterations with iterationsGT draws continue until a step
# is an unrestricted Geyer-Thompson step, i.e. the observed statistics are
# within the range of the simulated ones; the estimate is then taken from a
# final step based on finalIterationsGT draws, which reduces Monte Carlo error.
runPhase2GT <- function(dep.var,
                        state,
                        cache,
                        effects,
                        initialParameters,
                        observedStatistics,
                        burnIn,
                        thinning,
                        iterationsGT,
                        maxIterationsGT,
                        finalIterationsGT = 3 * iterationsGT,
                        minESS = 0.5,
                        parallel,
                        cpus,
                        multinomialProposal = FALSE,
                        allowLoops,
                        verbose = FALSE) {
  if (parallel && cpus > 1) {
    sfInit(parallel = TRUE, cpus = cpus)
    sfLibrary("MoNAn", character.only = TRUE)
    on.exit(sfStop(), add = TRUE)
  } else {
    parallel <- FALSE
    cpus <- 1
  }

  splitDraws <- function(n) {
    draws <- rep(n %/% cpus, cpus)
    rest <- n %% cpus
    if (rest > 0) {
      draws[1:rest] <- draws[1:rest] + 1
    }
    draws
  }

  runChain <- function(chain, dep.var, effects, parameters, burnIn, thinning,
                       allowLoops, multinomialProposal) {
    res <- simulateStatisticVectors(
      dep.var, chain$state, chain$cache, effects, parameters,
      burnIn = burnIn, iterations = chain$draws, thinning = thinning,
      allowLoops = allowLoops, multinomialProposal = multinomialProposal,
      returnState = TRUE
    )
    list(statistics = res$statistics, state = res$state, cache = res$cache)
  }

  # simulate a sample of statistics at the given parameters; chains continue
  # from where the previous sample ended
  chains <- lapply(seq_len(cpus), function(k) list(state = state, cache = cache))
  simulateSample <- function(parameters, n) {
    draws <- splitDraws(n)
    for (k in seq_len(cpus)) {
      chains[[k]]$draws <<- draws[k]
    }
    applyFun <- if (parallel) sfLapply else lapply
    chains <<- applyFun(chains, runChain, dep.var = dep.var, effects = effects,
                        parameters = parameters, burnIn = burnIn,
                        thinning = thinning, allowLoops = allowLoops,
                        multinomialProposal = multinomialProposal)
    Reduce(rbind, lapply(chains, "[[", "statistics"))
  }

  report <- function(label, statistics, step, parameters) {
    if (!verbose) {
      return(invisible(NULL))
    }
    tRatios <- (colMeans(statistics) - observedStatistics) / apply(statistics, 2, sd)
    cat(paste0("\n", label, ": max |t-ratio| ", round(max(abs(tRatios)), 3),
               ", effective sample size ",
               if (is.na(step$ess)) "-" else round(step$ess), " of ", nrow(statistics),
               ", ", step$method, " step",
               if (step$restricted) " (restricted)" else "", "\n"))
    cat(paste(names(parameters), "\n", parameters, "\n"))
  }

  parameters <- initialParameters
  names(parameters) <- effects$name

  for (iter in 1:maxIterationsGT) {
    statistics <- simulateSample(parameters, iterationsGT)
    step <- geyerThompsonStep(statistics, observedStatistics, minESS)
    parameters <- parameters + step$delta
    report(paste("Geyer-Thompson iteration", iter), statistics, step, parameters)
    if (!step$restricted) {
      break
    }
  }

  # final step from a larger sample
  statistics <- simulateSample(parameters, finalIterationsGT)
  step <- geyerThompsonStep(statistics, observedStatistics, minESS)
  parameters <- parameters + step$delta
  report("Geyer-Thompson final step", statistics, step, parameters)

  parameters
}


# geyerThompsonStep
# One parameter update from a sample of statistics simulated at theta0.
# If the observed statistics lie within the range of the sample (measured by
# Mahalanobis distance), the Geyer-Thompson approximation of
# l(theta0 + delta) - l(theta0), -log mean_k exp(delta . (s_k - s_obs)), is
# maximised with Newton steps, keeping the effective sample size of the
# importance weights >= minESS * M. Otherwise, importance sampling is
# unreliable and a step based on the normal approximation of the statistics
# is taken towards a target between the sample mean and the observed
# statistics that lies within the sample range (similar to the step length
# of Hummel, Hunter & Handcock 2012).
geyerThompsonStep <- function(statistics, observedStatistics, minESS = 0.5,
                              maxNewton = 100, tol = 1e-8) {
  M <- nrow(statistics)
  p <- ncol(statistics)
  mu <- colMeans(statistics)
  C <- cov(statistics)
  solveC <- function(A, b) {
    tryCatch(solve(A, b), error = function(e) {
      solve(A + diag(1e-8 * max(1, max(diag(A))), p), b)
    })
  }

  # range of the sample: Mahalanobis distances of the draws from their mean
  centred <- sweep(statistics, 2, mu)
  sampleDistances <- sqrt(pmax(0, rowSums(centred * t(solveC(C, t(centred))))))
  range <- as.numeric(quantile(sampleDistances, 0.9))
  obsDistance <- sqrt(max(0, sum((observedStatistics - mu) *
                                   solveC(C, observedStatistics - mu))))

  if (obsDistance > range) {
    gamma <- range / obsDistance
    delta <- as.vector(solveC(C, gamma * (observedStatistics - mu)))
    return(list(delta = delta, ess = NA, restricted = TRUE, method = "normal"))
  }

  d <- sweep(statistics, 2, observedStatistics)
  weightsOf <- function(delta) {
    eta <- as.vector(d %*% delta)
    w <- exp(eta - max(eta))
    w / sum(w)
  }
  # objective, computed stably: -log mean exp(eta)
  objective <- function(delta) {
    eta <- as.vector(d %*% delta)
    -(max(eta) + log(mean(exp(eta - max(eta)))))
  }
  essOf <- function(w) 1 / sum(w^2)

  delta <- rep(0, p)
  restricted <- FALSE
  for (k in 1:maxNewton) {
    w <- weightsOf(delta)
    m <- colSums(w * d)
    Cw <- crossprod(d * sqrt(w)) - tcrossprod(m)
    newtonStep <- -solveC(Cw, m)

    # backtrack until the objective does not decrease and the effective
    # sample size is large enough
    a <- 1
    repeat {
      candidate <- delta + a * newtonStep
      if (objective(candidate) >= objective(delta) - 1e-12 &&
          essOf(weightsOf(candidate)) >= minESS * M) {
        break
      }
      a <- a / 2
      if (a < 1e-6) {
        candidate <- delta
        break
      }
    }
    if (a < 1) {
      restricted <- restricted ||
        essOf(weightsOf(candidate)) < minESS * M * 1.01
    }
    change <- max(abs(candidate - delta))
    delta <- candidate
    if (change < tol) {
      break
    }
  }

  list(delta = delta, ess = essOf(weightsOf(delta)), restricted = restricted,
       method = "Geyer-Thompson")
}


# getMPLE
# Maximum pseudo-likelihood estimate: a conditional logit of each resource's
# observed destination against all possible destinations, given the
# destinations of all other resources. Used as starting values for the
# Geyer-Thompson estimation. Returns NULL if the estimate does not exist
# (e.g. separation) or the optimisation fails.
getMPLE <- function(dep.var, state, cache, effects, allowLoops,
                    maxNewton = 50, tol = 1e-8) {
  state <- toMutableState(state, dep.var)
  cache <- toMutableCache(cache)
  edgelist <- state[[dep.var]]
  nodeName <- edgelist$nodeSet[1]
  resourceName <- edgelist$nodeSet[3]
  nodes <- state[[nodeName]]$ids[state[[nodeName]]$considerWhenSampling]
  resources <- which(state[[resourceName]]$considerWhenSampling)
  effectFunctions <- effects$effectFormulas
  p <- length(effectFunctions)

  # for each resource: change statistics of adding it to every possible
  # destination after removing it, and which destination is observed
  choiceSets <- lapply(resources, function(e) {
    i <- edgelist$data[e, 1]
    j <- edgelist$data[e, 2]
    destinations <- if (allowLoops) nodes else setdiff(nodes, i)

    cells <- getCacheCells(cache[[dep.var]], i, j)
    updateCacheInPlace(cache[[dep.var]], i, j, resourceID = e, state = state, update = -1)
    setEdgeInPlace(edgelist, e, rep(NA, 2))
    X <- t(vapply(destinations, function(k) {
      vapply(effectFunctions, function(f) {
        f(dep.var = dep.var, state = state, cache = cache,
          i = i, j = k, update = 1, edge = e)
      }, numeric(1))
    }, numeric(p)))
    restoreCacheCells(cache[[dep.var]], i, j, cells)
    setEdgeInPlace(edgelist, e, c(i, j))

    list(X = matrix(X, ncol = p), observed = match(j, destinations))
  })

  theta <- rep(0, p)
  logPL <- function(theta) {
    sum(vapply(choiceSets, function(cs) {
      eta <- as.vector(cs$X %*% theta)
      eta[cs$observed] - max(eta) - log(sum(exp(eta - max(eta))))
    }, numeric(1)))
  }
  for (k in 1:maxNewton) {
    gradient <- rep(0, p)
    hessian <- matrix(0, p, p)
    for (cs in choiceSets) {
      eta <- as.vector(cs$X %*% theta)
      prob <- exp(eta - max(eta))
      prob <- prob / sum(prob)
      expected <- colSums(prob * cs$X)
      gradient <- gradient + cs$X[cs$observed, ] - expected
      hessian <- hessian - (crossprod(cs$X * sqrt(prob)) - tcrossprod(expected))
    }
    step <- tryCatch(-solve(hessian, gradient), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      return(NULL)
    }
    # step halving to keep the pseudo-likelihood increasing
    a <- 1
    current <- logPL(theta)
    while (logPL(theta + a * step) < current - 1e-10 && a > 1e-6) {
      a <- a / 2
    }
    theta <- theta + a * step
    if (max(abs(a * step)) < tol) {
      break
    }
  }
  if (any(!is.finite(theta)) || max(abs(theta)) > 1e3) {
    return(NULL)
  }
  names(theta) <- effects$name
  theta
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
                      verbose = FALSE,
                      returnDeps = FALSE,
                      multinomialProposal = FALSE,
                      fish = fish) {
  # simulate statistic matrix
  # if parallel computing, initialize several simulation chains and rbind the results at the end
  if (parallel) {
    iterationsPerCPU <- rep(iterationsN3 %/% cpus, cpus)
    rest <- iterationsN3 %% cpus
    if (rest > 0) {
      iterationsPerCPU[1:rest] <- iterationsPerCPU[1:rest] + 1
    }

    sfInit(parallel = TRUE, cpus = cpus)
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
          fish = FALSE
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
        }), recursive = FALSE)
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
  
  # extract only the dependent var to return from deps
  if (returnDeps) {
    deps <- lapply(stats[[2]], function(x){
      x$state[[dep.var]]
    })
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
                         multinomialProposal = FALSE,
                         allowLoops,
                         verbose = FALSE) {
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

  ### return mean of the parameters discarding the 
  ### initial dk that are clearly before convergence
  # find the smaller of two numbers: 
  # 25% of the number of iterations OR
  # when the expected residual of the step-size is less than 0.05
  dk1 <- round(iterations/4)
  
  dk2 <- 1
  res <- (1-gain)
  
  repeat{
    if(dk2 > 100) break
    if(any(res < 0.05)) break
    dk2 <- dk2 + 1
    res <- res*(1-gain)
  }
  
  dk <- min(dk1, dk2)
  
  if(ncol(parameters) == 1){
    return(list(
      parameters = mean(parameters[dk:iterations,]),
      state = state,
      cache = cache
    ))
  }
  
  return(list(
    parameters = colMeans(parameters[dk:iterations,]),
    state = state,
    cache = cache
  ))
}


# simulateNSteps
# returns state and cache in their in-place form (see toMutableState);
# convert with toListState / toListCache before storing them
simulateNSteps <-
  function(dep.var,
           state,
           cache,
           effects,
           parameters,
           allowLoops = TRUE,
           senderFixed = TRUE,
           receiverFixed = FALSE,
           multinomialProposal = FALSE,
           debug = FALSE,
           n = 1) {
    if (senderFixed && !debug && useCppSampler(effects, state)) {
      res <- simulateChain(dep.var, state, cache, effects, parameters,
                           burnIn = n, iterations = 0, thinning = 0,
                           allowLoops = allowLoops,
                           multinomialProposal = multinomialProposal)
      return(list(state = res$state, cache = res$cache))
    }

    state <- toMutableState(state, dep.var)
    cache <- toMutableCache(cache)
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
# changes state and cache in place; both must be in their in-place form
simulateOneStep <-
  function(dep.var,
           state,
           cache,
           effects,
           parameters,
           allowLoops,
           senderFixed = TRUE,
           receiverFixed = FALSE,
           multinomialProposal = FALSE,
           debug = FALSE) {
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

    effectFunctions <- effects$effectFormulas

    ## MULTINOMIAL PROPOSAL ## CONSIDER REFACTOR TO FUNCTION

    if (multinomialProposal) {
      if (debug) {
        print(paste("Proposed multinomial change from", i, j))
      }

      # remove resource from process state
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
      updateCacheInPlace(cache[[dep.var]], i, j,
                         resourceID = randomEdge, state = state, update = -1)
      setEdgeInPlace(state[[dep.var]], randomEdge, rep(NA, 2))

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
        simplify = TRUE
        )
      # one row per effect, also if there is only one effect (sapply then
      # returns a vector)
      statisticsCreate <- matrix(statisticsCreate, nrow = length(effectFunctions))
      objValues <-
        colSums((statisticsCreate + statisticsDrop) * parameters)
      pChange <- exp(objValues) / sum(exp(objValues))

      # draw receiver node
      draw <- which(cumsum(pChange) >= runif(1))[1]
      l <- nodeSet[draw]

      # update process state and cache
      updateCacheInPlace(cache[[dep.var]], i, l,
                         resourceID = randomEdge, state = state, update = 1)
      setEdgeInPlace(state[[dep.var]], randomEdge, c(i, l))
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

      # remove the resource from state and cache, keeping what is needed to
      # undo this if the change is rejected
      oldCells <- getCacheCells(cache[[dep.var]], i, j)
      updateCacheInPlace(cache[[dep.var]], i, j,
                         resourceID = randomEdge, state = state, update = -1)
      setEdgeInPlace(state[[dep.var]], randomEdge, rep(NA, 2))

      statisticsCreate <-
        unlist(lapply(effectFunctions, function(f) {
          f(
            dep.var = dep.var,
            state = state,
            cache = cache,
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
        updateCacheInPlace(cache[[dep.var]], k, l,
                           resourceID = randomEdge, state = state, update = 1)
        setEdgeInPlace(state[[dep.var]], randomEdge, c(k, l))
      } else {
        restoreCacheCells(cache[[dep.var]], i, j, oldCells)
        setEdgeInPlace(state[[dep.var]], randomEdge, c(i, j))
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
                                     verbose = FALSE,
                                     multinomialProposal = FALSE,
                                     fish = FALSE,
                                     returnDeps = FALSE,
                                     returnState = FALSE) {
  # the whole chain in C++ if all effects have a C++ version
  if (useCppSampler(effects, state)) {
    if (verbose) {
      cat(paste("Starting burn-in with", burnIn, "steps\n"))
    }
    res <- simulateChain(dep.var, state, cache, effects, initialParameters,
                         burnIn = burnIn, iterations = iterations, thinning = thinning,
                         allowLoops = allowLoops, multinomialProposal = multinomialProposal,
                         returnDraws = returnDeps)
    statisticsMatrix <- res$statistics
    rownames(statisticsMatrix) <- rep("stats", nrow(statisticsMatrix))
    if (returnDeps) {
      return(list(statisticsMatrix, res$draws))
    }
    if (returnState) {
      return(list(statistics = statisticsMatrix, state = res$state, cache = res$cache))
    }
    return(statisticsMatrix)
  }

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
      deps[[i]] <- list(state = toListState(state, dep.var),
                        cache = toListCache(cache))
    }
  }

  if (fish) {
    cat("\n")
  } else if (verbose) {
    cat("\n")
  }


  if (returnDeps) {
    return(list(statisticsMatrix, deps))
  } else if (returnState) {
    # the end of the chain (in-place form), to continue it later
    return(list(statistics = statisticsMatrix, state = state, cache = cache))
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


########## in-place state and cache for MCMC chains
# Within one chain, the dependent variable of the state and the cache are
# stored in environments and updated in place. This avoids copying the edgelist
# and all cache matrices in every simulation step. Effects read them as before
# (state[[dep.var]]$data, cache[[dep.var]]$valuedNetwork[i, j], ...).
# Updates are done by functions whose environment is that environment, using
# <<-; e$x[i, j] <- v would copy x, and evalq() is not byte-compiled.
# Anything that leaves a chain is converted back with toListState/toListCache.

# toMutableState
toMutableState <- function(state, dep.var) {
  if (is.environment(state[[dep.var]])) {
    return(state)
  }
  e <- list2env(unclass(state[[dep.var]]), envir = new.env(parent = baseenv()))
  e$.names <- names(state[[dep.var]])
  e$.class <- class(state[[dep.var]])

  # placeholder for the variable the function below updates in e (its
  # environment is set to e below)
  data <- NULL
  e$.setEdge <- function(edge, value) {
    data[edge, ] <<- value
  }
  environment(e$.setEdge) <- e

  state[[dep.var]] <- e
  state
}

# toListState
toListState <- function(state, dep.var) {
  e <- state[[dep.var]]
  if (!is.environment(e)) {
    return(state)
  }
  l <- mget(e$.names, envir = e)
  class(l) <- e$.class
  state[[dep.var]] <- l
  state
}

# toMutableCache
toMutableCache <- function(cache) {
  lapply(cache, function(x) {
    if (is.environment(x)) {
      return(x)
    }
    e <- list2env(x, envir = new.env(parent = baseenv()))
    e$.names <- names(x)
    if (is.null(x$resourceNetworks)) {
      e$resourceNetworks <- NULL
    }

    # placeholders for the variables the functions below update in e (their
    # environment is set to e below)
    valuedNetwork <- netFlowsNetwork <- minNetwork <- resourceNetworks <- NULL

    # same logic as updateWeightedCache
    e$.update <- function(sender, receiver, resourceID, state, update) {
      # Updates of weighted resource caches
      for (ressCovar in names(resourceNetworks)) {
        resourceNetworks[[ressCovar]][sender, receiver] <<-
          resourceNetworks[[ressCovar]][sender, receiver] +
          update * state[[ressCovar]]$data[resourceID]
      }

      # The following uses information about the net flow network before the update
      if (update > 0 && netFlowsNetwork[sender, receiver] < 0) {
        minNetwork[sender, receiver] <<- minNetwork[sender, receiver] + 1
      }
      if (update < 0 && netFlowsNetwork[sender, receiver] <= 0 && sender != receiver) {
        minNetwork[sender, receiver] <<- minNetwork[sender, receiver] - 1
      }
      if (minNetwork[sender, receiver] < 0) {
        stop(paste("Error in cache update, negative min tie.",
                   sender, receiver, update, netFlowsNetwork[sender, receiver]))
      }
      minNetwork[receiver, sender] <<- minNetwork[sender, receiver]

      # weighted updates
      valuedNetwork[sender, receiver] <<- valuedNetwork[sender, receiver] + update
      netFlowsNetwork[sender, receiver] <<- netFlowsNetwork[sender, receiver] + update
      netFlowsNetwork[receiver, sender] <<- netFlowsNetwork[receiver, sender] - update
    }

    # the cells that an update of sender-receiver changes, to restore them
    # exactly if a proposed change is rejected (adding back a numeric resource
    # covariate value is not always exact in floating point)
    e$.getCells <- function(sender, receiver) {
      list(
        valued = valuedNetwork[sender, receiver],
        netFlows = c(netFlowsNetwork[sender, receiver], netFlowsNetwork[receiver, sender]),
        min = c(minNetwork[sender, receiver], minNetwork[receiver, sender]),
        resource = lapply(resourceNetworks, function(m) m[sender, receiver])
      )
    }

    e$.restoreCells <- function(sender, receiver, cells) {
      for (ressCovar in names(resourceNetworks)) {
        resourceNetworks[[ressCovar]][sender, receiver] <<- cells$resource[[ressCovar]]
      }
      valuedNetwork[sender, receiver] <<- cells$valued
      netFlowsNetwork[sender, receiver] <<- cells$netFlows[1]
      netFlowsNetwork[receiver, sender] <<- cells$netFlows[2]
      minNetwork[sender, receiver] <<- cells$min[1]
      minNetwork[receiver, sender] <<- cells$min[2]
    }

    environment(e$.update) <- e
    environment(e$.getCells) <- e
    environment(e$.restoreCells) <- e
    e
  })
}

# toListCache
toListCache <- function(cache) {
  lapply(cache, function(e) {
    if (!is.environment(e)) {
      return(e)
    }
    mget(e$.names, envir = e)
  })
}

# updateCacheInPlace
# same as updateWeightedCache, but changes the cache environment in place
updateCacheInPlace <- function(cacheEnv, sender, receiver, resourceID, state, update) {
  cacheEnv$.update(sender, receiver, resourceID, state, update)
  invisible(NULL)
}

# getCacheCells
getCacheCells <- function(cacheEnv, sender, receiver) {
  cacheEnv$.getCells(sender, receiver)
}

# restoreCacheCells
restoreCacheCells <- function(cacheEnv, sender, receiver, cells) {
  cacheEnv$.restoreCells(sender, receiver, cells)
  invisible(NULL)
}

# setEdgeInPlace
# sets the sender and receiver of one resource in the edgelist environment
setEdgeInPlace <- function(edgelistEnv, edge, value) {
  edgelistEnv$.setEdge(edge, value)
  invisible(NULL)
}
