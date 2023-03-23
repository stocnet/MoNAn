########## coreFunctions


# as.nodeVariable
#' Title
#'
#' @param values
#'
#' @return
#' @export
#'
#' @examples
as.nodeVariable <- function(values) {
  createNodeVariable(values)
}


# createEdgelist
#' Title
#'
#' @param el
#' @param nodeSet
#'
#' @return
#' @export
#'
#' @examples
createEdgelist <-
  function(el, nodeSet = c("actors", "actors", "edges")) {
    if (dim(el)[2] != 2)
      stop("Two columns expected in edge list creation.")
    if (length(nodeSet) == 1)
      stop("Three node sets need to be specified for edge lists: nodes / nodes / edges")
    l <- list(data = el,
              nodeSet = nodeSet,
              size = dim(el))
    class(l) <- "edgelist.netdist"
    l
}


# createEffectsObject
#' Title
#'
#' @param effectInit
#' @param checkProcessState
#'
#' @return
#' @export
#'
#' @examples
createEffectsObject <-
  function(effectInit, checkProcessState = NULL) {
    # TODO add default parameters to effect names
    effectNames <-
      lapply(effectInit, function(x)
        ifelse(is.list(x), do.call(paste, x), x))

    effects <-
      lapply(effectInit, function(x)
        eval(parse(text = x[[1]])))

    # Update signatures of the effects based on default parameters and above specified parameters
    signatures <- lapply(effects, formals)
    setParams <- lapply(effectInit, function(x)
      x[-1])
    signatures <- mapply(function(s, p) {
      s[names(p)] <- p
      s["cache"] <- alist(cache = NULL)
      s
    } ,
    signatures, setParams, SIMPLIFY = F)
    if ("matrix" %in% class(signatures))
      signatures <- apply(signatures, 2, invisible)

    # Assign signatures with default values to generic functions
    for (i in 1:length(effects))
      formals(effects[[i]]) <- unlist(signatures[i], recursive = F)

    # If a process state is provided, check whether all params refer to existing objects
    if (!is.null(checkProcessState)) {
      stateNames <- names(checkProcessState)
      params <- unlist(signatures)
      refs <- unique(params[names(params) == "attribute.index"])
      for (r in refs)
        if (!(r %in% stateNames))
          stop(paste("Unknown process state reference:", r))
    }

    effectslist <- list(effectFormulas = effects,
                        name = unlist(effectNames))
    class(effectslist) <- "effectsList.tssm"
    effectslist
}


# createNetwork
#' Title
#'
#' @param m
#' @param isSymmetric
#' @param isBipartite
#' @param nodeSet
#'
#' @return
#' @export
#'
#' @examples
createNetwork <-
  function(m,
           isSymmetric = FALSE,
           isBipartite = FALSE,
           nodeSet = c("actors", "actors")) {
    if (!is.matrix(m))
      stop("Not a matrix.")
    if (length(nodeSet) == 1)
      nodeSet <- c(nodeSet, nodeSet)
    l <- list(
      data = m,
      isSymmetric = isSymmetric,
      isBipartite = isBipartite,
      nodeSet = nodeSet,
      size = dim(m)
    )
    class(l) <- "network.tssm"
    l
  }


# createNodeSet
#' Title
#'
#' @param x
#' @param isPresent
#' @param considerWhenSampling
#'
#' @return
#' @export
#'
#' @examples
createNodeSet <-
  function(x = NULL,
           isPresent = NULL,
           considerWhenSampling = NULL) {
    if (is.null(x) && is.null(isPresent))
      stop("Only null parameters.")
    if (is.null(x))
      x <- length(isPresent)
    if (length(x) == 1 &&
        !is.numeric(x))
      stop("non-numeric value of single valued x.")
    ids <- if (length(x) == 1)
      1:x
    else
      x
    if (is.null(isPresent))
      isPresent <- rep(T, length(ids))
    if (is.null(considerWhenSampling))
      considerWhenSampling <- rep(T, length(ids))
    if (length(isPresent) != length(ids))
      stop("Wrong length of presence vector.")
    l <- list(
      isPresent = isPresent,
      ids = ids,
      considerWhenSampling = considerWhenSampling
    )
    class(l) <- "nodeSet.tssm"
    l
}


# createNodeVariable
#' Title
#'
#' @param values
#' @param range
#' @param nodeSet
#' @param addSame
#' @param addSim
#'
#' @return
#' @export
#'
#' @examples
createNodeVariable <-
  function(values,
           range = NULL,
           nodeSet = "actors",
           addSame = F,
           addSim = F) {
    if (!is.numeric(values) &&
        !all(is.na(values)))
      stop("Values not numeric.")
    l <- list(
      data = values,
      range = range,
      nodeSet = nodeSet,
      size = length(values)
    )

    if (addSame)
      l[["same"]] <- outer(values, values, "==") * 1
    if (addSim)
      l[["sim"]] <-
      1 - abs(outer(values, values, "-")) / (range(values)[2] - range(values)[1])

    class(l) <- "nodeVar.tssm"
    l
}


# createProcessState
#' Title
#'
#' @param elements
#'
#' @return
#' @export
#'
#' @examples
createProcessState <- function(elements) {
  if (!is.list(elements))
    stop("Expecting a list.")
  if (length(elements) == 0)
    stop("List must have at least one element.")
  if (is.null(names(elements)))
    stop("List must be named.")

  nodeSetIDs <- c()
  linkedElementIDs <- c()
  sizes <- c()
  for (i in 1:length(elements)) {
    e <- elements[[i]]
    if (!(
      class(e) %in% c(
        "edgelist.netdist",
        "nodeSet.tssm",
        "nodeVar.tssm",
        "network.tssm"
      )
    ))
      stop(paste("Unknown element of class", class(e)))

    # TODO CLEANUP from here. What is necessary, hat should be extended?

    if (class(e) == "nodeSet.tssm")
      nodeSetIDs <- c(nodeSetIDs, i)
    if (class(e) %in% c("network.tssm", "nodeVar.tssm")) {
      linkedElementIDs <- c(linkedElementIDs, i)
      sizes <- c(sizes, e$size)
    }
  }

  # if no node sets were found, create a default
  # TODO: the node set check should be based on the nodeSet names, not their size
  if (length(nodeSetIDs) == 0) {
    if (length(unique(sizes)) != 1)
      stop("Differing element sizes without defining node sets.")
    elements <-
      append(elements, list(actors = createNodeSet(sizes[1])))
  }

  # check if all elements in 'linkedElementIDs' have a corresponding node set
  # TODO implement

  class(elements) <- "processState.tssm"
  elements
}


# createWeightedCache
#' Title
#'
#' @param processState 
#' @param cacheObjectNames 
#' @param ressourceCovariates 
#'
#' @return
#' @export
#'
#' @examples
createWeightedCache <-
  function(processState,
           cacheObjectNames,
           ressourceCovariates = NULL) {
    cache <- list()

    for (name in cacheObjectNames) {
      if (!(class(processState[[name]]) %in% c("network.tssm", "edgelist.netdist")))
        stop(paste(name, "is not a network or edgelist."))

      nodeSet1 <- processState[[processState[[name]]$nodeSet[1]]]$ids
      nodeSet2 <- processState[[processState[[name]]$nodeSet[2]]]$ids
      nActors1 <- length(nodeSet1)
      nActors2 <- length(nodeSet2)

      cache[[name]] <- list()
      if (class(processState[[name]]) == "network.tssm")
        cache[[name]]$valuedNetwork <-  processState[[name]]$data
      if (class(processState[[name]]) == "edgelist.netdist") {
        # create valued network from edge list
        m <- matrix(0, nActors1, nActors2)

        # create weighted ressource networks
        m.ressource <-
          lapply(ressourceCovariates, function(v)
            matrix(0, nrow = nActors1, ncol = nActors2))
        names(m.ressource) <- ressourceCovariates

        for (i in 1:processState[[name]]$size[1]) {
          sender <- processState[[name]]$data[i, 1]
          receiver <- processState[[name]]$data[i, 2]
          v <- m[sender, receiver]
          m[sender, receiver] <- v + 1

          # cache network * ressource covariate matrices
          for (ressCovar in ressourceCovariates) {
            v <- m.ressource[[ressCovar]] [sender, receiver]
            m.ressource[[ressCovar]][sender, receiver] <- v +
              processState[[ressCovar]]$data[i]
          }

        }
        cache[[name]]$valuedNetwork <- m
        cache[[name]]$ressourceNetworks <- m.ressource

        # edge ids
        edgeSet <- processState[[processState[[name]]$nodeSet[3]]]$ids
        nEdges <- length(edgeSet)
      }
      if (nActors1 == nActors2) {
        cache[[name]]$netFlowsNetwork <-
          cache[[name]]$valuedNetwork - t(cache[[name]]$valuedNetwork)
        cache[[name]]$minNetwork <-
          matrix(mapply(min, cache[[name]]$valuedNetwork, t(cache[[name]]$valuedNetwork)) ,
                 nActors1,
                 nActors2)
        diag(cache[[name]]$minNetwork) <- 0
      }

    } # end for loop

    cache
}


# estimateMobilityNetwork
#' Title
#'
#' @aliases estimateDistributionNetwork
#' @param dep.var
#' @param state
#' @param cache
#' @param effects
#' @param initialParameters
#' @param prevAns
#' @param burnInN1
#' @param iterationsN1
#' @param thinningN1
#' @param gainN1
#' @param burnInN2
#' @param nsubN2
#' @param initGain
#' @param thinningN2
#' @param initialIterationsN2
#' @param iterationsN3
#' @param burnInN3
#' @param thinningN3
#' @param allowLoops
#' @param parallel
#' @param cpus
#' @param verbose
#' @param returnDeps
#' @param multinomialProposal
#' @param fish
#'
#' @return
#' @export
#'
#' @examples
estimateMobilityNetwork <-
  function(dep.var,
           state,
           cache,
           effects,
           initialParameters = NULL,
           prevAns = NULL,
           burnInN1 = 500,
           iterationsN1 = 200,
           thinningN1 = 20,
           gainN1 = 1,
           burnInN2 = burnInN1,
           nsubN2 = 6,
           initGain = 0.02,
           thinningN2 = 100,
           initialIterationsN2 = 50,
           iterationsN3 = 1000,
           burnInN3 = 2000,
           thinningN3 = 30,
           allowLoops = NULL,
           parallel = F,
           cpus = 3,
           verbose = F,
           returnDeps = F,
           multinomialProposal = F,
           fish = F) {
    # set parameters to default values if not defined explicitly
    if (is.null(initialParameters))
      initialParameters <- rep(0, length(effects$name))
    if (is.null(allowLoops)) {
      allowLoops <- T
      if (all(diag(cache[[dep.var]]$valuedNetwork) == 0))
        allowLoops <- F
    }

    observedStatistics <-
      getNetworkStatistics(dep.var, state, cache, effects)

    # check wether the prevAns object has the same effects, if not cause error
    if (!is.null(prevAns)) {
      if (!(length(effects$name) == length(names(prevAns$estimates)))) {
        stop("prevAns object has different number of effects as the effects object")
      }
      if (!all(effects$name == names(prevAns$estimates))) {
        stop("effect names of prevAns object and effects object do not match")
      }
    }

    # decide whether to run phase 1 based on whether a prevAns object is specified
    if (!is.null(prevAns)) {
      # skip phase 1 as all information is already there
      cat("skipping phase 1 and taking values from prevAns\n")
      initialParameters <- prevAns$estimates
      sensitivityVector <- 1 / diag(solve(prevAns$covarianceMatrix))
    } else {
      # run phase 1 to get initial estimates and a sensitivity vector
      if (verbose)
        cat("Starting phase 1\n")
      resPhase1 <-
        runPhase1(
          dep.var,
          state,
          cache,
          effects,
          initialParameters,
          burnInN1,
          iterationsN1,
          thinningN1,
          gainN1,
          allowLoops,
          multinomialProposal = multinomialProposal,
          verbose
        )
      initialParameters <- resPhase1$estimates
      sensitivityVector <- resPhase1$sensitivityVector
    }

    # run phase 2 to get parameter estimates
    if (verbose)
      cat("Starting phase 2\n")
    resPhase2 <- runPhase2(
      dep.var,
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
      parallel = parallel,
      cpus = cpus,
      multinomialProposal = multinomialProposal,
      allowLoops,
      verbose = verbose
    )

    # run phase 3 to get final covariance matrix and convergence check
    if (verbose)
      cat("Starting phase 3\n")
    resPhase3 <- runPhase3(
      dep.var,
      state,
      cache,
      effects,
      resPhase2,
      observedStatistics,
      iterationsN3,
      burnInN3,
      thinningN3,
      parallel = parallel,
      cpus = cpus,
      allowLoops,
      verbose = verbose,
      returnDeps = returnDeps,
      multinomialProposal = multinomialProposal,
      fish = fish
    )

    if (!returnDeps) {
      result <- list(
        estimates = resPhase2,
        standardErrors = sqrt(diag(resPhase3$covarianceMatrix)),
        covarianceMatrix = resPhase3$covarianceMatrix,
        convergenceStatistics = resPhase3$convergenceStatistics,
        state = state,
        cache = cache
      )
    }

    if (returnDeps) {
      result <- list(
        estimates = resPhase2,
        standardErrors = sqrt(diag(resPhase3$covarianceMatrix)),
        covarianceMatrix = resPhase3$covarianceMatrix,
        convergenceStatistics = resPhase3$convergenceStatistics,
        state = state,
        cache = cache,
        deps = resPhase3$deps
      )
    }

    class(result) <- "result.netdist"

    return(result)

}


# estimateDistributionNetwork
estimateDistributionNetwork <- estimateMobilityNetwork


# simulateDistributionNetworks
#' Title
#'
#' @param cache
#' @param state
#' @param effects
#' @param ans
#' @param allowLoops
#' @param burnin
#' @param thinning
#' @param nSimulations
#' @param dep.var
#'
#' @return
#' @export
#'
#' @examples
simulateDistributionNetworks <-
  function(cache,
           state,
           effects,
           ans,
           allowLoops,
           burnin,
           thinning,
           nSimulations,
           dep.var) {
    # extract parameters from ans
    parameters <- ans$estimates

    # generate a state and cache after burnin with parameters
    r <-
      simulateNSteps(
        dep.var = dep.var,
        state = state,
        cache = cache,
        effects = effects,
        parameters = parameters,
        allowLoops = allowLoops,
        n = burnin
      )

    # extract state amd cache
    simState <- r$state
    simCache <- r$cache

    # generate list in which simulated states and caches will be stored
    simulatedList <- list()

    for (i in 1:nSimulations) {
      r <-
        simulateNSteps(
          dependentVariable,
          simState,
          simCache,
          effects,
          allowLoops = F,
          n = thinning,
          parameters = parameters
        )
      simState <- r$state
      simCache <- r$cache
      simulatedList[[i]] <- list(state = simState, cache = simCache)
      int <- i %% 11
      s <- "_,.-'``'-.,"
      cat(substr(s, int + 1, int + 1))
      if (runif(1) < 0.02)
        cat("<o_><")
    }
    class(simulatedList) <- "sims.NetDist"
    return(simulatedList)
}
