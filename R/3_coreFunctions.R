####### coreFunctions


#' createAlgorithm
#' 
#' Specifies the algorithm used in later model estimation based on characteristics 
#' of the state and the effects.
#'
#' @param state
#' @param effects
#' @param burnInN1 The number of simulation steps before the first draw in Phase 1.
#' A recommended value is at least n_Individuals * n_locations if
#' multinomialProposal = F, and at least n_Individuals if multinomialProposal = T 
#' which is set as default.
#' @param iterationsN1 The number of draws taken in Phase 1.
#' A recommended value is at least 4 * n_effects which is set as default.
#' If the value is too low, there will be an error in Phase 1.
#' @param thinningN1 The number of simulation steps between two draws in Phase 1.
#' A recommended value is at least 0.5 * n_Individuals * n_locations if
#' multinomialProposal = F, and at least n_Individuals if multinomialProposal = T
#' which is set as default.
#' @param gainN1 The size of the updating step after Phase 1. A conservative
#' value is 0, values higher than 0.25 are courageous. Defaults to 0.1.
#' @param burnInN2 The number of simulation steps before the first draw in Phase 1.
#' A recommended value is at least n_Individuals * n_locations if
#' multinomialProposal = F, and at least n_Individuals if multinomialProposal = T 
#' which is set as default.
#' @param nsubN2 Number of subphases in Phase 2. In case this is the first
#' estimation, 4 subphases are recommended and is set as default. If convergence 
#' in a previous estimation was close, then 1-2 subphases should be enough.
#' @param initGain The magnitude of parameter updates in the first subphase of
#' Phase 2. Values of around 0.2 (default) are recommended.
#' @param thinningN2 The number of simulation steps between two draws in Phase 2.
#' A recommended value is at least 0.5 * n_Individuals * n_locations if
#' multinomialProposal = F, and at least n_Individuals if multinomialProposal = T
#' which is set as default.
#' @param initialIterationsN2 The number of draws taken in subpphase 1 of Phase 2.
#' For first estimations, a recommended value is around 20-50 (default to 25).
#' Note that in later subphases, the number of iterations increases.
#' If this is a further estimation to improve convergence, higher values (100+)
#' are recommended.
#' @param iterationsN3 Number of draws in Phase 3. Recommended are at the very 
#' least 500 (default).
#' In case this value is too low, the outcome might erroneously indicate a lack
#' of convergence.
#' @param burnInN3 The number of simulation steps before the first draw in Phase 3.
#' A recommended value is at least 3 * n_Individuals * n_locations if
#' multinomialProposal = F, and at least 3 * n_Individuals if multinomialProposal = T
#' which is set as default.
#' @param thinningN3 The number of simulation steps between two draws in Phase 3.
#' A recommended value is at least n_Individuals * n_locations if
#' multinomialProposal = F, and at least 2 * n_Individuals if multinomialProposal = T
#' which is set as default.
#' In case this value is too low, the outcome might erroneously indicate a lack
#' of convergence.
#' @param allowLoops Logical: can individuals/resources stay in their origin?
#'
#' @return
#' @export
#'
#' @seealso [createProcessState()], [createEffectsObject()], [estimateMobilityNetwork()]
#'
#' @examples
#' # define algorithm based on state and effects characteristics
#' myAlg <- createAlgorithm(myState, myEffects, myCache, myDependentVariable)
createAlgorithm <-
  function(state,
           effects,
           cache,
           dep.var,
           burnInN1 = NULL,
           iterationsN1 = NULL,
           thinningN1 = NULL,
           gainN1 = 0.1,
           burnInN2 = NULL,
           nsubN2 = 4,
           initGain = 0.2,
           thinningN2 = NULL,
           initialIterationsN2 = 25,
           iterationsN3 = 500,
           burnInN3 = NULL,
           thinningN3 = NULL,
           allowLoops = NULL) {
    algorithm <- list()

    if (is.null(burnInN1)) {
      algorithm[["burnInN1"]] <- as.numeric(c(length(state[["people"]][["ids"]]) * length(state[["locations"]][["ids"]]), 
                                              length(state[["people"]][["ids"]])))
    } else {
      algorithm[["burnInN1"]] <- as.numeric(burnInN1)
    }

    if (is.null(iterationsN1)) {
      algorithm[["iterationsN1"]] <- as.numeric(length(effects[["effectFormulas"]]) * 4)
    } else {
      algorithm[["iterationsN1"]] <- as.numeric(iterationsN1)
    }

    if (is.null(thinningN1)) {
      algorithm[["thinningN1"]] <- as.numeric(c(length(state[["people"]][["ids"]]) * length(state[["locations"]][["ids"]]) * 0.5,
                                              length(state[["people"]][["ids"]])))
    } else {
      algorithm[["thinningN1"]] <- as.numeric(thinningN1)
    }

    algorithm[["gainN1"]] <- as.numeric(gainN1)

    if (is.null(burnInN2)) {
      algorithm[["burnInN2"]] <- as.numeric(c(length(state[["people"]][["ids"]]) * length(state[["locations"]][["ids"]]), 
                                              length(state[["people"]][["ids"]])))
    } else {
      algorithm[["burnInN2"]] <- as.numeric(burnInN2)
    }

    algorithm[["nsubN2"]] <- as.numeric(nsubN2)

    algorithm[["initGain"]] <- as.numeric(initGain)

    if (is.null(thinningN2)) {
      algorithm[["thinningN2"]] <- as.numeric(c(length(state[["people"]][["ids"]]) * length(state[["locations"]][["ids"]]) * 0.5,
                                              length(state[["people"]][["ids"]])))
    } else {
      algorithm[["thinningN2"]] <- as.numeric(thinningN2)
    }

    algorithm[["initialIterationsN2"]] <- as.numeric(initialIterationsN2)

    algorithm[["iterationsN3"]] <- as.numeric(iterationsN3)

    if (is.null(burnInN3)) {
      algorithm[["burnInN3"]] <- as.numeric(c(length(state[["people"]][["ids"]]) * length(state[["locations"]][["ids"]]) * 3,
                                            length(state[["people"]][["ids"]]) * 3))
    } else {
      algorithm[["burnInN3"]] <- as.numeric(burnInN3)
    }

    if (is.null(thinningN3)) {
      algorithm[["thinningN3"]] <- as.numeric(c(length(state[["people"]][["ids"]]) * length(state[["locations"]][["ids"]]),
                                              length(state[["people"]][["ids"]]) * 2))
    } else {
      algorithm[["thinningN3"]] <- as.numeric(thinningN3)
    }

    if (is.null(allowLoops)) {
      algorithm[["allowLoops"]] <- TRUE
      if (all(diag(cache[[dep.var]]$valuedNetwork) == 0)) {
        algorithm[["allowLoops"]] <- FALSE
      }
    } else {
      algorithm[["allowLoops"]] <- as.logical(allowLoops)
    }

    class(algorithm) <- "algorithm.monan"
    algorithm
  }


#' createEdgelist
#'
#' Creates an edgelist object, which is the standard format of the outcome to be modelled
#' by monan.
#'
#' @param el An edgelist in the form of a matrix with two columns and N rows.
#' The first column indicates the origin of a person/resource, the second row the destination.
#' Each row is one observation.
#' @param nodeSet The nodesets of the edgelists. For edgelists, this is a vector with three
#' entries describing between what locations who is mobile, i.e., c(location, location, individuals/resources).
#'
#' @return an object of type edgelist.monan
#' @export
#'
#' @seealso [createProcessState()]
#'
#' @examples
#' # create an object of class edgelist.monan
#' transfers <- createEdgelist(mobilityEdgelist, c("organisations", "organisations", "people"))
createEdgelist <-
  function(el, nodeSet = c("actors", "actors", "edges")) {
    if (dim(el)[2] != 2) {
      stop("Two columns expected in edge list creation.")
    }
    if (length(nodeSet) == 1) {
      stop("Three node sets need to be specified for edge lists: nodes / nodes / edges")
    }
    l <- list(
      data = el,
      nodeSet = nodeSet,
      size = dim(el)
    )
    class(l) <- "edgelist.monan"
    l
  }


#' createEffectsObject
#'
#' Specifies the model with endogenous and exogenous predictors.
#' The predictors in the model are called “effects”.
#'
#' @param effectInit A list of "effects", where each effect to be included is specified as a
#' further list that contains the effect name and the additional parameters it needs.
#' Effects without further parameters only contain the effect name (e.g., loops).
#' @param checkProcessState For internal use only.
#'
#' @return A monan effects object.
#' @export
#'
#' @examples
#' # create an effects object
#' myEffects <- createEffectsObject(
#'   list(
#'     list("loops"),
#'     list("min_reciprocity"),
#'     list("dyadic_covariate", attribute.index = "sameRegion"),
#'     list("alter_covariate", attribute.index = "size"),
#'     list("resource_covar_to_node_covar",
#'       attribute.index = "region",
#'       resource.attribute.index = "sex"
#'     ),
#'     list("loops_resource_covar", resource.attribute.index = "sex")
#'   )
#' )
createEffectsObject <-
  function(effectInit, checkProcessState = NULL) {
    # TODO add default parameters to effect names
    effectNames <-
      lapply(effectInit, function(x) {
        ifelse(is.list(x), do.call(paste, x), x)
      })

    effects <-
      lapply(effectInit, function(x) {
        eval(parse(text = x[[1]]))
      })

    # Update signatures of the effects based on default parameters and above specified parameters
    signatures <- lapply(effects, formals)
    setParams <- lapply(effectInit, function(x) {
      x[-1]
    })
    signatures <- mapply(
      function(s, p) {
        s[names(p)] <- p
        s["cache"] <- alist(cache = NULL)
        s
      },
      signatures, setParams,
      SIMPLIFY = FALSE
    )
    if ("matrix" %in% class(signatures)) {
      signatures <- apply(signatures, 2, invisible)
    }

    # Assign signatures with default values to generic functions
    for (i in 1:length(effects)) {
      formals(effects[[i]]) <- unlist(signatures[i], recursive = FALSE)
    }

    # If a process state is provided, check whether all params refer to existing objects
    if (!is.null(checkProcessState)) {
      stateNames <- names(checkProcessState)
      params <- unlist(signatures)
      refs <- unique(params[names(params) == "attribute.index"])
      for (r in refs) {
        if (!(r %in% stateNames)) {
          stop(paste("Unknown process state reference:", r))
        }
      }
    }

    effectslist <- list(
      effectFormulas = effects,
      name = unlist(effectNames)
    )
    class(effectslist) <- "effectsList.monan"
    effectslist
  }


#' createNetwork
#'
#' Defines a network between locations, generally to be used as a predictor in the model.
#' NOTE: The outcome variable of the model is not defined as a network, but as an edgelist!
#'
#' @param m A square matrix containing the network data.
#' @param isSymmetric Currently not in use, defaults to FALSE.
#' @param isBipartite Currently not in use, defaults to FALSE.
#' @param nodeSet Which nodeset are the nodes of the network. Usually this will
#' be the locations in the data.
#'
#' @return a monan network object.
#' @export
#'
#' @seealso [createProcessState()], [createEdgelist()]
#'
#' @examples
#' # create an object of class network.monan
#' sameRegion <- outer(orgRegion, orgRegion, "==") * 1
#' sameRegion <- createNetwork(sameRegion, nodeSet = c("locations", "locations"))
createNetwork <-
  function(m,
           isSymmetric = FALSE,
           isBipartite = FALSE,
           nodeSet = c("actors", "actors")) {
    if (!is.matrix(m)) {
      stop("Not a matrix.")
    }
    if (length(nodeSet) == 1) {
      nodeSet <- c(nodeSet, nodeSet)
    }
    l <- list(
      data = m,
      isSymmetric = isSymmetric,
      isBipartite = isBipartite,
      nodeSet = nodeSet,
      size = dim(m)
    )
    class(l) <- "network.monan"
    l
  }


#' createNodeSet
#'
#' Determines and names the nodesets of individuals and locations that make up the mobility network.
#'
#' @param x Either a single number indicating how many items are in this node-set
#' or a vector from 1:n_items.
#' @param isPresent Currently not in use.
#' @param considerWhenSampling A boolean/logical vector of the length of the nodeset.
#' Only in use in special cases.
#' If the nodeset indicates a location, considerWhenSampling indicates whether the
#' location is a possible destination, or is only an origin (e.g. a training facility).
#' Entries in the vector of locations that cannot be a destination are FALSE.
#' If the nodeset indicates mobile individuals, considerWhenSampling indicates whether
#' their mobility should be modelled or whether it is structurally determined, that
#' is, their mobility is exogenously defined and does not follow the same logic as
#' the mobility of everybody else.
#'
#' @return A monan nodeset object.
#' @export
#'
#' @seealso [createProcessState()]
#'
#' @examples
#' # create an object of class nodeSet.monan
#' people <- createNodeSet(1:nrow(mobilityEdgelist))
#' locations <- createNodeSet(length(orgRegion))
createNodeSet <-
  function(x = NULL,
           isPresent = NULL,
           considerWhenSampling = NULL) {
    if (is.null(x) && is.null(isPresent)) {
      stop("Only null parameters.")
    }
    if (is.null(x)) {
      x <- length(isPresent)
    }
    if (length(x) == 1 &&
      !is.numeric(x)) {
      stop("non-numeric value of single valued x.")
    }
    ids <- if (length(x) == 1) {
      1:x
    } else {
      x
    }
    if (is.null(isPresent)) {
      isPresent <- rep(T, length(ids))
    }
    if (is.null(considerWhenSampling)) {
      considerWhenSampling <- rep(T, length(ids))
    }
    if (length(isPresent) != length(ids)) {
      stop("Wrong length of presence vector.")
    }
    l <- list(
      isPresent = isPresent,
      ids = ids,
      considerWhenSampling = considerWhenSampling
    )
    class(l) <- "nodeSet.monan"
    l
  }


#' createNodeVariable
#'
#' Assigns an covariate to one nodeset, i.e., an exogenous characteristic of mobile
#' individuals/resources or locations.
#'
#' @param values A vector assigning the covariate value to each element of the nodeset.
#' @param range The range of possible values, currently not in use.
#' @param nodeSet The nodeset to which the covariate applies.
#' @param addSame Will the variable be used to model categorical homophily (e.g.,
#' with the same_covariate effect)? In this case, addSame needs to be set to TRUE.
#' @param addSim Will the variable be used to model continuous homophily (e.g.,
#' with the sim_covariate effect)? In this case, addSame needs to be set to TRUE.
#'
#' @return A monan node variable object
#' @export
#'
#' @seealso [createProcessState()]
#'
#' @examples
#' # create an object of class nodeVar.monan
#' region <- createNodeVariable(orgRegion, nodeSet = "locations")
#' size <- createNodeVariable(orgSize, nodeSet = "locations", addSim = TRUE)
#' sex <- createNodeVariable(indSex, nodeSet = "people")
createNodeVariable <-
  function(values,
           range = NULL,
           nodeSet = "actors",
           addSame = FALSE,
           addSim = FALSE) {
    if (!is.numeric(values) &&
      !all(is.na(values))) {
      stop("Values not numeric.")
    }
    l <- list(
      data = values,
      range = range,
      nodeSet = nodeSet,
      size = length(values)
    )

    if (addSame) {
      l[["same"]] <- outer(values, values, "==") * 1
    }
    if (addSim) {
      l[["sim"]] <-
        1 - abs(outer(values, values, "-")) / (range(values)[2] - range(values)[1])
    }

    class(l) <- "nodeVar.monan"
    l
  }


#' createProcessState
#'
#' Creates the "Process state", i.e., a MoNAn object that stores all information
#' about the data that will be used in the estimation. This includes the
#' outcome variable (edgelist), the nodesets, and all covariates.
#'
#' @param elements A named list of the outcome variable (edgelist), the nodesets,
#' and all covariates that contain the information
#' about the data that will be used in the estimation.
#'
#' @return A monan process state object.
#' @export
#'
#' @seealso [createEdgelist()], [createNodeSet()],
#' [createNodeVariable()], [createNetwork()]
#'
#' @examples
#' # Create a process state out of the mobility data objects:
#' # create objects (which are later combined to the process state)
#' transfers <- createEdgelist(mobilityEdgelist,
#'   nodeSet = c("organisations", "organisations", "people")
#' )
#' people <- createNodeSet(1:nrow(mobilityEdgelist))
#' organisations <- createNodeSet(1:length(orgRegion))
#' sameRegion <- outer(orgRegion, orgRegion, "==") * 1
#' sameRegion <- createNetwork(sameRegion,
#'   nodeSet = c("organisations", "organisations")
#' )
#' region <- createNodeVariable(orgRegion, nodeSet = "organisations")
#' size <- createNodeVariable(orgSize, nodeSet = "organisations", addSim = TRUE)
#' sex <- createNodeVariable(indSex, nodeSet = "people")
#'
#' # combine created objects to the process state
#' myState <- createProcessState(list(
#'   transfers = transfers,
#'   people = people,
#'   organisations = organisations,
#'   sameRegion = sameRegion,
#'   region = region,
#'   size = size,
#'   sex = sex
#' ))
createProcessState <- function(elements) {
  if (!is.list(elements)) {
    stop("Expecting a list.")
  }
  if (length(elements) == 0) {
    stop("List must have at least one element.")
  }
  if (is.null(names(elements))) {
    stop("List must be named.")
  }

  nodeSetIDs <- c()
  linkedElementIDs <- c()
  sizes <- c()
  for (i in 1:length(elements)) {
    e <- elements[[i]]
    if (!(
      class(e) %in% c(
        "edgelist.monan",
        "nodeSet.monan",
        "nodeVar.monan",
        "network.monan"
      )
    )) {
      stop(paste("Unknown element of class", class(e)))
    }

    # TODO CLEANUP from here. What is necessary, hat should be extended?

    if (is(e, "nodeSet.monan")) {
      nodeSetIDs <- c(nodeSetIDs, i)
    }
    if (class(e) %in% c("network.monan", "nodeVar.monan")) {
      linkedElementIDs <- c(linkedElementIDs, i)
      sizes <- c(sizes, e$size)
    }
  }

  # if no node sets were found, create a default
  # TODO: the node set check should be based on the nodeSet names, not their size
  if (length(nodeSetIDs) == 0) {
    if (length(unique(sizes)) != 1) {
      stop("Differing element sizes without defining node sets.")
    }
    elements <-
      append(elements, list(actors = createNodeSet(sizes[1])))
  }

  # check if all elements in 'linkedElementIDs' have a corresponding node set
  # TODO implement

  class(elements) <- "processState.monan"
  elements
}


#' createWeightedCache
#'
#' Creates a necessary internal object used in simulating the chains in the
#' simulation and estimation of the model.
#' In case variables of the individuals in the data are included in the state,
#' they need to be explicitly mentioned in the creation of the cache under
#' “resourceCovariates”.
#'
#' @param processState The processs state that provides the data basis for
#' creating the cache.
#' @param cacheObjectNames The name of the outcome variable, i.e., the
#' edgelist in the data that should be modelled as a character vector.
#' @param resourceCovariates A vector of resource covariates that will be
#' used in the model specification?
#'
#' @return A monan cache object.
#' @export
#'
#' @seealso [createProcessState()]
#'
#' @examples
#' # define dependent variable and create cache object
#' myDependentVariable <- "transfers"
#' myCache <- createWeightedCache(myState, myDependentVariable, resourceCovariates = c("sex"))
createWeightedCache <-
  function(processState,
           cacheObjectNames,
           resourceCovariates = NULL) {
    cache <- list()

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


#' estimateMobilityNetwork
#'
#' The core function of the package in which the model for the analysis of
#' mobility tables is estimated.
#'
#' @param dep.var The outcome variable that is modelled.
#' @param state A monan state object that contains all relevant information about
#' the outcome in the form of an edgelist, the nodesets, and covariates.
#' @param cache A monan cache object created from the same state object that is
#' used in the estimation
#' @param effects An effect object that specifies the model.
#' @param initialParameters Starting values for the parameters. Using starting
#' values e.g. from a multinomial logit model (see getMultinomialStatistics())
#' increases the chances of model convergence in the first run of the estimation
#' considerably.
#' @param prevAns If a previous estimation did not yield satisfactory convergence,
#' the outcome object of that estimation should be specified here to provide new
#' starting values for the estimation.
#' @param burnInN1 The number of simulation steps before the first draw in Phase 1.
#' A recommended value is at least n_Individuals * n_locations if
#' multinomialProposal = F, and at least n_Individuals if multinomialProposal = T
#' @param iterationsN1 The number of draws taken in Phase 1.
#' A recommended value is at least 4 * n_effects.
#' If the value is too low, there will be an error in Phase 1.
#' @param thinningN1 The number of simulation steps between two draws in Phase 1.
#' A recommended value is at least 0.5 * n_Individuals * n_locations if
#' multinomialProposal = F, and at least n_Individuals if multinomialProposal = T.
#' @param gainN1 The size of the updating step after Phase 1. A conservative
#' value is 0, values higher than 0.25 are courageous.
#' @param burnInN2 The number of simulation steps before the first draw in Phase 2.
#' A recommended value is at least n_Individuals * n_locations if
#' multinomialProposal = F, and at least n_Individuals if multinomialProposal = T
#' @param nsubN2 Number of subphases in Phase 2. In case this is the first
#' estimation, 4 subphases are recommended. If convergence in a previous estimation
#' was close, then 1-2 subphases should be enough.
#' @param initGain The magnitude of parameter updates in the first subphase of
#' Phase 2. Values of around 0.2 are recommended.
#' @param thinningN2 The number of simulation steps between two draws in Phase 2.
#' A recommended value is at least 0.5 * n_Individuals * n_locations if
#' multinomialProposal = F, and at least n_Individuals if multinomialProposal = T.
#' @param initialIterationsN2 The number of draws taken in subpphase 1 of Phase 2.
#' For first estimations, a recommended value is around 20-50.
#' Note that in later subphases, the number
#' of iterations increases.
#' If this is a further estimation to improve convergence, higher values (100+)
#' are recommended.
#' If the value is too low, the model might take more runs to converge.
#' @param iterationsN3 Number of draws in Phase 3.
#' Recommended are at the very least 500.
#' In case this value is too low, the outcome might erroneously indicate a lack
#' of convergence.
#' @param burnInN3 The number of simulation steps before the first draw in Phase 3.
#' A recommended value is at least 3 * n_Individuals * n_locations if
#' multinomialProposal = F, and at least 3 * n_Individuals if multinomialProposal = T
#' @param thinningN3 The number of simulation steps between two draws in Phase 3.
#' A recommended value is at least n_Individuals * n_locations if
#' multinomialProposal = F, and at least 2 * n_Individuals if multinomialProposal = T.
#' In case this value is too low, the outcome might erroneously indicate a lack
#' of convergence.
#' @param allowLoops Logical: can individuals/resources stay in their origin?
#' @param parallel Computation on multiple cores?
#' @param cpus Number of cores for computation in case parallel = T.
#' @param verbose Logical: display information about estimation progress in the console?
#' @param returnDeps Logical: should the simulated values of Phase 3 be stored and returned?
#' This is necessary to run GoF tests.
#' Note that this might result in very large objects.
#' @param multinomialProposal How should the next possible outcome in the simulation chains
#' be sampled? If TRUE, fewer simulation steps are needed, but each simulation
#' step takes considerably longer.
#' @param fish Logical: display a fish?
#'
#' @aliases estimateDistributionNetwork
#'
#' @return A monan results object that contains the estimates, standard errors,
#' and convergence statistics. Furthermore, the covariance matrix used to calculate
#' the standard errors is included, which also shows colinearity between effects.
#' In case returnDeps = TRUE, the simulations of Phase 3 are included, too.
#' @export
#'
#' @seealso [createProcessState()], [createWeightedCache()],
#' [createEffectsObject()]
#'
#' @examples
#' \dontrun{
#' # estimate mobility network
#' myResDN <- estimateMobilityNetwork(myDependentVariable,
#'   myState, myCache, myEffects, myAlg,
#'   initialParameters = NULL,
#'   parallel = TRUE, cpus = 4,
#'   verbose = TRUE,
#'   returnDeps = TRUE,
#'   multinomialProposal = FALSE,
#'   fish = FALSE
#' )
#'
#' myResDN
#'
#' # estimate mobility network again based on previous results to improve convergence
#' myResDN <- estimateMobilityNetwork(myDependentVariable,
#'   myState, myCache, myEffects, myAlg,
#'   prevAns = myResDN,
#'   parallel = TRUE, cpus = 4,
#'   verbose = TRUE,
#'   returnDeps = TRUE,
#'   multinomialProposal = FALSE,
#'   fish = FALSE
#' )
#' }
estimateMobilityNetwork <-
  function(dep.var,
           state,
           cache,
           effects,
           algorithm,
           initialParameters = NULL,
           prevAns = NULL,
           parallel = FALSE,
           cpus = 3,
           verbose = FALSE,
           returnDeps = FALSE,
           multinomialProposal = FALSE,
           fish = FALSE) {
    # set parameters to default values if not defined explicitly
    if (is.null(initialParameters)) {
      initialParameters <- rep(0, length(effects$name))
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
      if (verbose) {
        cat("Starting phase 1\n")
      }
      resPhase1 <-
        runPhase1(
          dep.var,
          state,
          cache,
          effects,
          initialParameters,
          if (multinomialProposal == FALSE) {
            algorithm$burnInN1[1]} else {algorithm$burnInN1[2]},
          algorithm$iterationsN1,
          if (multinomialProposal == FALSE) {
            algorithm$thinningN1[1]} else {algorithm$thinningN1[2]},
          algorithm$gainN1,
          algorithm$allowLoops,
          multinomialProposal = multinomialProposal,
          verbose
        )
      initialParameters <- resPhase1$estimates
      sensitivityVector <- resPhase1$sensitivityVector
    }

    # run phase 2 to get parameter estimates
    if (verbose) {
      cat("Starting phase 2\n")
    }
    resPhase2 <- runPhase2(
      dep.var,
      state,
      cache,
      effects,
      initialParameters,
      sensitivityVector,
      if (multinomialProposal == FALSE) {
        algorithm$burnInN2[1]} else {algorithm$burnInN2[2]},
      algorithm$nsubN2,
      algorithm$initGain,
      if (multinomialProposal == FALSE) {
        algorithm$thinningN2[1]} else {algorithm$thinningN2[2]},
      algorithm$initialIterationsN2,
      parallel = parallel,
      cpus = cpus,
      multinomialProposal = multinomialProposal,
      algorithm$allowLoops,
      verbose = verbose
    )

    # run phase 3 to get final covariance matrix and convergence check
    if (verbose) {
      cat("Starting phase 3\n")
    }
    resPhase3 <- runPhase3(
      dep.var,
      state,
      cache,
      effects,
      resPhase2,
      observedStatistics,
      algorithm$iterationsN3,
      if (multinomialProposal == FALSE) {
        algorithm$burnInN3[1]} else {algorithm$burnInN3[2]},
      if (multinomialProposal == FALSE) {
        algorithm$thinningN3[1]} else {algorithm$thinningN3[2]},
      parallel = parallel,
      cpus = cpus,
      algorithm$allowLoops,
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

    class(result) <- "result.monan"

    return(result)
  }


#' estimateDistributionNetwork
#'
#' @rdname estimateMobilityNetwork
estimateDistributionNetwork <- estimateMobilityNetwork


#' simulateMobilityNetworks
#'
#' Simulates mobility networks for a given data, effects, and parameters. This
#' function is mainly interesting to explore the behaviour of the model or to
#' do counter-factual simulations.
#'
#' @aliases simulateDistributionNetworks
#' @param dep.var The name of the outcome variable that is simulated.
#' @param state A monan state object that contains all relevant information about
#' nodesets, and covariates. Further, an edgelist of the dependent variable needs
#' to be specified with the initial mobility network as starting value for the
#' simulation. For a large enough burnin, any initial mobility network is allowed.
#' @param cache A monan cache object created from the same state object that is
#' used in the simulation.
#' @param effects An effect object that specifies the model.
#' @param parameters The parameters associated with the effects that shall be used
#' in the simulations.
#' @param allowLoops Logical: can individuals/resources stay in their origin?
#' @param burnin The number of simulation steps that are taken before the first draw of a
#' network is taken. A number too small will mean the first draw is influenced
#' by the initially specified network. A recommended value for the lower bound is 3 * n_Individuals *
#' n_locations.
#' @param thinning The number of simulation steps that are taken between two draws of a
#' network. A recommended value for the lower bound is n_Individuals * n_locations.
#' @param nSimulations The number of mobility networks to be simulated.
#'
#' @return A list with nSimulations entries, where each entry contains a further list with the
#' state and the cache of the current simulation stored.
#' @export
#'
#' @examples
#' \dontrun{
#' # simulate a mobility network
#' mySimDN <- simulateMobilityNetworks(myDependentVariable,
#'   myState,
#'   myCache,
#'   myEffects,
#'   parameters = c(2, 1, 1.5, 0.1, -1, -0.5),
#'   allowLoops = TRUE,
#'   burnin = 45000,
#'   thinning = 15000,
#'   nSimulations = 10
#' )
#'
#' mySimDN[[1]]
#' }
simulateMobilityNetworks <-
  function(dep.var,
           state,
           cache,
           effects,
           parameters,
           allowLoops,
           burnin,
           thinning,
           nSimulations) {
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
          dep.var,
          simState,
          simCache,
          effects,
          allowLoops = FALSE,
          n = thinning,
          parameters = parameters
        )
      simState <- r$state
      simCache <- r$cache
      simulatedList[[i]] <- list(state = simState, cache = simCache)
      int <- i %% 11
      s <- "_,.-'``'-.,"
      cat(substr(s, int + 1, int + 1))
      if (runif(1) < 0.02) {
        cat("<o_><")
      }
    }
    class(simulatedList) <- "sims.monan"
    return(simulatedList)
  }

#' simulateDistributionNetworks
#'
#' @rdname simulateMobilityNetworks
simulateDistributionNetworks <- simulateMobilityNetworks
