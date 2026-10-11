########## C++ sampler
# The MCMC chain can run in C++ (src/sampler.cpp) for models whose effects all
# have a C++ version (src/effects.cpp). For the same seed, it gives the same
# chain as the R sampler in 1_hiddenFunctions.R.

# effects with a C++ version
cppEffects <- c("loops", "concentration_AC", "reciprocity_AC", "dyadic_covariate",
                "alter_covariate", "resource_covar_to_node_covar",
                "loops_resource_covar",
                "loops_AC", "loops_GW", "loops_additional_origin",
                "loops_x_loops_additional_origin", "loops_node_covar",
                "loops_resource_covar_node_covar",
                "concentration_GW", "concentration_AC_dyad_covar",
                "concentration_GW_dyad_covar", "concentration_AC_resource_covar_bin",
                "concentration_GW_resource_covar_bin",
                "reciprocity_GW", "reciprocity_AC_dyad_covar",
                "reciprocity_AC_dyad_covar_bin", "reciprocity_GW_dyad_covar",
                "reciprocity_GW_dyad_covar_bin",
                "in_weights_AC", "in_weights_GW",
                "transitivity_basic", "transitivity_min", "transitivity_netflow",
                "transitivity_AC", "transitivity_GW",
                "triad120C", "triad120D", "triad120U",
                "concentration_basic", "concentration_basic_cube",
                "concentration_prop", "concentration_prop_AC", "concentration_max",
                "reciprocity_basic", "reciprocity_min",
                "reciprocity_min_dyad_covar", "reciprocity_min_resource_covar",
                "same_covariate", "sim_covariate", "dyadic_covariate_resource_attribute",
                "in_ties_loops", "present_relations",
                "staying_by_prop_bin_inflow", "crowding_out_prop_covar_bin",
                "joining_similar_avoiding_dissimilar_covar_bin",
                "avoiding_dissimilar_covar_bin",
                "joining_similar_avoiding_dissimilar_covar_cont",
                "avoiding_dissimilar_covar_cont")

# cppEffectSpec
# translates an effects object into the specification used by the C++ effects:
# for each effect its name and arguments, with covariates given by their data.
# Returns NULL if an effect has no C++ version.
cppEffectSpec <- function(effects, state) {
  specs <- vector("list", length(effects$effectFormulas))
  for (k in seq_along(effects$effectFormulas)) {
    effectName <- strsplit(effects$name[k], " ")[[1]][1]
    if (!(effectName %in% cppEffects)) {
      return(NULL)
    }
    args <- as.list(formals(effects$effectFormulas[[k]]))
    spec <- list(name = effectName)
    if (!is.null(args$alpha)) {
      spec$alpha <- as.numeric(args$alpha)
    }
    if (!is.null(args$loop.contribution)) {
      spec$loopContribution <- as.logical(args$loop.contribution)
    }
    if (!is.null(args$attribute.index) && !identical(args$attribute.index, quote(expr = ))) {
      spec$attribute <- as.numeric(state[[args$attribute.index]]$data)
      spec$attributeIsMatrix <- is.matrix(state[[args$attribute.index]]$data)
      # same_covariate and sim_covariate use the 'same' / 'sim' matrix of the
      # covariate
      for (m in c("same", "sim")) {
        if (effectName == paste0(m, "_covariate")) {
          if (is.null(state[[args$attribute.index]][[m]])) {
            stop(paste0("Effect ", effectName, " expects covariates with a '", m, "' attribute."))
          }
          spec$attribute <- as.numeric(state[[args$attribute.index]][[m]])
          spec$attributeIsMatrix <- TRUE
        }
      }
    }
    if (!is.null(args$resource.attribute.index) &&
        !identical(args$resource.attribute.index, quote(expr = ))) {
      spec$resourceAttribute <- args$resource.attribute.index
    }
    specs[[k]] <- spec
  }
  specs
}

# useCppSampler
# whether the C++ sampler is used for these effects; options(MoNAn.useCpp =
# FALSE) switches to the R sampler
useCppSampler <- function(effects, state) {
  isTRUE(getOption("MoNAn.useCpp", TRUE)) && !is.null(cppEffectSpec(effects, state))
}

# cppChainInputs
# the arguments of the C++ functions that describe state and cache
cppChainInputs <- function(dep.var, state, cache) {
  state <- toListState(state, dep.var)
  cache <- toListCache(cache)
  edgelist <- state[[dep.var]]
  nodeName <- edgelist$nodeSet[1]
  resourceName <- edgelist$nodeSet[3]
  resourceNames <- names(cache[[dep.var]]$resourceNetworks)
  list(
    state = state,
    edgelist = matrix(as.integer(edgelist$data), ncol = 2),
    cache = cache[[dep.var]],
    resourceCovariates = stats::setNames(
      lapply(resourceNames, function(v) as.numeric(state[[v]]$data)), resourceNames
    ),
    resourceSet = as.integer(which(state[[resourceName]]$considerWhenSampling)),
    nodeSet = as.integer(state[[nodeName]]$ids[state[[nodeName]]$considerWhenSampling])
  )
}

# simulateChain
# runs one chain in C++; returns the statistics of each draw and the state and
# cache (as lists) at the end of the chain, and optionally of each draw
simulateChain <- function(dep.var, state, cache, effects, parameters,
                          burnIn, iterations, thinning, allowLoops,
                          multinomialProposal = FALSE, returnDraws = FALSE) {
  inputs <- cppChainInputs(dep.var, state, cache)
  res <- simulateChainCpp(
    inputs$edgelist, inputs$cache, inputs$resourceCovariates,
    inputs$resourceSet, inputs$nodeSet, cppEffectSpec(effects, inputs$state),
    as.numeric(parameters), as.integer(burnIn), as.integer(iterations),
    as.integer(thinning), allowLoops, multinomialProposal, returnDraws
  )

  toState <- function(edgelistData) {
    s <- inputs$state
    s[[dep.var]]$data <- edgelistData
    s
  }
  toCache <- function(cacheDepVar) {
    cacheDepVar <- cacheDepVar[names(inputs$cache)]
    if ("resourceNetworks" %in% names(inputs$cache) &&
        length(inputs$cache$resourceNetworks) == 0) {
      cacheDepVar$resourceNetworks <- inputs$cache$resourceNetworks
    }
    c <- list()
    c[[dep.var]] <- cacheDepVar
    c
  }

  statistics <- res$statistics
  colnames(statistics) <- effects$name
  out <- list(statistics = statistics,
              state = toState(res$edgelist),
              cache = toCache(res$cache))
  if (returnDraws) {
    out$draws <- lapply(res$draws, function(d) {
      list(state = toState(d$edgelist), cache = toCache(d$cache))
    })
  }
  out
}

# cppNetworkStatistics
# target statistics of a state in C++, as getNetworkStatistics
cppNetworkStatistics <- function(dep.var, state, cache, effects) {
  inputs <- cppChainInputs(dep.var, state, cache)
  statistics <- networkStatisticsCpp(inputs$edgelist, inputs$cache, inputs$resourceCovariates,
                                     cppEffectSpec(effects, inputs$state))
  names(statistics) <- effects$name
  statistics
}

# mpleObjectiveCpp
# as mpleObjectiveR, but calculated in C++ (mpleNewtonCpp in src/sampler.cpp),
# which recomputes the change statistics in each call instead of storing them
mpleObjectiveCpp <- function(dep.var, state, cache, effects, allowLoops) {
  inputs <- cppChainInputs(dep.var, state, cache)
  specs <- cppEffectSpec(effects, inputs$state)
  function(theta, derivatives = TRUE) {
    mpleNewtonCpp(inputs$edgelist, inputs$cache, inputs$resourceCovariates,
                  inputs$resourceSet, inputs$nodeSet, specs,
                  as.numeric(theta), allowLoops, derivatives)
  }
}

# cppEffectValues
# change statistic and target contribution of one C++ effect, for tests
cppEffectValues <- function(effects, k, state, cache, i, j, edge, update) {
  dep.var <- state$dep.var
  inputs <- cppChainInputs(dep.var, state, cache)
  effectValuesCpp(inputs$edgelist, inputs$cache, inputs$resourceCovariates,
                  cppEffectSpec(effects, inputs$state)[[k]],
                  as.integer(i), as.integer(j), as.integer(edge), as.integer(update))
}
