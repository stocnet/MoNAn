########## effectFunctions: effects concerning transitivity (Neighbourhood dependence)


#' transitivity_basic
#' 
#' Is mobility clustered in groups? This is represented by the count of
#' transitive triads among three nodes. This effect is prone to degeneracy.
#'
#' @param dep.var 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#'
#' 
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
transitivity_basic <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    if (i == j) {
      return(0)
    }
    
    net <- cache[[dep.var]]$valuedNetwork
    diag(net) <- 0
    
    if (getTargetContribution) {
      twoPaths <- apply(cbind(net[i, ], net[, j]), 1, prod)
      triadValues <- mapply(prod, net[i, j], twoPaths)
      return(sum(triadValues))
    }
    
    ### change contribution in three parts
    
    # part 1: i,j closes two-path
    
    twoPaths <- apply(cbind(net[i, ], net[, j]), 1, prod)
    ind_cont <- update * sum(twoPaths)
    
    # part 2: i,j closes instar
    
    inStars <- apply(cbind(net[i, ], net[j, ]), 1, prod)
    ind_cont <- ind_cont + update * sum(inStars)
    
    # part 3: i,j closes outstar
    
    outStars <- apply(cbind(net[, i], net[, j]), 1, prod)
    ind_cont <- ind_cont + update * sum(outStars)
    
    return(ind_cont)
}


#' transitivity_min
#' 
#' Is mobility clustered in groups? This is represented by the minimum of reciprocated 
#' mobility being present among three nodes. Using the minimum ensures that the effect 
#' is not degenerate and it is sample size consistent.
#'
#' @param dep.var 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#'
#' 
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
transitivity_min <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    if (i == j) {
      return(0)
    }
    
    m.min <- cache[[dep.var]]$minNetwork
    diag(m.min) <- 0
    
    if (getTargetContribution) {
      twoPaths <- apply(cbind(m.min[i, ], m.min[, j]), 1, min)
      triadValues <- mapply(min, m.min[i, j], twoPaths)
      v <- sum(triadValues) / 6
      return(v)
    }
    
    # cases without contributions
    if (update == 0) {
      return(0)
    }
    
    # change contributions are only possible if the minimum tie in the dyad is changed
    if (update > 0 &&
        cache[[dep.var]]$netFlowsNetwork[i, j] >= 0) {
      return(0)
    }
    if (update < 0 &&
        cache[[dep.var]]$netFlowsNetwork[i, j] > 0) {
      return(0)
    }
    
    # (+/-)1 change contributions only for those triads where the current tie is the weakest / unique weakest for negative / positive updates
    dyadValue <- m.min[i, j]
    twoPathValues <- apply(cbind(m.min[i, ], m.min[, j]), 1, min)
    if (update > 0) {
      v <- sum(dyadValue < twoPathValues)
    }
    if (update < 0) {
      v <- -sum(dyadValue <= twoPathValues)
    }
    return(v)
  }


#' transitivity_netflow
#' 
#' Do individuals move in one direction in locally ordered triads? E.g., is there 
#' a local hierarchy that individuals follow when moving between locations? The 
#' effect is sample size consistent.
#'
#' @param dep.var 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#'
#' 
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
transitivity_netflow <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    if (i == j) {
      return(0)
    }
    
    getPathWeights <- function(ik, kj) {
      if (ik < 0 && kj < 0) {
        return(0)
      }
      return(min(abs(ik), abs(kj)))
    }
    
    if (getTargetContribution) {
      if (cache[[dep.var]]$netFlowsNetwork[i, j] <= 0) {
        return(0)
      }
      pathWeights <-
        mapply(
          getPathWeights,
          cache[[dep.var]]$netFlowsNetwork[i, ],
          cache[[dep.var]]$netFlowsNetwork[, j]
        )
      v <-
        sum(mapply(min, pathWeights, cache[[dep.var]]$netFlowsNetwork[i, j]))
      return(v / 3)
    }
    
    # calculate contribution before change and relevant path weights (those relating to transitive structures)
    if (cache[[dep.var]]$netFlowsNetwork[i, j] > 0) {
      pathWeights <-
        mapply(
          getPathWeights,
          cache[[dep.var]]$netFlowsNetwork[i, ],
          cache[[dep.var]]$netFlowsNetwork[, j]
        )
      contributionBefore <-
        sum(mapply(min, pathWeights, cache[[dep.var]]$netFlowsNetwork[i, j]))
    }
    if (cache[[dep.var]]$netFlowsNetwork[i, j] < 0) {
      pathWeights <-
        mapply(
          getPathWeights,
          cache[[dep.var]]$netFlowsNetwork[j, ],
          cache[[dep.var]]$netFlowsNetwork[, i]
        )
      contributionBefore <-
        sum(mapply(min, pathWeights, cache[[dep.var]]$netFlowsNetwork[j, i]))
    }
    if (cache[[dep.var]]$netFlowsNetwork[i, j] == 0) {
      contributionBefore <- 0
      if (update > 0) {
        pathWeights <-
          mapply(
            getPathWeights,
            cache[[dep.var]]$netFlowsNetwork[i, ],
            cache[[dep.var]]$netFlowsNetwork[, j]
          )
      }
      if (update < 0) {
        pathWeights <-
          mapply(
            getPathWeights,
            cache[[dep.var]]$netFlowsNetwork[j, ],
            cache[[dep.var]]$netFlowsNetwork[, i]
          )
      }
    }
    
    # calculate contribution after the tie change
    tieWeightAfter <- cache[[dep.var]]$netFlowsNetwork[i, j] + update
    contributionAfter <-
      sum(mapply(min, pathWeights, abs(tieWeightAfter)))
    
    return(contributionAfter - contributionBefore)
  }


