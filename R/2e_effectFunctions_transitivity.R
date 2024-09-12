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


#' transitivity_AC
#'
#' Is mobility clustered in groups? This is represented by the count of
#' transitive triads among three nodes, where the number of two-paths is
#' geometrically weighted down to avoid degeneracy.
#' 
#' @param dep.var 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#' @param alpha 
#'
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
transitivity_AC <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE,
           alpha = 1.1) {
    
    if (i == j) {
      return(0)
    }
    
    g_mar <- function(y, a){
      contr <- 0
      if(y>0) {
        contr <-  (1 - (1-1/a)^(y)) 
      } else {
        contr <- 0
      }
      return(contr)
    }
    
    net <- cache[[dep.var]]$valuedNetwork
    diag(net) <- 0
    
    if (getTargetContribution) {
      twoPaths <- sum(apply(cbind(net[i, ], net[, j]), 1, prod))
      return(net[i,j] * g_mar(twoPaths, alpha))
    }
    
    ### change contribution in three parts
    
    # part 1: i,j closes two-path
    
    twoPaths <- sum(apply(cbind(net[i, ], net[, j]), 1, prod))
    ind_cont <- update * g_mar(twoPaths, alpha)
    
    # part 2: i,j closes instar
    net_u <- net
    net_u[i,j] <- net_u[i,j] + update
    
    nTP_before <- (net[i,] %*% net)
    con_be <- sum(net[i,] * unlist(lapply(nTP_before, g_mar, alpha))) - (net[i,j] * g_mar(nTP_before[j], alpha))
    nTP_after <- (net_u[i,] %*% net_u)
    con_af <- sum(net_u[i,] * unlist(lapply(nTP_after, g_mar, alpha))) - (net_u[i,j] * g_mar(nTP_after[j], alpha))
    ind_cont <- ind_cont + (con_af - con_be)
    
    # part 3: i,j closes outstar
    
    nTP_before <- (net %*% net[,j])
    con_be <- sum(net[,j] * unlist(lapply(nTP_before, g_mar, alpha)))
    nTP_after <- (net_u %*% net_u[,j])
    con_af <- sum(net[,j] * unlist(lapply(nTP_after, g_mar, alpha)))
    ind_cont <- ind_cont + (con_af - con_be)
    
    return(ind_cont)
  }


#' transitivity_GW
#'
#' Is mobility clustered in groups? This is represented by the count of
#' transitive triads among three nodes, where the number of two-paths is
#' geometrically weighted down to avoid degeneracy.
#' 
#' @param dep.var 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#' @param alpha 
#'
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
transitivity_GW <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE,
           alpha = 1.1) {
    
    if (i == j) {
      return(0)
    }
    
    g_mar <- function(y, a){
      contr <- 0
      for(k in 0:y){
        contr <- contr + exp(-log(a)*k)
      }
      contr - 1
    }
    
    net <- cache[[dep.var]]$valuedNetwork
    diag(net) <- 0
    
    if (getTargetContribution) {
      twoPaths <- sum(apply(cbind(net[i, ], net[, j]), 1, prod))
      return(net[i,j] * g_mar(twoPaths, alpha))
    }
    
    ### change contribution in three parts
    
    # part 1: i,j closes two-path
    
    twoPaths <- sum(apply(cbind(net[i, ], net[, j]), 1, prod))
    ind_cont <- update * g_mar(twoPaths, alpha)
    
    # part 2: i,j closes instar
    net_u <- net
    net_u[i,j] <- net_u[i,j] + update
    
    nTP_before <- (net[i,] %*% net)
    con_be <- sum(net[i,] * unlist(lapply(nTP_before, g_mar, alpha))) - (net[i,j] * g_mar(nTP_before[j], alpha))
    nTP_after <- (net_u[i,] %*% net_u)
    con_af <- sum(net_u[i,] * unlist(lapply(nTP_after, g_mar, alpha))) - (net_u[i,j] * g_mar(nTP_after[j], alpha))
    ind_cont <- ind_cont + (con_af - con_be)
    
    # part 3: i,j closes outstar
    
    nTP_before <- (net %*% net[,j])
    con_be <- sum(net[,j] * unlist(lapply(nTP_before, g_mar, alpha)))
    nTP_after <- (net_u %*% net_u[,j])
    con_af <- sum(net[,j] * unlist(lapply(nTP_after, g_mar, alpha)))
    ind_cont <- ind_cont + (con_af - con_be)
    
    return(ind_cont)
  }


#' triad120D
#' 
#' Models the prevalence of the 120D triad
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#'
#' @keywords internal
triad120D <- function(dep.var = 1,
                      state,
                      cache,
                      i,
                      j,
                      edge,
                      update,
                      getTargetContribution = FALSE){
  
  if (i == j) {
    return(0)
  }
  
  # the target contribution is calculated for the i-j pair that is reciprocated
  net <- cache[[dep.var]]$valuedNetwork
  diag(net) <- 0
  
  if (getTargetContribution) {
    n_down_stars <- apply(cbind(net[, i], net[, j]), 1, prod)
    triadValues <- mapply(prod, (net[i, j] * net[j, i]), n_down_stars)
    return(sum(triadValues) / 2)
  }
  
  # for the change stat, four contributions are needed, one for each tie var,
  # but two are isomorphic
  
  # tie on the reciprocated dyad
  n_down_stars <- apply(cbind(net[, i], net[, j]), 1, prod)
  ind_cont <- update * sum(n_down_stars) * net[j,i]
  
  # tie on the down star
  recip_net <- net * t(net)
  n_recip_in_stars <- apply(cbind(net[i, ], recip_net[, j]), 1, prod)
  ind_cont <- ind_cont + update * sum(n_recip_in_stars)
  
  return(ind_cont)
}


#' triad120U
#' 
#' Models the prevalence of the 120U triad
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#'
#' @keywords internal
triad120U <- function(dep.var = 1,
                      state,
                      cache,
                      i,
                      j,
                      edge,
                      update,
                      getTargetContribution = FALSE){
  
  if (i == j) {
    return(0)
  }
  
  # the target contribution is calculated for the i-j pair that is reciprocated
  net <- cache[[dep.var]]$valuedNetwork
  diag(net) <- 0
  
  if (getTargetContribution) {
    n_up_stars <- apply(cbind(net[i, ], net[j, ]), 1, prod)
    triadValues <- mapply(prod, (net[i, j] * net[j, i]), n_up_stars)
    return(sum(triadValues) / 2)
  }
  
  # for the change stat, four contributions are needed, one for each tie var,
  # but two are isomorphic
  
  # tie on the reciprocated dyad
  n_up_stars <- apply(cbind(net[i, ], net[j, ]), 1, prod)
  ind_cont <- update * sum(n_up_stars) * net[j,i]
  
  # tie on the down star
  recip_net <- net * t(net)
  n_recip1_two_paths <- apply(cbind(net[, j], recip_net[i, ]), 1, prod)
  ind_cont <- ind_cont + update * sum(n_recip1_two_paths)
  
  return(ind_cont)
}


#' triad120C
#' 
#' Models the prevalence of the 120C triad
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#'
#' @keywords internal
triad120C <- function(dep.var = 1,
                      state,
                      cache,
                      i,
                      j,
                      edge,
                      update,
                      getTargetContribution = FALSE){
  
  if (i == j) {
    return(0)
  }
  
  # the target contribution is calculated for the i-j pair that is reciprocated
  net <- cache[[dep.var]]$valuedNetwork
  diag(net) <- 0
  
  if (getTargetContribution) {
    n_two_paths <- apply(cbind(net[i, ], net[, j]), 1, prod)
    triadValues <- mapply(prod, (net[i, j] * net[j, i]), n_two_paths)
    return(sum(triadValues))
  }
  
  # for the change stat, four contributions are needed, one for each tie var

  # tie on the reciprocated dyad; transitive closure
  n_two_paths <- apply(cbind(net[i, ], net[, j]), 1, prod)
  ind_cont <- update * sum(n_two_paths) * net[j,i]

  # tie on the reciprocated dyad; cyclic closure
  n_rev_two_paths <- apply(cbind(net[, i], net[j, ]), 1, prod)
  ind_cont <- ind_cont + update * sum(n_rev_two_paths) * net[j,i]
  
  # tie on the first part of two-path
  recip_net <- net * t(net)
  n_recip1_two_paths <- apply(cbind(net[j, ], recip_net[i, ]), 1, prod)
  ind_cont <- ind_cont + update * sum(n_recip1_two_paths)
  
  # tie on the second part of two-path
  n_recip2_two_paths <- apply(cbind(net[, i], recip_net[, j]), 1, prod)
  ind_cont <- ind_cont + update * sum(n_recip2_two_paths)
  
  return(ind_cont)
}


