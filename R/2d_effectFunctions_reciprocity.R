########## effectFunctions: effects concerning reciprocity (Edge dependence)


#' reciprocity_basic
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
#' @keywords internal
reciprocity_basic <- function(dep.var = 1, state, cache, i, j, edge, update, getTargetContribution = F){
  if (i == j) return(0)
  if (getTargetContribution){
    return (cache[[dep.var]]$valuedNetwork[i, j] * cache[[dep.var]]$valuedNetwork[j, i] / 2)
  } else {
    return(update*cache[[dep.var]]$valuedNetwork[j, i])
  }
}


#' reciprocity_min
#' 
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
#' @keywords internal
reciprocity_min <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
    if (i == j) {
      return(0)
    }
    
    if (getTargetContribution) {
      return(cache[[dep.var]]$minNetwork[i, j] / 2)
    }
    
    # simplified version that assumes that update are 1 or -1 and the network is an integer network
    if (update > 0 &&
        cache[[dep.var]]$netFlowsNetwork[i, j] < 0) {
      return(1)
    }
    if (update < 0 &&
        cache[[dep.var]]$netFlowsNetwork[i, j] <= 0) {
      return(-1)
    }
    return(0)
  }


#' reciprocity_min_resource_covar
#' 
#'
#' @param dep.var 
#' @param resource.attribute.index 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#'
#' 
#' @keywords internal
reciprocity_min_resource_covar <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
    if (i == j) {
      return(0)
    }
    
    if (getTargetContribution) {
      return(min(
        cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j],
        cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, i]
      ) / 2)
    }
    
    # simplified version that assumes that update are 1 or -1 and the network is an integer network
    if (update > 0 &&
        (cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j] <
         cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, i])) {
      return(1)
    }
    if (update < 0 &&
        (cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j] <=
         cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, i])) {
      return(-1)
    }
    return(0)
  }


#' reciprocity_GW
#'
#' @param dep.var 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#' @param lambda 
#'
#' @keywords internal
reciprocity_GW <- function(dep.var = 1, state, cache, i, j, edge, update, 
                           getTargetContribution = F, lambda = 2){
  if(lambda <= 0) stop("lambda parameter in reciprocity_GW function must be positive")
  if(i==j) return(0)
  if(getTargetContribution){
    nOutgoing <- cache[[dep.var]]$valuedNetwork[i, j]
    nIncoming <- cache[[dep.var]]$valuedNetwork[j, i]
    if(nOutgoing == 0) return(0)
    if(nIncoming == 0) return(0)
    v <- list()
    for(turn in 1:nIncoming){
      ind_cont <- 0
      for(k in 1:turn){
        ind_cont <- ind_cont + (1 / (lambda)^(k-1))
      }
      v[[turn]] <- ind_cont
    }
    return(nOutgoing * sum(unlist(v)))
  }
  ### calculate change statistic
  # first the part of the change statistic of the actor that moves
  nIncoming <- cache[[dep.var]]$valuedNetwork[j, i]
  if(nIncoming == 0) return(0)
  ind_cont <- 0
  for(k in 1:nIncoming){
    ind_cont <- ind_cont + (1 / (lambda)^(k-1))
  }
  # now the part of the change statistic for those that were already there
  nOutgoing <- cache[[dep.var]]$valuedNetwork[i, j]
  if(update == 1){
    ind_cont <- ind_cont + nIncoming * (1 / (lambda ^ (nOutgoing) ))
  }
  if(update == -1){
    ind_cont <- ind_cont + nIncoming * (1 / (lambda ^ (nOutgoing - 1) ))
  }
  return(update*ind_cont)
}

#' reciprocity_GW_dyad_covar_bin
#'
#' @param dep.var 
#' @param attribute.index 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param getTargetContribution 
#' @param lambda 
#'
#' @keywords internal
reciprocity_GW_dyad_covar_bin <- function(dep.var = 1, attribute.index, state, cache, i, j, edge, update, 
                                          getTargetContribution = F, lambda = 2){
  if(lambda <= 0) stop("lambda parameter in reciprocity_GW_dyad_covar_bin function must be positive")
  if(!all(state[[attribute.index]]$data == t(state[[attribute.index]]$data))) stop("attribute.index in reciprocity_GW_dyad_covar_bin function must be symmetric")
  if(!all(state[[attribute.index]]$data %in% c(0,1))) stop("all values of attribute.index in reciprocity_GW_dyad_covar_bin function must be 0 or 1")
  if(dim(state[[attribute.index]]$data)[1] != dim(state[[attribute.index]]$data)[2]) stop("attribute.index in reciprocity_GW_dyad_covar_bin function must be a square matrix")
  if(length(dim(state[[attribute.index]]$data)) != 2) stop("attribute.index in reciprocity_GW_dyad_covar_bin function must be a square matrix")
  if(i==j) return(0)
  if(getTargetContribution){
    nOutgoing <- cache[[dep.var]]$valuedNetwork[i, j] * state[[attribute.index]]$data[i,j]
    nIncoming <- cache[[dep.var]]$valuedNetwork[j, i] * state[[attribute.index]]$data[j,i]
    if(nOutgoing == 0) return(0)
    if(nIncoming == 0) return(0)
    v <- list()
    for(turn in 1:nIncoming){
      ind_cont <- 0
      for(k in 1:turn){
        ind_cont <- ind_cont + (1 / (lambda)^(k-1))
      }
      v[[turn]] <- ind_cont
    }
    return(nOutgoing * sum(unlist(v)))
  }
  ### calculate change statistic
  # first the part of the change statistic of the actor that moves
  nIncoming <- cache[[dep.var]]$valuedNetwork[j, i] * state[[attribute.index]]$data[j,i]
  if(nIncoming == 0) return(0)
  ind_cont <- 0
  for(k in 1:nIncoming){
    ind_cont <- ind_cont + (1 / (lambda)^(k-1))
  }
  # now the part of the change statistic for those that were already there
  nOutgoing <- cache[[dep.var]]$valuedNetwork[i, j] * state[[attribute.index]]$data[i,j]
  if(update == 1){
    ind_cont <- ind_cont + nIncoming * (1 / (lambda ^ (nOutgoing) ))
  }
  if(update == -1){
    ind_cont <- ind_cont + nIncoming * (1 / (lambda ^ (nOutgoing - 1) ))
  }
  return(update*ind_cont)
}




