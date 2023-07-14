########## effectFunctions: effects concerning concentration (Arc/tie dependence)

#' concentration_sq
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
concentration_sq <- function(dep.var = 1, state, cache, i, j, edge, update, getTargetContribution = F){
  if(getTargetContribution){
    return( (cache[[dep.var]]$valuedNetwork[i, j])^2 )
  } else {
    v <- (cache[[dep.var]]$valuedNetwork[i, j] + update)^2 - 
      (cache[[dep.var]]$valuedNetwork[i, j])^2
    return(v)
  }
}


#' concentration_GW
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
concentration_GW <- function(dep.var = 1, state, cache, i, j, edge, update, 
                             getTargetContribution = F, lambda = 2){
  if(lambda <= 0) stop("lambda parameter in concentration_GW function must be positive")
  if(getTargetContribution){
    nRessources <- cache[[dep.var]]$valuedNetwork[i, j]
    if(nRessources == 0) return(0)
    v <- list()
    for(turn in 1:nRessources){
      ind_cont <- 0
      for(k in 1:turn){
        ind_cont <- ind_cont + 2*(1 / (lambda)^(k))
      }
      v[[turn]] <- ind_cont
    }
    return(sum(unlist(v)))
  }
  ### calculate change statistic
  if (update > 0){
    ind_cont <- 0
    for(k in 1:(cache[[dep.var]]$valuedNetwork[i, j] + update)){
      ind_cont <- ind_cont + 2*(1 / (lambda)^(k))
    }
    return(ind_cont)
  }
  if (update < 0){
    ind_cont <- 0
    for(k in 1:(cache[[dep.var]]$valuedNetwork[i, j])){
      ind_cont <- ind_cont + 2*(1 / (lambda)^(k))
    }
    return(-ind_cont)
  }
}

#' concentration_GW_dyad_covar_bin
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
concentration_GW_dyad_covar_bin <- function(dep.var = 1, attribute.index, state, cache, i, j, edge, update, 
                                            getTargetContribution = F, lambda = 2){
  if(lambda <= 0) stop("lambda parameter in concentration_GW function must be positive")
  if(!all(state[[attribute.index]]$data == t(state[[attribute.index]]$data))) stop("attribute.index in concentration_GW_dyad_covar_bin function must be symmetric")
  if(!all(state[[attribute.index]]$data %in% c(0,1))) stop("all values of attribute.index in concentration_GW_dyad_covar_bin function must be 0 or 1")
  if(dim(state[[attribute.index]]$data)[1] != dim(state[[attribute.index]]$data)[2]) stop("attribute.index in concentration_GW_dyad_covar_bin function must be a square matrix")
  if(length(dim(state[[attribute.index]]$data)) != 2) stop("attribute.index in concentration_GW_dyad_covar_bin function must be a square matrix")
  
  if(getTargetContribution){
    nRessources <- cache[[dep.var]]$valuedNetwork[i, j] * state[[attribute.index]]$data[i, j]
    if(nRessources == 0) return(0)
    v <- list()
    for(turn in 1:nRessources){
      ind_cont <- 0
      for(k in 1:turn){
        ind_cont <- ind_cont + 2*(1 / (lambda)^(k))
      }
      v[[turn]] <- ind_cont
    }
    return(sum(unlist(v)))
  }
  ### calculate change statistic
  
  if(state[[attribute.index]]$data[i, j] == 0) return(0)
  
  if (update > 0){
    ind_cont <- 0
    for(k in 1:(cache[[dep.var]]$valuedNetwork[i, j] + update)){
      ind_cont <- ind_cont + 2*(1 / (lambda)^(k))
    }
    return(ind_cont)
  }
  if (update < 0){
    ind_cont <- 0
    for(k in 1:(cache[[dep.var]]$valuedNetwork[i, j])){
      ind_cont <- ind_cont + 2*(1 / (lambda)^(k))
    }
    return(-ind_cont)
  }
}

#' concentration_exponent
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
#' @param exponent 
#'
#' 
#' @keywords internal
concentration_exponent <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           exponent = .5) {
    if (i == j) {
      return(0)
    }
    
    if (getTargetContribution) {
      return((cache[[dep.var]]$valuedNetwork[i, j])^exponent)
    }
    
    v <- (cache[[dep.var]]$valuedNetwork[i, j] + update)^exponent -
      (cache[[dep.var]]$valuedNetwork[i, j])^exponent
    
    return(v)
  }


#' concentration_sigmoid
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
#' @param alpha 
#'
#' 
#' @keywords internal
concentration_sigmoid <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           alpha) {
    if (i == j) {
      return(0)
    }
    if (alpha <= 0) {
      stop("Alpha parameter in sigmoid tie weights function must be positive")
    }
    
    if (getTargetContribution) {
      nResources <- cache[[dep.var]]$valuedNetwork[i, j]
      v <- 0
      for (zz in 0:nResources) {
        v <- v + zz / (zz + alpha)
      }
      return(v)
    }
    
    value.old <- cache[[dep.var]]$valuedNetwork[i, j]
    value.new <- cache[[dep.var]]$valuedNetwork[i, j] + update
    if (update > 0) {
      v <- value.new / (value.new + alpha)
    }
    if (update < 0) {
      v <- -value.old / (value.old + alpha)
    }
    
    return(v)
  }

#' concentration_sigmoid_dyad_covar
#' 
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
#' @param alpha 
#'
#' 
#' @keywords internal
concentration_sigmoid_dyad_covar <-
  function(dep.var = 1,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           alpha) {
    if (i == j) {
      return(0)
    }
    if (alpha <= 0) {
      stop("Alpha parameter in sigmoid tie weights function must be positive")
    }
    
    if (getTargetContribution) {
      nResources <- cache[[dep.var]]$valuedNetwork[i, j]
      v <- 0
      for (zz in 0:nResources) {
        v <- v + zz / (zz + alpha)
      }
      return(v * state[[attribute.index]]$data[i, j])
    }
    
    value.old <- cache[[dep.var]]$valuedNetwork[i, j]
    value.new <- cache[[dep.var]]$valuedNetwork[i, j] + update
    if (update > 0) {
      v <- value.new / (value.new + alpha)
    }
    if (update < 0) {
      v <- -value.old / (value.old + alpha)
    }
    
    return(v * state[[attribute.index]]$data[i, j])
  }

