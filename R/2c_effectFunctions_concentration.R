########## effectFunctions: effects concerning concentration (Arc/tie dependence)

#' concentration_basic
#' 
#' Is there a bandwagon effect in mobility, i.e. do mobile individuals move to 
#' locations that are the destination of many others from their origin?
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
#' @keywords internal
concentration_basic <- function(dep.var = 1, state, cache, i, j, edge, update, getTargetContribution = FALSE){
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
#' Is there a bandwagon effect in mobility, i.e. do mobile individuals move to locations 
#' that are the destination of many others from their origin? The functional form of this 
#' statistic assumes that there are decreasing additional returns to more others on the 
#' same mobility path. For example, the probability to choose a mobility path that already 
#' contains 20 other individuals is hardly different from a path with 25 other individuals; 
#' however, there is a substantial difference in the comparison of paths with 2 or 7 other 
#' individuals.
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
concentration_GW <- function(dep.var = 1, state, cache, i, j, edge, update, 
                             getTargetContribution = FALSE, lambda = 2){
  if(lambda <= 0) stop("lambda parameter in concentration_GW function must be positive")
  
  if(i == j) return(0)
  
  if(getTargetContribution){
    
    g_cum <- function(y, a){
      contr <- 0
      for(k in 0:y){
        contr <- contr + (y-k) * exp(-log(a)*k)
      }
      contr - y
    }
    
    nResources <- cache[[dep.var]]$valuedNetwork[i, j]
    
    return(g_cum(y = nResources, a = alpha))
  }
  
  ### calculate change statistic
  
  g_mar <- function(y, a){
    contr <- 0
    for(k in 0:y){
      contr <- contr + exp(-log(a)*k)
    }
    contr - 1
  }
  
  tie_val <- cache[[dep.var]]$valuedNetwork[i, j]

  if(update < 0){
    return(update * g_mar(y = (tie_val + update), a = alpha))
  }
  if(update > 0){
    return(update * g_mar(y = tie_val, a = alpha))
  }
}


#' concentration_GW_dyad_covar_bin
#' 
#' Are bandwagon effects (concentration) particularly prevalent between locations 
#' that share characteristics as encoded in a binary dyadic covariate? E.g., do 
#' workers follow the moves of other workers mainly in case they go to organisations 
#' in the same region?
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
concentration_GW_dyad_covar_bin <- function(dep.var = 1, attribute.index, state, cache, i, j, edge, update, 
                                            getTargetContribution = FALSE, lambda = 2){
  if(lambda <= 0) stop("lambda parameter in concentration_GW function must be positive")
  if(!all(state[[attribute.index]]$data == t(state[[attribute.index]]$data))) stop("attribute.index in concentration_GW_dyad_covar_bin function must be symmetric")
  if(!all(state[[attribute.index]]$data %in% c(0,1))) stop("all values of attribute.index in concentration_GW_dyad_covar_bin function must be 0 or 1")
  if(dim(state[[attribute.index]]$data)[1] != dim(state[[attribute.index]]$data)[2]) stop("attribute.index in concentration_GW_dyad_covar_bin function must be a square matrix")
  if(length(dim(state[[attribute.index]]$data)) != 2) stop("attribute.index in concentration_GW_dyad_covar_bin function must be a square matrix")
  
  if(i == j) return(0)
  
  if(getTargetContribution){
    
    g_cum <- function(y, a){
      contr <- 0
      for(k in 0:y){
        contr <- contr + (y-k) * exp(-log(a)*k)
      }
      contr - y
    }
    
    nResources <- cache[[dep.var]]$valuedNetwork[i, j] * state[[attribute.index]]$data[i, j]
    
    return(g_cum(y = nResources, a = alpha))
  }
  
  ### calculate change statistic
  
  if(state[[attribute.index]]$data[i, j] == 0) return(0)
  
  g_mar <- function(y, a){
    contr <- 0
    for(k in 0:y){
      contr <- contr + exp(-log(a)*k)
    }
    contr - 1
  }
  
  tie_val <- cache[[dep.var]]$valuedNetwork[i, j] * state[[attribute.index]]$data[i, j]
  
  if(update < 0){
    return(update * g_mar(y = (tie_val + update), a = alpha))
  }
  if(update > 0){
    return(update * g_mar(y = tie_val, a = alpha))
  }
}


