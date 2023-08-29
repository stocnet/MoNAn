########## effectFunctions: effects concerning loops

#' loops
#' 
#' Do individuals stay in their location of origin, compared to going to a different location?
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
loops <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    if (getTargetContribution) {
      return(cache[[dep.var]]$valuedNetwork[i, j] * (i == j))
    }
    
    if (i == j) {
      return(update)
    }
    return(0)
  }


#' loops_GW
#' 
#' Do individuals stay in their current location, in case many other from their 
#' current location also stay? This effect tests whether the ‘benefit’ of staying 
#' in the origin location depends on the number of others also staying. Note that 
#' this effect should be modelled alongside the loops effect.
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
loops_GW <- function(dep.var = 1,
                     state,
                     cache,
                     i,
                     j,
                     edge,
                     update,
                     getTargetContribution = FALSE,
                     alpha = 2) {
  if (alpha <= 0) {
    stop("Alpha parameter in GW loops weights function must be positive")
  }
  
  if (getTargetContribution) {
    if (i == j) {
      nResources <- cache[[dep.var]]$valuedNetwork[i, i]
      
      if (nResources == 0) {
        return(0)
      }
      
      v <- list()
      for (turn in 1:nResources) {
        ind_cont <- 0
        for (k in 1:turn) {
          ind_cont <- ind_cont + (1 / (alpha)^(k))
        }
        v[[turn]] <- ind_cont
      }
      return(sum(unlist(v)))
    }
    
    return(0)
  }
  
  if (i != j) {
    return(0)
  }
  
  if (i == j) {
    value.old <- cache[[dep.var]]$valuedNetwork[i, i]
    value.new <- cache[[dep.var]]$valuedNetwork[i, i] + update
    
    if (update > 0) {
      ind_cont <- 0
      for (k in 1:value.new) {
        ind_cont <- ind_cont + (1 / (alpha)^(k))
      }
      v <- ind_cont
    }
    
    if (update < 0) {
      ind_cont <- 0
      for (k in 1:value.old) {
        ind_cont <- ind_cont + (1 / (alpha)^(k))
      }
      v <- -ind_cont
    }
    
    return(v)
  }
}


#' loops_node_covar
#' 
#' Are locations with specific attributes ‘stickier’ than others, i.e., do individuals 
#' have a higher propensity to stay in some locations? E.g., are individuals working 
#' in organisations in one region less likely to change their employer?
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
#'
#' 
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
loops_node_covar <-
  function(dep.var = 1,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    if (i != j) {
      return(0)
    }
    
    if (getTargetContribution) {
      return(cache[[dep.var]]$valuedNetwork[i, j] * state[[attribute.index]]$data[i])
    }
    
    return(update * state[[attribute.index]]$data[i])
  }


#' loops_resource_covar_node_covar
#' 
#' This is an interaction of the previous two effects: Do individuals with certain 
#' characteristics have a tendency to stay in locations of certain types? Note that 
#' this effect should be included alongside the main effects of ‘loops by individual 
#' covariate’ and ‘loops by location covariate’.
#'
#' @param dep.var 
#' @param resource.attribute.index 
#' @param attribute.index 
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
loops_resource_covar_node_covar <-
  function(dep.var = 1,
           resource.attribute.index,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    # the count in the resource network on the diagonal times the node covar value
    if (getTargetContribution) {
      return((cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j] * (i == j)) * state[[attribute.index]]$data[i])
    }
    
    # for loops; the resource attribute; times the node attribute
    if (i == j) {
      return(update * state[[resource.attribute.index]]$data[edge] * state[[attribute.index]]$data[j])
    }
    return(0)
  }


#' loops_resource_covar
#'
#' Are individuals with certain characteristics more likely to remain in their current 
#' location? For example, are men more likely to remain in their current organisation, 
#' while women are more likely to move employer? 
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
loops_resource_covar <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    if (i != j) {
      return(0)
    }
    
    if (getTargetContribution) {
      return(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j])
    }
    
    return(update * state[[resource.attribute.index]]$data[edge])
  }






