########## effectFunctions: effects concerning loops

#' loops
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
#' @return None.
#' @keywords internal
loops <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
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
#' @return None.
#' @keywords internal
loops_GW <- function(dep.var = 1,
                     state,
                     cache,
                     i,
                     j,
                     edge,
                     update,
                     getTargetContribution = F,
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


#' loops_sigmoid
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
#' @return None.
#' @keywords internal
loops_sigmoid <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           alpha) {
    if (alpha <= 0) {
      stop("Alpha parameter in sigmoid loops weights function must be positive")
    }
    
    if (getTargetContribution) {
      if (i == j) {
        nResources <- cache[[dep.var]]$valuedNetwork[i, i]
        v <- 0
        for (zz in 0:nResources) {
          v <- v + zz / (zz + alpha)
        }
        return(v)
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
        v <- value.new / (value.new + alpha)
      }
      if (update < 0) {
        v <- -value.old / (value.old + alpha)
      }
      
      return(v)
    }
  }


#' loops_GW_prop
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
#' @param half.cont 
#'
#' 
#' @return None.
#' @keywords internal
loops_GW_prop <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           half.cont = 0.1) {
    if (half.cont <= 0) {
      stop("half.cont parameter in GW loopsNegExp function must be positive")
    }
    
    getExpChange <- function(count, sum) {
      alpha <- -(1 / half.cont) * log(0.5)
      return(1 - exp(-alpha * count / sum))
    }
    
    if (getTargetContribution) {
      if (i == j) {
        nLoops <- cache[[dep.var]]$valuedNetwork[i, i]
        if (nLoops == 0) {
          return(0)
        }
        
        resInRow <- sum(cache[[dep.var]]$valuedNetwork[i, ])
        return(sum(getExpChange(1:nLoops, resInRow)))
      }
      return(0)
    }
    
    if (i != j) {
      return(0)
    }
    
    if (i == j) {
      nLoops <- cache[[dep.var]]$valuedNetwork[i, i]
      resInRow <- sum(cache[[dep.var]]$valuedNetwork[i, ])
      
      if (update == 1) {
        return(getExpChange((nLoops + 1), resInRow))
      }
      if (update == -1) {
        return(-getExpChange((nLoops), resInRow))
      }
    }
  }



#' loops_node_covar
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
#'
#' 
#' @return None.
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
           getTargetContribution = F) {
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
#' @return None.
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
           getTargetContribution = F) {
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
#' @return None.
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
           getTargetContribution = F) {
    if (i != j) {
      return(0)
    }
    
    if (getTargetContribution) {
      return(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j])
    }
    
    return(update * state[[resource.attribute.index]]$data[edge])
  }






