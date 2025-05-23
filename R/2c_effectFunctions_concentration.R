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

#' concentration_max
#' 
#' Is there a bandwagon effect in mobility following the most common path?
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
concentration_max <- function(dep.var = 1, state, cache, i, j, edge, update, getTargetContribution = FALSE){
  if(getTargetContribution){
    # each tie contributes the nth proportion of the row.max
    return( max(cache[[dep.var]]$valuedNetwork[i, ]) / 
              (dim(cache[[dep.var]]$valuedNetwork)[1]) )
  } else {
    v <- (cache[[dep.var]]$valuedNetwork[i, j])
    if(update == 1){
      if(max(cache[[dep.var]]$valuedNetwork[i, ]) == v){
        return(update)
      } 
    }
    if(update == -1){
      if(v > sort(cache[[dep.var]]$valuedNetwork[i, ], decreasing = T)[2]){
        return(update)
      }
    }
    return(0)
  }
}

#' concentration_basic_cube
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
concentration_basic_cube <- function(dep.var = 1, state, cache, i, j, edge, update, getTargetContribution = FALSE){
  if(getTargetContribution){
    return( (cache[[dep.var]]$valuedNetwork[i, j])^3 )
  } else {
    v <- (cache[[dep.var]]$valuedNetwork[i, j] + update)^3 - 
      (cache[[dep.var]]$valuedNetwork[i, j])^3
    return(v)
  }
}

#' concentration_AC
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
#' @param alpha 
#'
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
concentration_AC <- function(dep.var = 1, state, cache, i, j, edge, update, 
                             getTargetContribution = FALSE, alpha = 2){
  if(alpha < 1) stop("alpha parameter in concentration_AC function must be 1 or larger")
  
  if(i == j) return(0)
  
  if(getTargetContribution){
    
    g_cum <- function(y, a){
      contr <- 0
      
      #for(k in 0:y){
      #  contr <- contr + (y-k) * exp(-log(a)*k)
      #}
      #contr - y
      
      if(y>1){
        for(k in 1:(y-1)){
          contr <- contr + (1 - (1-1/a)^(k))
        }
      }
      return(contr)
    }
    
    nResources <- cache[[dep.var]]$valuedNetwork[i, j]
    
    return(g_cum(y = nResources, a = alpha))
  }
  
  ### calculate change statistic
  
  g_mar <- function(y, a){
    contr <- 0
    
    #for(k in 0:y){
    #  contr <- contr + exp(-log(a)*k)
    #}
    #contr - 1
    
    if(y>0) {
      contr <-  (1 - (1-1/a)^(y)) 
    } else {
      contr <- 0
    }
    return(contr)
  }
  
  tie_val <- cache[[dep.var]]$valuedNetwork[i, j]

  if(update < 0){
    return(update * g_mar(y = (tie_val + update), a = alpha))
  }
  if(update > 0){
    return(update * g_mar(y = tie_val, a = alpha))
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
#' @param alpha 
#'
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
concentration_GW <- function(dep.var = 1, state, cache, i, j, edge, update, 
                             getTargetContribution = FALSE, alpha = 2){
  if(alpha < 1) stop("alpha parameter in concentration_GW function must be 1 or larger")
  
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


#' concentration_AC_resource_covar_bin
#' 
#' Is there a bandwaggon effect in mobility, akin to the concentration_AC
#' effect, but where people of group 1 only consider others of group 1
#' in their decision to move?
#' 
#' @param dep.var 
#' @param state 
#' @param resource.attribute.index 
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
concentration_AC_resource_covar_bin <- function(dep.var = 1, 
                                                state,
                                                resource.attribute.index,
                                                cache, i, j, 
                                                edge,
                                                update, 
                                                getTargetContribution = FALSE, 
                                                alpha = 2){
  if(alpha < 1) stop("alpha parameter in concentration_AC_resource_covar_bin function must be 1 or larger")
  
  if(i == j) return(0)
  
  if(getTargetContribution){
    
    g_cum <- function(y, a){
      contr <- 0
      
      #for(k in 0:y){
      #  contr <- contr + (y-k) * exp(-log(a)*k)
      #}
      #contr - y
      
      if(y>1){
        for(k in 1:(y-1)){
          contr <- contr + (1 - (1-1/a)^(k))
        }
      }
      return(contr)
    }
    
    nResources <- cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j]
    
    return(g_cum(y = nResources, a = alpha))
  }
  
  ### calculate change statistic
  
  # if resource is of type 0, it cannot contribute
  
  if(state[[resource.attribute.index]]$data[edge] == 0){
    return(0)
  }
  
  g_mar <- function(y, a){
    contr <- 0
    
    #for(k in 0:y){
    #  contr <- contr + exp(-log(a)*k)
    #}
    #contr - 1
    
    if(y>0) {
      contr <-  (1 - (1-1/a)^(y)) 
    } else {
      contr <- 0
    }
    return(contr)
  }
  
  tie_val <- cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j]
  
  if(update < 0){
    return(update * g_mar(y = (tie_val + update), a = alpha))
  }
  if(update > 0){
    return(update * g_mar(y = tie_val, a = alpha))
  }
}

#' concentration_GW_resource_covar_bin
#' 
#' Is there a bandwaggon effect in mobility, akin to the concentration_GW
#' effect, but where people of group 1 only consider others of group 1
#' in their decision to move?
#' 
#' @param dep.var 
#' @param state 
#' @param resource.attribute.index 
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
concentration_GW_resource_covar_bin <- function(dep.var = 1, 
                                                state,
                                                resource.attribute.index,
                                                cache, i, j, 
                                                edge,
                                                update, 
                                                getTargetContribution = FALSE, 
                                                alpha = 2){
  if(alpha < 1) stop("alpha parameter in concentration_GW_resource_covar_bin function must be 1 or larger")
  
  if(i == j) return(0)
  
  if(getTargetContribution){
    
    g_cum <- function(y, a){
      contr <- 0
      for(k in 0:y){
        contr <- contr + (y-k) * exp(-log(a)*k)
      }
      contr - y
    }
    
    nResources <- cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j]
    
    return(g_cum(y = nResources, a = alpha))
  }
  
  ### calculate change statistic
  
  # if resource is of type 0, it cannot contribute
  
  if(state[[resource.attribute.index]]$data[edge] == 0){
    return(0)
  }
  
  g_mar <- function(y, a){
    contr <- 0
    for(k in 0:y){
      contr <- contr + exp(-log(a)*k)
    }
    contr - 1
  }
  
  tie_val <- cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j]
  
  if(update < 0){
    return(update * g_mar(y = (tie_val + update), a = alpha))
  }
  if(update > 0){
    return(update * g_mar(y = tie_val, a = alpha))
  }
}


#' concentration_AC_dyad_covar
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
#' @param alpha 
#'
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
concentration_AC_dyad_covar <- function(dep.var = 1, attribute.index, state, cache, i, j, edge, update, 
                                            getTargetContribution = FALSE, alpha = 2){
  if(alpha < 1) stop("alpha parameter in concentration_AC_dyad_covar function must be 1 or larger")
  # if(!all(state[[attribute.index]]$data == t(state[[attribute.index]]$data))) stop("attribute.index in concentration_AC_dyad_covar_bin function must be symmetric")
  # if(!all(state[[attribute.index]]$data %in% c(0,1))) stop("all values of attribute.index in concentration_AC_dyad_covar_bin function must be 0 or 1")
  if(dim(state[[attribute.index]]$data)[1] != dim(state[[attribute.index]]$data)[2]) stop("attribute.index in concentration_AC_dyad_covar_bin function must be a square matrix")
  if(length(dim(state[[attribute.index]]$data)) != 2) stop("attribute.index in concentration_AC_dyad_covar_bin function must be a square matrix")
  
  if(i == j) return(0)
  
  if(getTargetContribution){
    
    g_cum <- function(y, a){
      contr <- 0
      
      #for(k in 0:y){
      #  contr <- contr + (y-k) * exp(-log(a)*k)
      #}
      #contr - y
      
      if(y>1){
        for(k in 1:(y-1)){
          contr <- contr + (1 - (1-1/a)^(k))
        }
      }
      return(contr)
    }
    
    nResources <- cache[[dep.var]]$valuedNetwork[i, j] 
    
    return(g_cum(y = nResources, a = alpha) * state[[attribute.index]]$data[i, j])
  }
  
  ### calculate change statistic
  
  if(state[[attribute.index]]$data[i, j] == 0) return(0)
  
  g_mar <- function(y, a){
    contr <- 0
    
    #for(k in 0:y){
    #  contr <- contr + exp(-log(a)*k)
    #}
    #contr - 1
    
    if(y>0) {
      contr <-  (1 - (1-1/a)^(y)) 
    } else {
      contr <- 0
    }
    return(contr)
  }
  
  tie_val <- cache[[dep.var]]$valuedNetwork[i, j] 
  
  if(update < 0){
    return(update * g_mar(y = (tie_val + update), a = alpha) * state[[attribute.index]]$data[i, j])
  }
  if(update > 0){
    return(update * g_mar(y = tie_val, a = alpha) * state[[attribute.index]]$data[i, j])
  }
}


#' concentration_GW_dyad_covar
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
#' @param alpha 
#'
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
concentration_GW_dyad_covar <- function(dep.var = 1, attribute.index, state, cache, i, j, edge, update, 
                                        getTargetContribution = FALSE, alpha = 2){
  if(alpha < 1) stop("alpha parameter in concentration_GW function must be 1 or larger")
  # if(!all(state[[attribute.index]]$data == t(state[[attribute.index]]$data))) stop("attribute.index in concentration_GW_dyad_covar_bin function must be symmetric")
  # if(!all(state[[attribute.index]]$data %in% c(0,1))) stop("all values of attribute.index in concentration_GW_dyad_covar_bin function must be 0 or 1")
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
    
    nResources <- cache[[dep.var]]$valuedNetwork[i, j] 
    
    return(g_cum(y = nResources, a = alpha) * state[[attribute.index]]$data[i, j])
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
  
  tie_val <- cache[[dep.var]]$valuedNetwork[i, j] 
  
  if(update < 0){
    return(update * g_mar(y = (tie_val + update), a = alpha) * state[[attribute.index]]$data[i, j])
  }
  if(update > 0){
    return(update * g_mar(y = tie_val, a = alpha) * state[[attribute.index]]$data[i, j])
  }
}

## Currently under testing ##

#' concentration_prop
#'
#' Is there a bandwagon effect in mobility, i.e. do mobile individuals move to locations
#' that are the destination of many others from their origin? The functional form of this
#' statistic assumes that individuals consider the proportions of individuals (coming from
#' the same origin) going to a certain destination, instead of the total number.
#' This specification uses the idea that each individual has a contribution equal to 
#' the proportion of individuals from its origin going on the same tie as him/her.
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
concentration_prop <- function(dep.var = 1, state, cache, i, j, edge, update,
                               getTargetContribution = FALSE){
  ### if only one person is in the origin,
  ### this origin contributes 0 to the statistic
  if(sum(cache[[dep.var]]$valuedNetwork[i,]) < 1){
    return(0)
  }

  ### calculate target statistic
  if(getTargetContribution){
    stat <- (cache[[dep.var]]$valuedNetwork[i, j])^2 / 
      sum(cache[[dep.var]]$valuedNetwork[i,]) 
    return(stat)
  }

  ### calculate change statistic
  cont <- ((cache[[dep.var]]$valuedNetwork[i, j] + update)^2 - 
           (cache[[dep.var]]$valuedNetwork[i, j])^2) / 
          sum(cache[[dep.var]]$valuedNetwork[i,]) 
  return(cont)
}


#' concentration_prop_AC
#'
#' Is there a bandwagon effect in mobility, i.e. do mobile individuals move to locations
#' that are the destination of many others from their origin? The functional form of this
#' statistic assumes that individuals consider the proportions of individuals (coming from
#' the same origin) going to a certain destination, instead of the total number.
#' This specification uses the idea that each individual's contribution is a function of
#' the proportion of individuals from its origin going on the same tie as him/her.
#' The function is defined as an increasing function with decreasing positive gradient,
#' which means that the benefits of having additional individuals on a tie decrease with
#' more individuals. It's a similar logic to the previous AC statistics.
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
concentration_prop_AC <- function(dep.var = 1, state, cache, i, j, edge, update,
                               getTargetContribution = FALSE){
  
  if(alpha < 1) stop("alpha parameter in concentration_AC function must be 1 or larger")

  g <- function(k, a){
    return((1 - (1-1/a)^(k)))
  }
  
  ### calculate target statistic
  if(getTargetContribution){
    cont <- g(k = cache[[dep.var]]$valuedNetwork[i, j] / 
              sum(cache[[dep.var]]$valuedNetwork[i,]) ,
              a = alpha)
    stat <- cache[[dep.var]]$valuedNetwork[i, j] * cont
    return(stat)
  }
  
  ### calculate change statistic
  before <- (cache[[dep.var]]$valuedNetwork[i, j]) *
    g(k = cache[[dep.var]]$valuedNetwork[i, j] / 
      sum(cache[[dep.var]]$valuedNetwork[i,]) ,
      a = alpha) 
  if(update == -1){
    after <- (cache[[dep.var]]$valuedNetwork[i, j] - 1) *
              g(k = (cache[[dep.var]]$valuedNetwork[i, j]-1) / 
                sum(cache[[dep.var]]$valuedNetwork[i,]) ,
                a = alpha) 
    cont <- after - before
  }
  if(update == 1){
    after <- (cache[[dep.var]]$valuedNetwork[i, j] + 1) *
      g(k = (cache[[dep.var]]$valuedNetwork[i, j]+1) / 
          sum(cache[[dep.var]]$valuedNetwork[i,]) ,
        a = alpha) 
    cont <- after - before
  }
  return(cont)
}



## Deprecated ##
## If de-commenting, put back the names in the namespace ##

#' #' concentration_goodbad
#' #' 
#' #' Is there a bandwagon effect in mobility, i.e. do mobile individuals move to locations 
#' #' that are the destination of many others from their origin? The functional form of this 
#' #' statistic assumes that individuals consider the proportions of individuals (coming from
#' #' the same origin) going to a certain destination, instead of the total number.
#' #' This specification uses the idea of counting "good" and "bad" cliques.
#' #' 
#' #' @param dep.var 
#' #' @param state 
#' #' @param cache 
#' #' @param i 
#' #' @param j 
#' #' @param edge 
#' #' @param update 
#' #' @param getTargetContribution 
#' #'
#' #' @return Returns the change statistic or target statistic of the effect for 
#' #' internal use by the estimation algorithm. 
#' #' @keywords internal
#' concentration_goodbad <- function(dep.var = 1, state, cache, i, j, edge, update, 
#'                                getTargetContribution = FALSE){
#'   ### if only one person is in the origin, 
#'   ### this origin contributes 0 to the statistic
#'   if(sum(cache[[dep.var]]$valuedNetwork[i,]) < 2){
#'     return(0)
#'   }
#'   
#'   ### calculate target statistic
#'   if(getTargetContribution){
#'     numerator <- 2 * (cache[[dep.var]]$valuedNetwork[i, j])^2 - 
#'       cache[[dep.var]]$valuedNetwork[i, j] * sum(cache[[dep.var]]$valuedNetwork[i,]) - 
#'       cache[[dep.var]]$valuedNetwork[i, j]
#'     denominator <- sum(cache[[dep.var]]$valuedNetwork[i,]) - 1
#'     return( numerator/denominator )
#'   }
#'   
#'   ### calculate change statistic
#'   if(update == -1){
#'     cont <- 2*(-2*cache[[dep.var]]$valuedNetwork[i, j] + 1 + 
#'                  sum(cache[[dep.var]]$valuedNetwork[i,]) ) /
#'       (sum(cache[[dep.var]]$valuedNetwork[i,]) - 1)
#'   }
#'   if(update == 1){
#'     cont <- 2*(2*cache[[dep.var]]$valuedNetwork[i, j] - 
#'                  sum(cache[[dep.var]]$valuedNetwork[i,])) / 
#'       sum(cache[[dep.var]]$valuedNetwork[i,])
#'   }
#'   return(cont)
#' }
#' 
#' #' concentration_goodbad_orig_cov
#' #' 
#' #' Is there a bandwagon effect in mobility, i.e. do mobile individuals move to locations 
#' #' that are the destination of many others from their origin? The functional form of this 
#' #' statistic assumes that individuals consider the proportions of individuals (coming from
#' #' the same origin) going to a certain destination, instead of the total number.
#' #' This specification uses the idea of counting "good" and "bad" cliques.
#' #' 
#' #' This is weighted by an attribute of the origin, to model differences in 
#' #' concentration by origin characteristic.
#' #' 
#' #' @param dep.var 
#' #' @param attribute.index 
#' #' @param state 
#' #' @param cache 
#' #' @param i 
#' #' @param j 
#' #' @param edge 
#' #' @param update 
#' #' @param getTargetContribution 
#' #'
#' #' @return Returns the change statistic or target statistic of the effect for 
#' #' internal use by the estimation algorithm. 
#' #' @keywords internal
#' concentration_goodbad_orig_cov <- function(dep.var = 1, attribute.index, state, cache, i, j, edge, update, 
#'                                         getTargetContribution = FALSE){
#'   if(length(state[[attribute.index]]$nodeSet) != 1) stop("attribute in concentration_prop_orig_cov must be monadic")
#'   if(getTargetContribution){
#'     numerator <- 2 * (cache[[dep.var]]$valuedNetwork[i, j])^2 - 
#'       cache[[dep.var]]$valuedNetwork[i, j] * sum(cache[[dep.var]]$valuedNetwork[i,]) - 
#'       cache[[dep.var]]$valuedNetwork[i, j]
#'     denominator <- sum(cache[[dep.var]]$valuedNetwork[i,]) - 1
#'     if(denominator <= 0) {
#'       numerator <- 0
#'       denominator <- 1
#'     }
#'     return(state[[attribute.index]]$data[i] * numerator / denominator )
#'   }
#'   
#'   ### calculate change statistic
#'   
#'   if(update == -1){
#'     cont <- 2*(-2*cache[[dep.var]]$valuedNetwork[i, j] + 1 + sum(cache[[dep.var]]$valuedNetwork[i,]) ) /
#'       (sum(cache[[dep.var]]$valuedNetwork[i,]) - 1)
#'   }
#'   if(update == 1){
#'     cont <- 2*(2*cache[[dep.var]]$valuedNetwork[i, j] - sum(cache[[dep.var]]$valuedNetwork[i,])) / 
#'       sum(cache[[dep.var]]$valuedNetwork[i,])
#'   }
#'   return(state[[attribute.index]]$data[i] * cont)
#' }
#' 
#' #' concentration_rankGW
#' #' 
#' #' Is there a bandwagon effect in mobility, i.e. do mobile individuals move to locations 
#' #' that are the destination of many others from their origin? The functional form of this 
#' #' statistic assumes that there are increasing returns to choosing more populated paths
#' #' but the return of using a path only depends on the ranking of this path in terms of number 
#' #' of people using it among all paths (and not the actual numbers of people on this path).
#' #' 
#' #' @param dep.var 
#' #' @param state 
#' #' @param cache 
#' #' @param i 
#' #' @param j 
#' #' @param edge 
#' #' @param update 
#' #' @param getTargetContribution 
#' #' @param lambda 
#' #'
#' #' @return Returns the change statistic or target statistic of the effect for 
#' #' internal use by the estimation algorithm.
#' #' @keywords internal
#' concentration_rankGW <- function(dep.var = 1, state, cache, i, j, edge, update, 
#'                                  getTargetContribution = FALSE, lambda = 2){
#'   if(lambda <= 0) stop("lambda parameter in concentration_GW function must be positive")
#'   if(getTargetContribution){
#'     nRessources <- cache[[dep.var]]$valuedNetwork[i, j]
#'     if(nRessources == 0) return(0)
#'     rank <- rank(-cache[[dep.var]]$valuedNetwork[i,])[j]
#'     v <- nRessources * exp(rank*log(1 / (lambda)))
#'     return(v)
#'   }
#'   ### calculate change statistic
#'   ranks1 <- rank(-cache[[dep.var]]$valuedNetwork[i,])
#'   newnet <- cache[[dep.var]]$valuedNetwork
#'   newnet[i,j] <- newnet[i,j] + update
#'   ranks2 <- rank(-newnet[i,])
#'   if(ranks1[j] == ranks2[j]){
#'     v <- newnet[i, j] * exp(ranks2[j] * log(1 / (lambda))) -
#'       cache[[dep.var]]$valuedNetwork[i, j] * exp(ranks1[j] * log(1 / (lambda)))
#'   } else {
#'     changedranks <- which( (ranks1-ranks2) !=0 )
#'     v <- 0
#'     for(changed in changedranks) {
#'       v <- v + newnet[i, changed] * exp(ranks2[changed] * log(1 / (lambda))) -
#'         cache[[dep.var]]$valuedNetwork[i, changed] * exp(ranks1[changed] * log(1 / (lambda)))
#'     }
#'     return(v)
#'   }
#' }
#' 
#' 
#' #' concentration_basic_squared
#' #' 
#' #' Is there a bandwagon effect in mobility, i.e. do mobile individuals move to 
#' #' locations that are the destination of many others from their origin?
#' #' 
#' #' @param dep.var 
#' #' @param state 
#' #' @param cache 
#' #' @param i 
#' #' @param j 
#' #' @param edge 
#' #' @param update 
#' #' @param getTargetContribution 
#' #'
#' #' @return Returns the change statistic or target statistic of the effect for 
#' #' internal use by the estimation algorithm.
#' #' @keywords internal
#' concentration_basic_squared <- function(dep.var = 1, state, cache, i, j, edge, update, getTargetContribution = FALSE){
#'   if(getTargetContribution){
#'     return( sum(cache[[dep.var]]$valuedNetwork[i,j]^2 * 
#'                   cache[[dep.var]]$valuedNetwork[i,]^2) )
#'   } else {
#'     previous_rows <- cache[[dep.var]]$valuedNetwork[i,]
#'     new_rows <- previous_rows
#'     new_rows[j] <- new_rows[j] + update
#'     v <- (sum(new_rows^2))^2 - (sum(previous_rows^2))^2
#'     return(v)
#'   }
#' }
#' 
#' 
#' 
#' #' concentration_norm
#' #' 
#' #' Is there a bandwagon effect in mobility, i.e. do mobile individuals move to 
#' #' locations that are the destination of many others from their origin?
#' #' 
#' #' @param dep.var 
#' #' @param state 
#' #' @param cache 
#' #' @param i 
#' #' @param j 
#' #' @param edge 
#' #' @param update 
#' #' @param getTargetContribution 
#' #'
#' #' @return Returns the change statistic or target statistic of the effect for 
#' #' internal use by the estimation algorithm.
#' #' @keywords internal
#' concentration_norm <- function(dep.var = 1, state, cache, i, j, edge, update, getTargetContribution = FALSE){
#'   if(getTargetContribution){
#'     return( cache[[dep.var]]$valuedNetwork[i,j]^2 / 
#'               (sum(cache[[dep.var]]$valuedNetwork[i,]) - 1) )
#'   } else {
#'     if(update == -1){
#'       cont <- ((cache[[dep.var]]$valuedNetwork[i,j]-1)^2 -
#'                  cache[[dep.var]]$valuedNetwork[i,j]^2 ) /
#'         (sum(cache[[dep.var]]$valuedNetwork[i,]) - 1)
#'     }
#'     if(update == 1){
#'       cont <- ((cache[[dep.var]]$valuedNetwork[i,j]+1)^2 -
#'                  cache[[dep.var]]$valuedNetwork[i,j]^2 ) /
#'         sum(cache[[dep.var]]$valuedNetwork[i,])
#'     }
#'     return(cont)
#'   }
#' }
#' 
#' #' concentration_norm_squared
#' #' 
#' #' Is there a bandwagon effect in mobility, i.e. do mobile individuals move to 
#' #' locations that are the destination of many others from their origin?
#' #' 
#' #' @param dep.var 
#' #' @param state 
#' #' @param cache 
#' #' @param i 
#' #' @param j 
#' #' @param edge 
#' #' @param update 
#' #' @param getTargetContribution 
#' #'
#' #' @return Returns the change statistic or target statistic of the effect for 
#' #' internal use by the estimation algorithm.
#' #' @keywords internal
#' concentration_norm_squared <- function(dep.var = 1, state, cache, i, j, edge, update, getTargetContribution = FALSE){
#'   if(getTargetContribution){
#'     return( sum(cache[[dep.var]]$valuedNetwork[i,j]^2 * cache[[dep.var]]$valuedNetwork[i,]^2 / 
#'                   (sum(cache[[dep.var]]$valuedNetwork[i,]) - 1)^2 ))
#'   } else {
#'     previous_rows <- cache[[dep.var]]$valuedNetwork[i,]
#'     new_rows <- previous_rows
#'     if(update == -1){
#'       new_rows[j] <- new_rows[j] - 1
#'       cont <- ((sum(new_rows^2 / (sum(cache[[dep.var]]$valuedNetwork[i,]) - 1)))^2 - 
#'                  (sum(previous_rows^2 / (sum(cache[[dep.var]]$valuedNetwork[i,]) - 1)))^2) 
#'     }
#'     if(update == 1){
#'       new_rows[j] <- new_rows[j] + 1
#'       cont <- ((sum(new_rows^2 / (sum(cache[[dep.var]]$valuedNetwork[i,]))))^2 - 
#'                  (sum(previous_rows^2 / (sum(cache[[dep.var]]$valuedNetwork[i,]))))^2) 
#'     }
#'     return(cont)
#'   }
#' }
#' 
#' 
