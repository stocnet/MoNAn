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
  
  g_cum <- function(y, a){
    contr <- 0
    if(y>1){
      for(k in 1:(y-1)){
        contr <- contr + (1 - (1-1/a)^(k))
      }
    }
    return(contr)
  }
  
  if (getTargetContribution) {
    if (i == j) {
      nResources <- cache[[dep.var]]$valuedNetwork[i, i]
      
      return(g_cum(y = nResources, a = alpha))
    }
    
    return(0)
  }
  
  if (i != j) {
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
  
  if (i == j) {
    tie_val <- cache[[dep.var]]$valuedNetwork[i, i]
    if(update < 0){
      return(update * g_mar(y = (tie_val + update), a = alpha))
    }
    if(update > 0){
      return(update * g_mar(y = tie_val, a = alpha))
    }
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

#' loops_additional_origin
#' 
#' This effect models loops for cases in which individuals have more than one 
#' origin. The additional origin not specified in the mobility data is included
#' as a resource.attribute.index.
#' The question modeled is: Do individuals stay in the additional location of 
#' origin, compared to going to a different location?
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
loops_additional_origin <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    
    if (getTargetContribution) {
      res_index <- (state[[dep.var]]$data[,1] == i) * (state[[dep.var]]$data[,2] == j)
      cont <- sum(res_index * (state[[resource.attribute.index]]$data == j))
      return(cont)
    }
    
    dest.of.res <- j
    additional.orig.of.res <- state[[resource.attribute.index]]$data[edge]
    
    if(dest.of.res == additional.orig.of.res){
      return(update)
    }
    
    return(0)
    
  }


#' loops_x_loops_additional_origin
#' 
#' This effect is specified for cases in which individuals have more than one 
#' origin. The additional origin not specified in the mobility data is included
#' as a resource.attribute.index.
#' The question modeled is: Do individuals stay in the additional location of 
#' origin if this is additionally their origin as specified in the mobility data, 
#' compared to going to a different location?
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
loops_x_loops_additional_origin <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    
    if (getTargetContribution) {
      if(i != j){
        return(0)
      }
      res_index <- (state[[dep.var]]$data[,1] == i) * (state[[dep.var]]$data[,2] == j)
      cont <- sum(res_index * (state[[resource.attribute.index]]$data == j))
      return(cont)
    }
    
    orig.of.res <- state[[dep.var]]$data[edge,1]
    dest.of.res <- j
    additional.orig.of.res <- state[[resource.attribute.index]]$data[edge]
    
    if(dest.of.res == additional.orig.of.res){
      if(dest.of.res == orig.of.res){
        return(update)
      }
    }
    
    return(0)
  }

