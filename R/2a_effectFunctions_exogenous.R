########## effectFunctions: exogenous effects

#' alter_covariate
#'
#' @param dep.var 
#' @param attribute.index 
#' @param state 
#' @param cache 
#' @param i 
#' @param j 
#' @param edge 
#' @param update 
#' @param loop.contribution 
#' @param getTargetContribution 
#'
#' 
#' @return None.
#' @keywords internal
alter_covariate <-
  function(dep.var = 1,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           loop.contribution = T,
           getTargetContribution = F) {
    if (!loop.contribution) {
      if (i == j) {
        return(0)
      }
    }
    
    if (getTargetContribution) {
      return(state[[attribute.index]]$data[j] * cache[[dep.var]]$valuedNetwork[i, j])
    }
    
    return(state[[attribute.index]]$data[j] * update)
  }


#' dyadic_covariate
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
dyadic_covariate <-
  function(dep.var = 1,
           attribute.index,
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
      return(cache[[dep.var]]$valuedNetwork[i, j] * state[[attribute.index]]$data[i, j])
    }
    
    return(update * state[[attribute.index]]$data[i, j])
  }


#' dyadic_covariate_resource_attribute
#' 
#'
#' @param dep.var 
#' @param attribute.index 
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
dyadic_covariate_resource_attribute <-
  function(dep.var = 1,
           attribute.index,
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
      return(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j] * state[[attribute.index]]$data[i, j])
    }
    
    return(update * state[[resource.attribute.index]]$data[edge] * state[[attribute.index]]$data[i, j])
  }

#' resource_covar_to_node_covar
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
#' @param loop.contribution 
#' @param getTargetContribution 
#'
#' 
#' @return None.
#' @keywords internal
resource_covar_to_node_covar <-
  function(dep.var = 1,
           resource.attribute.index,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           loop.contribution = F,
           getTargetContribution = F) {
    if (!loop.contribution) {
      if (i == j) {
        return(0)
      }
    }
    
    if (getTargetContribution) {
      return(state[[attribute.index]]$data[j] *
               cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j])
    }
    
    return(update * state[[resource.attribute.index]]$data[edge] * state[[attribute.index]]$data[j])
  }


#' same_covariate
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
same_covariate <-
  function(dep.var = 1,
           attribute.index,
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
    
    if (is.null(state[[attribute.index]]$same)) {
      stop("Effect same_covariate expects covariates with a 'same' attribute.")
    }
    
    if (getTargetContribution) {
      return(cache[[dep.var]]$valuedNetwork[i, j] * state[[attribute.index]]$same[i, j])
    }
    
    return(update * state[[attribute.index]]$same[i, j])
  }


#' sim_covariate
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
sim_covariate <-
  function(dep.var = 1,
           attribute.index,
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
    
    if (is.null(state[[attribute.index]]$sim)) {
      stop("Effect sim_covariate expects covariates with a 'sim' attribute.")
    }
    
    if (getTargetContribution) {
      return(cache[[dep.var]]$valuedNetwork[i, j] * state[[attribute.index]]$sim[i, j])
    }
    
    return(update * state[[attribute.index]]$sim[i, j])
  }

