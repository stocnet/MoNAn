########## effectFunctions: exogenous effects

#' alter_covariate
#' 
#' Are locations higher on some attribute v more popular targets of mobility? 
#' E.g., do workers have a tendency to move to larger organisations?
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
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
#' Is mobility between locations predicted by the dyadic covariate u? E.g., is mobility 
#' likely between organisations that are in the same region? Note that in many cases 
#' dyadic covariates can convey the same information as
#' the ‘same covariate’ or the ‘covariate similarity’ effects.
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
#' Is mobility between locations predicted by the dyadic covariate u weighted by 
#' the individual covariate w? E.g., is mobility of women more likely between 
#' organisations that are in the same region? Note that this effect can be used 
#' to also model the interaction between the ‘same covariate’/‘covariate similarity’ 
#' effect and individual attributes, since sameness and similarity between locations 
#' can be translated into dyadic covariates.
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
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
#' Do individuals with some individual attribute w tend to move to locations with 
#' some location characteristic v? E.g., do women move to larger organisations than men?
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
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
#' Is mobility more likely between locations that are identical on some attribute v? 
#' E.g., is mobility more likely between organisations that are located in the same region?
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
#' Is mobility more likely between locations that are similar on some attribute v? 
#' E.g., is mobility more likely between organisations that have a similar size?
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

