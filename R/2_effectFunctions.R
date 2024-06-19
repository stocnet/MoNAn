########## effectFunctions


#' in_ties_loops
#' 
#' Are individuals that are in locations with a large inflow more likely to stay 
#' in their current location?
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
in_ties_loops <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE) {
    if (getTargetContribution) {
      return((cache[[dep.var]]$valuedNetwork[i, j] * (i == j)) *
        (sum(cache[[dep.var]]$valuedNetwork[, j]) - cache[[dep.var]]$valuedNetwork[j, j]))
    }

    if (i == j) {
      return(update * (sum(cache[[dep.var]]$valuedNetwork[, j]) - cache[[dep.var]]$valuedNetwork[j, j]))
    }

    if (i != j) {
      return(update * cache[[dep.var]]$valuedNetwork[j, j])
    }
  }


#' in_weights_exponent
#' 
#' Is there a preferential attachment in the mobility network, i.e., do individuals 
#' move particularly to popular destinations?
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
#' @keywords internal
in_weights_exponent <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE,
           exponent = .5) {
    if (i == j) {
      return(0)
    }

    # a target contribution is calculated even for unconnected i-j pairs
    if (getTargetContribution) {
      return(sum(cache[[dep.var]]$valuedNetwork[, j])^exponent / length(cache[[dep.var]]$valuedNetwork[, j]))
    }

    v <-
      (sum(cache[[dep.var]]$valuedNetwork[, j]) + update)^exponent -
      (sum(cache[[dep.var]]$valuedNetwork[, j]))^exponent

    return(v)
  }


#' in_weights_GW
#' 
#' Is there a preferential attachment in the mobility network, i.e., do individuals 
#' move particularly to popular destinations?
#' The geometrically weighted version avoids degeneracy.
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
in_weights_GW <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE,
           alpha = 2) {
    if (alpha <= 0) {
      stop("Alpha parameter in in_weights_GW weights function must be positive")
    }
    
    in_weight <- (sum(cache[[dep.var]]$valuedNetwork[,j]))
    
    # a target contribution is calculated even for unconnected i-j pairs
    if (getTargetContribution) {
      g_cum <- function(y, a){
        contr <- 0
        if(y>1){
          for(k in 1:(y-1)){
            contr <- contr + (1 - (1-1/a)^(k))
          }
        }
        return(contr)
      }
      
      return(g_cum(y = in_weight, a = alpha) / (length(cache[[dep.var]]$valuedNetwork[, j])))
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
    
    if(update < 0){
      return(update * g_mar(y = (in_weight + update), a = alpha))
    }
    if(update > 0){
      return(update * g_mar(y = in_weight, a = alpha))
    }
  }

#' present_relations
#' 
#' Do individuals move along many or few paths out of their origin? This models 
#' whether individuals have a tendency against being the only one
#' moving to a particular destination from their origin.
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
present_relations <-
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

    if (getTargetContribution) {
      return((cache[[dep.var]]$valuedNetwork[i, j] > 0) * 1)
    }

    if (cache[[dep.var]]$valuedNetwork[i, j] == 0 &&
      update > 0) {
      return(1)
    }
    if (cache[[dep.var]]$valuedNetwork[i, j] == -update &&
      update < 0) {
      return(-1)
    }

    return(0)
  }

