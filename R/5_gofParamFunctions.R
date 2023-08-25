########## gofParamFunctions


#' getIndegree
#' 
#' Calculates the weighted indegree distribution of all locations in the 
#' network. 
#' The weighted indegree is simply the column sum of the mobility table.
#' 
#' @param cache Current Cache.
#' @param dep.var Dependent Variable.
#' @param lvls Levels for which the function calculates values.
#' @param ... Additional parameters.
#'
#' @return The degree distribution of simulated networks.
#' @export
#' 
#' @seealso [gofMobilityNetwork()]
#'  
#' @keywords internal
getIndegree <- function(cache, dep.var, lvls, ...) {
  m <- cache[[dep.var]]$valuedNetwork
  v <- cbind(lvls, 0)
  existing_v <- as.numeric(names(table(colSums(m)))) + 1
  v[existing_v[existing_v %in% lvls], 2] <- table(colSums(m))[existing_v %in% lvls]
  v[, 2]
}


#' getTieWeights
#' 
#' Extracts the distribution of tie weights in the mobility network.
#' 
#' @param cache Current Cache.
#' @param dep.var Dependent Variable.
#' @param lvls Levels for which the function calculates values.
#' @param ... Additional parameters.
#'
#' @return The tie weight distribution of simulated networks.
#' @export
#' 
#' @seealso [gofMobilityNetwork()]
#'
#' @keywords internal
getTieWeights <- function(cache, dep.var, lvls, ...) {
  m <- cache[[dep.var]]$valuedNetwork
  if (is.null(lvls)) {
    lvls <- 0:max(m)
  }
  v <- cbind(lvls, 0)
  existing_v <- as.numeric(names(table(m))) + 1
  v[existing_v[existing_v %in% lvls], 2] <-
    table(m)[existing_v %in% lvls]
  v[, 2]
}
