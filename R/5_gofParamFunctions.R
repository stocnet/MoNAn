########## gofParamFunctions


#' getIndegree
#' 
#' Calculates the weighted indegree distribution of all locations in the 
#' network. 
#' The weighted indegree is simply the column sum of the mobility table.
#' 
#' @param cache Current Cache.
#' @param dep.var Dependent Variable.
#' @param lvls Levels for which the function calculates values. Level k is the
#' number of locations with a weighted indegree of k - 1, i.e., level 1 counts
#' locations with indegree 0. If NULL, all levels from 1 to the maximum
#' indegree + 1 are used.
#' @param ... Additional parameters.
#'
#' @return The degree distribution of simulated networks.
#' @export
#' 
#' @seealso [gofMobilityNetwork()]
#'  
#' @keywords internal
getIndegree <- function(cache, dep.var, lvls, ...) {
  indegrees <- colSums(cache[[dep.var]]$valuedNetwork)
  if (is.null(lvls)) {
    lvls <- 1:(max(indegrees) + 1)
  }
  counts <- tabulate(indegrees + 1, nbins = max(c(lvls, indegrees + 1)))
  counts[lvls]
}


#' getTieWeights
#' 
#' Extracts the distribution of tie weights in the mobility network.
#' 
#' @param cache Current Cache.
#' @param dep.var Dependent Variable.
#' @param lvls Levels for which the function calculates values. Level k is the
#' number of ties with weight k - 1, i.e., level 1 counts the empty dyads.
#' If NULL, all levels from 1 to the maximum tie weight + 1 are used.
#' @param ... Additional parameters.
#'
#' @return The tie weight distribution of simulated networks.
#' @export
#' 
#' @seealso [gofMobilityNetwork()]
#'
#' @keywords internal
getTieWeights <- function(cache, dep.var, lvls, ...) {
  weights <- as.vector(cache[[dep.var]]$valuedNetwork)
  if (is.null(lvls)) {
    lvls <- 1:(max(weights) + 1)
  }
  counts <- tabulate(weights + 1, nbins = max(c(lvls, weights + 1)))
  counts[lvls]
}
