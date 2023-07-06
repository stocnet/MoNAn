########## gofParamFunctions


#' getIndegree
#'
#' @param cache Current Cache
#' @param dep.var Dependent Variable
#' @param lvls Levels for which the function calculates values
#' @param ... Additional parameters
#'
#' @export
#' 
#' @seealso [gofDistributionNetwork()]
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
#' @param cache Current Cache
#' @param dep.var Dependent Variable
#' @param lvls Levels for which the function calculates values
#' @param ... Additional parameters
#'
#' @export
#' 
#' @seealso [gofDistributionNetwork()]
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
