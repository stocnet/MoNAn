########## gofParamFunctions


# getTieWeights
getTieWeights <- function(cache, dep.var, lvls, ...) {
  m <- cache[[dep.var]]$valuedNetwork
  if (is.null(lvls))
    lvls <- 0:max(m)
  v <- cbind(lvls, 0)
  existing_v <- as.numeric(names(table(m))) + 1
  v[existing_v[existing_v %in% lvls], 2] <-
    table(m)[existing_v %in% lvls]
  v[, 2]
}
