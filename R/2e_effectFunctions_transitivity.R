########## effectFunctions: effects concerning transitivity (Neighbourhood dependence)


#' transitivity_basic
#' 
#' Is mobility clustered in groups? This is represented by the count of
#' transitive triads among three nodes. This effect is prone to degeneracy.
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
transitivity_basic <-
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

    # rows and columns of the network without loops
    net <- cache[[dep.var]]$valuedNetwork
    ri <- net[i, ]
    ri[i] <- 0
    rj <- net[j, ]
    rj[j] <- 0
    ci <- net[, i]
    ci[i] <- 0
    cj <- net[, j]
    cj[j] <- 0

    if (getTargetContribution) {
      return(net[i, j] * sum(ri * cj))
    }

    # two-paths, in-stars and out-stars closed by the i-j tie
    ind_cont <- update * (sum(ri * cj) + sum(ri * rj) + sum(ci * cj))

    return(ind_cont)
  }

#' transitivity_min
#' 
#' Is mobility clustered in groups? This is represented by the minimum of reciprocated 
#' mobility being present among three nodes. Using the minimum ensures that the effect 
#' is not degenerate and it is sample size consistent.
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
transitivity_min <-
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

    m.min <- cache[[dep.var]]$minNetwork
    # row i and column j of the min network without loops
    mi <- m.min[i, ]
    mi[i] <- 0
    mj <- m.min[, j]
    mj[j] <- 0

    if (getTargetContribution) {
      twoPaths <- pmin(mi, mj)
      v <- sum(pmin(m.min[i, j], twoPaths)) / 6
      return(v)
    }

    # cases without contributions
    if (update == 0) {
      return(0)
    }

    # change contributions are only possible if the minimum tie in the dyad is changed
    if (update > 0 &&
        cache[[dep.var]]$netFlowsNetwork[i, j] >= 0) {
      return(0)
    }
    if (update < 0 &&
        cache[[dep.var]]$netFlowsNetwork[i, j] > 0) {
      return(0)
    }

    # (+/-)1 change contributions only for those triads where the current tie is the weakest / unique weakest for negative / positive updates
    dyadValue <- m.min[i, j]
    twoPathValues <- pmin(mi, mj)
    if (update > 0) {
      v <- sum(dyadValue < twoPathValues)
    }
    if (update < 0) {
      v <- -sum(dyadValue <= twoPathValues)
    }
    return(v)
  }

#' transitivity_netflow
#' 
#' Do individuals move in one direction in locally ordered triads? E.g., is there 
#' a local hierarchy that individuals follow when moving between locations? The 
#' effect is sample size consistent.
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
transitivity_netflow <-
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

    nf <- cache[[dep.var]]$netFlowsNetwork

    # weight of the paths via each k: the smaller absolute net flow, unless
    # both net flows are negative
    getPathWeights <- function(ik, kj) {
      w <- pmin(abs(ik), abs(kj))
      w[ik < 0 & kj < 0] <- 0
      w
    }

    if (getTargetContribution) {
      if (nf[i, j] <= 0) {
        return(0)
      }
      pathWeights <- getPathWeights(nf[i, ], nf[, j])
      v <- sum(pmin(pathWeights, nf[i, j]))
      return(v / 3)
    }

    # calculate contribution before change and relevant path weights (those relating to transitive structures)
    if (nf[i, j] > 0) {
      pathWeights <- getPathWeights(nf[i, ], nf[, j])
      contributionBefore <- sum(pmin(pathWeights, nf[i, j]))
    }
    if (nf[i, j] < 0) {
      pathWeights <- getPathWeights(nf[j, ], nf[, i])
      contributionBefore <- sum(pmin(pathWeights, nf[j, i]))
    }
    if (nf[i, j] == 0) {
      contributionBefore <- 0
      if (update > 0) {
        pathWeights <- getPathWeights(nf[i, ], nf[, j])
      }
      if (update < 0) {
        pathWeights <- getPathWeights(nf[j, ], nf[, i])
      }
    }

    # calculate contribution after the tie change
    tieWeightAfter <- nf[i, j] + update
    contributionAfter <- sum(pmin(pathWeights, abs(tieWeightAfter)))

    return(contributionAfter - contributionBefore)
  }

#' transitivity_AC
#'
#' Is mobility clustered in groups? This is represented by the count of
#' transitive triads among three nodes, where the number of two-paths is
#' geometrically weighted down to avoid degeneracy.
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
transitivity_AC <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE,
           alpha = 1.1) {

    if (i == j) {
      return(0)
    }

    # vectorised; for y = 0 the contribution is 0
    g_mar <- function(y, a) {
      1 - (1 - 1 / a)^y
    }

    # rows and columns of the network without loops; the network itself is
    # not copied
    net <- cache[[dep.var]]$valuedNetwork
    # two-paths between all pairs of locations, kept up to date in the cache
    # (see getTwoPathNetwork); they replace the products of the network with
    # its rows and columns, ri %*% net and net %*% cj, which gave the same
    # (exact, integer) values
    twoPaths <- getTwoPathNetwork(cache[[dep.var]])
    ri <- net[i, ]
    ri[i] <- 0
    rj <- net[j, ]
    rj[j] <- 0
    ci <- net[, i]
    ci[i] <- 0
    cj <- net[, j]
    cj[j] <- 0

    if (getTargetContribution) {
      # two-paths from i to j: sum(ri * cj)
      return(ri[j] * g_mar(twoPaths[i, j], alpha))
    }

    ### change contribution in three parts

    # part 1: i,j closes two-path

    ind_cont <- update * g_mar(twoPaths[i, j], alpha)

    # part 2: i,j closes instar
    # two-paths from i to every k (without loops; row i of the two-path
    # network); adding update to the i-j tie adds update * (j's row) to them
    ri_u <- ri
    ri_u[j] <- ri_u[j] + update

    nTP_before <- twoPaths[i, ]
    con_be <- sum(ri * g_mar(nTP_before, alpha)) - (ri[j] * g_mar(nTP_before[j], alpha))
    nTP_after <- nTP_before + update * rj
    con_af <- sum(ri_u * g_mar(nTP_after, alpha)) - (ri_u[j] * g_mar(nTP_after[j], alpha))
    ind_cont <- ind_cont + (con_af - con_be)

    # part 3: i,j closes outstar
    # two-paths from every k to j (without loops; column j of the two-path
    # network); adding update to the i-j tie adds update * (i's column) to them

    nTP_before <- twoPaths[, j]
    con_be <- sum(cj * g_mar(nTP_before, alpha))
    nTP_after <- nTP_before + update * ci
    con_af <- sum(cj * g_mar(nTP_after, alpha))
    ind_cont <- ind_cont + (con_af - con_be)

    return(ind_cont)
  }

#' transitivity_GW
#'
#' Is mobility clustered in groups? This is represented by the count of
#' transitive triads among three nodes, where the number of two-paths is
#' geometrically weighted down to avoid degeneracy.
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
transitivity_GW <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = FALSE,
           alpha = 1.1) {

    if (i == j) {
      return(0)
    }

    # sum_{k = 1}^{y} a^-k, vectorised
    g_mar <- function(y, a) {
      if (a == 1) {
        return(y)
      }
      (1 - a^(-y)) / (a - 1)
    }

    # rows and columns of the network without loops; the network itself is
    # not copied
    net <- cache[[dep.var]]$valuedNetwork
    # two-paths between all pairs of locations, kept up to date in the cache
    # (see getTwoPathNetwork); they replace the products of the network with
    # its rows and columns, ri %*% net and net %*% cj, which gave the same
    # (exact, integer) values
    twoPaths <- getTwoPathNetwork(cache[[dep.var]])
    ri <- net[i, ]
    ri[i] <- 0
    rj <- net[j, ]
    rj[j] <- 0
    ci <- net[, i]
    ci[i] <- 0
    cj <- net[, j]
    cj[j] <- 0

    if (getTargetContribution) {
      # two-paths from i to j: sum(ri * cj)
      return(ri[j] * g_mar(twoPaths[i, j], alpha))
    }

    ### change contribution in three parts

    # part 1: i,j closes two-path

    ind_cont <- update * g_mar(twoPaths[i, j], alpha)

    # part 2: i,j closes instar
    # two-paths from i to every k (without loops; row i of the two-path
    # network); adding update to the i-j tie adds update * (j's row) to them
    ri_u <- ri
    ri_u[j] <- ri_u[j] + update

    nTP_before <- twoPaths[i, ]
    con_be <- sum(ri * g_mar(nTP_before, alpha)) - (ri[j] * g_mar(nTP_before[j], alpha))
    nTP_after <- nTP_before + update * rj
    con_af <- sum(ri_u * g_mar(nTP_after, alpha)) - (ri_u[j] * g_mar(nTP_after[j], alpha))
    ind_cont <- ind_cont + (con_af - con_be)

    # part 3: i,j closes outstar
    # two-paths from every k to j (without loops; column j of the two-path
    # network); adding update to the i-j tie adds update * (i's column) to them

    nTP_before <- twoPaths[, j]
    con_be <- sum(cj * g_mar(nTP_before, alpha))
    nTP_after <- nTP_before + update * ci
    con_af <- sum(cj * g_mar(nTP_after, alpha))
    ind_cont <- ind_cont + (con_af - con_be)

    return(ind_cont)
  }

#' triad120D
#' 
#' Models the prevalence of the 120D triad
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
#'
#' @keywords internal
triad120D <- function(dep.var = 1,
                      state,
                      cache,
                      i,
                      j,
                      edge,
                      update,
                      getTargetContribution = FALSE){
  if (i == j) {
    return(0)
  }

  # rows and columns of the network without loops
  net <- cache[[dep.var]]$valuedNetwork
  ri <- net[i, ]
  ri[i] <- 0
  rj <- net[j, ]
  rj[j] <- 0
  ci <- net[, i]
  ci[i] <- 0
  cj <- net[, j]
  cj[j] <- 0

  # the target contribution is calculated for the i-j pair that is reciprocated
  if (getTargetContribution) {
    n_down_stars <- sum(ci * cj)
    return(net[i, j] * net[j, i] * n_down_stars / 2)
  }

  # for the change stat, four contributions are needed, one for each tie var,
  # but two are isomorphic

  # tie on the reciprocated dyad
  ind_cont <- update * sum(ci * cj) * net[j, i]

  # tie on the down star (reciprocated ties of j: cj * rj)
  ind_cont <- ind_cont + update * sum(ri * cj * rj)

  return(ind_cont)
}

#' triad120U
#' 
#' Models the prevalence of the 120U triad
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
#'
#' @keywords internal
triad120U <- function(dep.var = 1,
                      state,
                      cache,
                      i,
                      j,
                      edge,
                      update,
                      getTargetContribution = FALSE){
  if (i == j) {
    return(0)
  }

  # rows and columns of the network without loops
  net <- cache[[dep.var]]$valuedNetwork
  ri <- net[i, ]
  ri[i] <- 0
  rj <- net[j, ]
  rj[j] <- 0
  ci <- net[, i]
  ci[i] <- 0
  cj <- net[, j]
  cj[j] <- 0

  # the target contribution is calculated for the i-j pair that is reciprocated
  if (getTargetContribution) {
    n_up_stars <- sum(ri * rj)
    return(net[i, j] * net[j, i] * n_up_stars / 2)
  }

  # for the change stat, four contributions are needed, one for each tie var,
  # but two are isomorphic

  # tie on the reciprocated dyad
  ind_cont <- update * sum(ri * rj) * net[j, i]

  # tie on the up star (reciprocated ties of i: ri * ci)
  ind_cont <- ind_cont + update * sum(cj * ri * ci)

  return(ind_cont)
}

#' triad120C
#' 
#' Models the prevalence of the 120C triad
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
#'
#' @keywords internal
triad120C <- function(dep.var = 1,
                      state,
                      cache,
                      i,
                      j,
                      edge,
                      update,
                      getTargetContribution = FALSE){
  if (i == j) {
    return(0)
  }

  # rows and columns of the network without loops
  net <- cache[[dep.var]]$valuedNetwork
  ri <- net[i, ]
  ri[i] <- 0
  rj <- net[j, ]
  rj[j] <- 0
  ci <- net[, i]
  ci[i] <- 0
  cj <- net[, j]
  cj[j] <- 0

  # the target contribution is calculated for the i-j pair that is reciprocated
  if (getTargetContribution) {
    n_two_paths <- sum(ri * cj)
    return(net[i, j] * net[j, i] * n_two_paths)
  }

  # for the change stat, four contributions are needed, one for each tie var

  # tie on the reciprocated dyad; transitive closure
  ind_cont <- update * sum(ri * cj) * net[j, i]

  # tie on the reciprocated dyad; cyclic closure
  ind_cont <- ind_cont + update * sum(ci * rj) * net[j, i]

  # tie on the first part of two-path (reciprocated ties of i: ri * ci)
  ind_cont <- ind_cont + update * sum(rj * ri * ci)

  # tie on the second part of two-path (reciprocated ties of j: cj * rj)
  ind_cont <- ind_cont + update * sum(ci * cj * rj)

  return(ind_cont)
}
