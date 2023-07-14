########## effectFunctions


#' crowding_out_by_resource_inflow
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
#' @keywords internal
crowding_out_by_resource_inflow <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
    # proportion needs binary coding with only 0 and 1
    if (!all(state[[resource.attribute.index]]$data %in% c(0, 1))) {
      stop(
        "effect crowding_out_by_resource_inflow only defined for binary covariates coded 0 1"
      )
    }

    # get the target contribution
    if (getTargetContribution) {
      # if a dyad not on the diagonal is checked, return 0
      if (i != j) {
        return(0)
      }

      # get the diagonal value
      loops <-
        cache[[dep.var]]$valuedNetwork[j, j] - cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]

      # get the number of X'ers that move in to the node

      in.resources <-
        sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) -
        cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]

      return(loops * in.resources)
    }

    # get the change statistics
    # change statistics depends on either a new loop that is changed
    # or on the inflow proportion that the target node has
    # both need to know the proportion before

    inflow.res.before <-
      sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) -
      cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]

    # first if a loop is formed
    if (i == j) {
      return(inflow.res.before * update * (1 - state[[resource.attribute.index]]$data[edge]))
    }

    # now if no loop is formed and the inflow number changes
    loops <-
      cache[[dep.var]]$valuedNetwork[j, j] - cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]

    inflow.res.after <-
      sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) -
      cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j] +
      update * state[[resource.attribute.index]]$data[edge]

    change <- inflow.res.after - inflow.res.before
    return(change * loops)
  }


#' crowding_out_prop_covar_bin
#' 
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
#' @keywords internal
crowding_out_prop_covar_bin <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
    # proportion needs binary coding with only 0 and 1
    if (!all(state[[resource.attribute.index]]$data %in% c(0, 1))) {
      stop("effect crowding_out_prop_covar_bin only defined for binary covariates coded 0 1")
    }

    # get the target contribution
    if (getTargetContribution) {
      # if a dyad not on the diagonal is checked, return 0
      if (i != j) {
        return(0)
      }

      # get the diagonal value
      loops <-
        cache[[dep.var]]$valuedNetwork[i, j] - cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j]

      # get the proportion of X'ers that move in to the node
      # if nobody comes in, the prop equals 0.5
      if ((sum(cache[[dep.var]]$valuedNetwork[, j]) - cache[[dep.var]]$valuedNetwork[i, j]) == 0) {
        prop <- 0.5
      } else {
        prop <-
          (sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) -
            cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j]) /
            (sum(cache[[dep.var]]$valuedNetwork[, j]) - cache[[dep.var]]$valuedNetwork[i, j])
      }

      return(loops * prop)
    }

    # get the change statistics
    # change statistics depends on either a new loop that is changed
    # or on the inflow proportion that the target node has
    # both need to know the proportion before

    if ((sum(cache[[dep.var]]$valuedNetwork[, j]) - cache[[dep.var]]$valuedNetwork[j, j]) == 0) {
      propBefore <- 0.5
    } else {
      propBefore <-
        (sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) -
          cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]) /
          (sum(cache[[dep.var]]$valuedNetwork[, j]) - cache[[dep.var]]$valuedNetwork[j, j])
    }

    # first if a loop is formed
    if (i == j) {
      # if the res is of attribute == 1, then there are no change stats
      if (state[[resource.attribute.index]]$data[edge] == 1) {
        return(0)
      }
      # if the res attribute == 0, the change stat is the proportion
      return(propBefore * update)
    }

    # now if no loop is formed and the proportion changes
    # number of loops of type 0 people in target occ
    loops <-
      cache[[dep.var]]$valuedNetwork[j, j] - cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]

    # if the last one leaves, the proportion after becomes 0

    if (((sum(cache[[dep.var]]$valuedNetwork[, j]) + update) - cache[[dep.var]]$valuedNetwork[j, j]) == 0) {
      propAfter <- 0.5
    } else {
      propAfter <-
        ((sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) +
          update * state[[resource.attribute.index]]$data[edge]) -
          cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]) /
          ((sum(cache[[dep.var]]$valuedNetwork[, j]) + update) - cache[[dep.var]]$valuedNetwork[j, j])
    }

    change <- propAfter - propBefore
    return(change * loops)
  }




#' in_proportion_exponent_covar_bin
#' 
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
#' @param exponent 
#'
#' 
#' @keywords internal
in_proportion_exponent_covar_bin <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           exponent = 2) {
    if (!all(state[[resource.attribute.index]]$data %in% c(0, 1))) {
      stop("effect in_proportion_exponent_covar only defined for binary covariates coded 0 1")
    }

    # a target contribution is calculated even for unconnected i-j pairs
    if (getTargetContribution) {
      if (sum(cache[[dep.var]]$valuedNetwork[, j]) == 0) {
        return(0)
      }
      # calculate proportion that have covar = 1
      prop <-
        (sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) /
          sum(cache[[dep.var]]$valuedNetwork[, j]))
      # divide by number of rows, as target stats goes through each dyad
      return(prop^exponent / length(cache[[dep.var]]$valuedNetwork[, j]))
    }

    # if a node has nobdy in it, then the old Proportion must be 0, not NaN
    if (sum(cache[[dep.var]]$valuedNetwork[, j]) == 0) {
      oldProp <- 0
    } else {
      oldProp <-
        (sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) /
          sum(cache[[dep.var]]$valuedNetwork[, j]))^exponent
    }

    # if a resource is the last to leave a node, then the new proportion must be 0, not NaN
    if ((sum(cache[[dep.var]]$valuedNetwork[, j]) + update) == 0) {
      newProp <- 0
    } else {
      newProp <-
        ((sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) +
          update * state[[resource.attribute.index]]$data[edge]) /
          (sum(cache[[dep.var]]$valuedNetwork[, j]) + update))^
          exponent
    }
    return(newProp - oldProp)
  }


#' in_ties_loops
#' 
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
#' @keywords internal
in_ties_loops <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
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
#' @keywords internal
in_weights_exponent <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
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


#' in_weights_exponent_covar
#' 
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
#' @param exponent 
#'
#' 
#' @keywords internal
in_weights_exponent_covar <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           exponent = .5) {
    if (i == j) {
      return(0)
    }

    # a target contribution is calculated even for unconnected i-j pairs
    if (getTargetContribution) {
      return(sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j])^
        exponent /
        length(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]))
    }

    v <-
      (sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) + update * state[[resource.attribute.index]]$data[edge])^
      exponent -
      (sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]))^
        exponent

    return(v)
  }




#' min_reciprocity
#' 
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
#' @keywords internal
min_reciprocity <-
  function(dep.var = 1,
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
      return(cache[[dep.var]]$minNetwork[i, j] / 2)
    }

    # simplified version that assumes that update are 1 or -1 and the network is an integer network
    if (update > 0 &&
      cache[[dep.var]]$netFlowsNetwork[i, j] < 0) {
      return(1)
    }
    if (update < 0 &&
      cache[[dep.var]]$netFlowsNetwork[i, j] <= 0) {
      return(-1)
    }
    return(0)
  }


#' min_reciprocity_resource_covar
#' 
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
#' @keywords internal
min_reciprocity_resource_covar <-
  function(dep.var = 1,
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
      return(min(
        cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j],
        cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, i]
      ) / 2)
    }

    # simplified version that assumes that update are 1 or -1 and the network is an integer network
    if (update > 0 &&
      (cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j] <
        cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, i])) {
      return(1)
    }
    if (update < 0 &&
      (cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j] <=
        cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, i])) {
      return(-1)
    }
    return(0)
  }


#' min_transitivity
#' 
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
#' @keywords internal
min_transitivity <-
  function(dep.var = 1,
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

    m.min <- cache[[dep.var]]$minNetwork
    diag(m.min) <- 0

    if (getTargetContribution) {
      twoPaths <- apply(cbind(m.min[i, ], m.min[, j]), 1, min)
      triadValues <- mapply(min, m.min[i, j], twoPaths)
      v <- sum(triadValues) / 6
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
    twoPathValues <- apply(cbind(m.min[i, ], m.min[, j]), 1, min)
    if (update > 0) {
      v <- sum(dyadValue < twoPathValues)
    }
    if (update < 0) {
      v <- -sum(dyadValue <= twoPathValues)
    }
    return(v)
  }


#' netflow_transitivity
#' 
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
#' @keywords internal
netflow_transitivity <-
  function(dep.var = 1,
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

    getPathWeights <- function(ik, kj) {
      if (ik < 0 && kj < 0) {
        return(0)
      }
      return(min(abs(ik), abs(kj)))
    }

    if (getTargetContribution) {
      if (cache[[dep.var]]$netFlowsNetwork[i, j] <= 0) {
        return(0)
      }
      pathWeights <-
        mapply(
          getPathWeights,
          cache[[dep.var]]$netFlowsNetwork[i, ],
          cache[[dep.var]]$netFlowsNetwork[, j]
        )
      v <-
        sum(mapply(min, pathWeights, cache[[dep.var]]$netFlowsNetwork[i, j]))
      return(v / 3)
    }

    # calculate contribution before change and relevant path weights (those relating to transitive structures)
    if (cache[[dep.var]]$netFlowsNetwork[i, j] > 0) {
      pathWeights <-
        mapply(
          getPathWeights,
          cache[[dep.var]]$netFlowsNetwork[i, ],
          cache[[dep.var]]$netFlowsNetwork[, j]
        )
      contributionBefore <-
        sum(mapply(min, pathWeights, cache[[dep.var]]$netFlowsNetwork[i, j]))
    }
    if (cache[[dep.var]]$netFlowsNetwork[i, j] < 0) {
      pathWeights <-
        mapply(
          getPathWeights,
          cache[[dep.var]]$netFlowsNetwork[j, ],
          cache[[dep.var]]$netFlowsNetwork[, i]
        )
      contributionBefore <-
        sum(mapply(min, pathWeights, cache[[dep.var]]$netFlowsNetwork[j, i]))
    }
    if (cache[[dep.var]]$netFlowsNetwork[i, j] == 0) {
      contributionBefore <- 0
      if (update > 0) {
        pathWeights <-
          mapply(
            getPathWeights,
            cache[[dep.var]]$netFlowsNetwork[i, ],
            cache[[dep.var]]$netFlowsNetwork[, j]
          )
      }
      if (update < 0) {
        pathWeights <-
          mapply(
            getPathWeights,
            cache[[dep.var]]$netFlowsNetwork[j, ],
            cache[[dep.var]]$netFlowsNetwork[, i]
          )
      }
    }

    # calculate contribution after the tie change
    tieWeightAfter <- cache[[dep.var]]$netFlowsNetwork[i, j] + update
    contributionAfter <-
      sum(mapply(min, pathWeights, abs(tieWeightAfter)))

    return(contributionAfter - contributionBefore)
  }


#' present_relations
#' 
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
#' @keywords internal
present_relations <-
  function(dep.var = 1,
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




#' staying_by_prop_bin_inflow
#' 
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
#' @keywords internal
staying_by_prop_bin_inflow <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
    # proportion needs binary coding with only 0 and 1
    if (!all(state[[resource.attribute.index]]$data %in% c(0, 1))) {
      stop("effect staying_by_prop_bin_inflow only defined for binary covariates coded 0 1")
    }

    # get the target contribution
    if (getTargetContribution) {
      # if a dyad not on the diagonal is checked, return 0
      if (i != j) {
        return(0)
      }

      # get the diagonal value, i.e. the amount of people that stay here
      loops <- cache[[dep.var]]$valuedNetwork[i, j]

      # get the proportion of X'ers that move in to the node
      # if nobody comes in, the prop equals 0.5
      if ((sum(cache[[dep.var]]$valuedNetwork[, j]) - cache[[dep.var]]$valuedNetwork[j, j]) == 0) {
        prop <- 0.5
      } else {
        prop <-
          (sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) -
            cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]) /
            (sum(cache[[dep.var]]$valuedNetwork[, j]) - cache[[dep.var]]$valuedNetwork[j, j])
      }

      return(loops * prop)
    }

    # get the change statistics
    # change statistics depends on either a new loop that is changed
    # or on the inflow proportion that the target node has
    # both need to know the proportion before

    if ((sum(cache[[dep.var]]$valuedNetwork[, j]) - cache[[dep.var]]$valuedNetwork[j, j]) == 0) {
      propBefore <- 0.5
    } else {
      propBefore <-
        (sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) -
          cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]) /
          (sum(cache[[dep.var]]$valuedNetwork[, j]) - cache[[dep.var]]$valuedNetwork[j, j])
    }

    # first if a loop is formed
    if (i == j) {
      # the change stat is the proportion of res with attribute that move in
      return(propBefore * update)
    }

    # now if no loop is formed and the proportion changes
    # number of loops in target occ
    loops <- cache[[dep.var]]$valuedNetwork[j, j]

    # if the last one leaves, the proportion after becomes 0.5

    if (((sum(cache[[dep.var]]$valuedNetwork[, j]) + update) - cache[[dep.var]]$valuedNetwork[j, j]) == 0) {
      propAfter <- 0.5
    } else {
      propAfter <-
        ((sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) +
          update * state[[resource.attribute.index]]$data[edge]) -
          cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]) /
          ((sum(cache[[dep.var]]$valuedNetwork[, j]) + update) - cache[[dep.var]]$valuedNetwork[j, j])
    }

    change <- propAfter - propBefore
    return(change * loops)
  }


#' staying_by_resource_inflow
#' 
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
#' @keywords internal
staying_by_resource_inflow <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
    # proportion needs binary coding with only 0 and 1
    if (!all(state[[resource.attribute.index]]$data %in% c(0, 1))) {
      stop("effect staying_by_resource_inflow only defined for binary covariates coded 0 1")
    }

    # get the target contribution
    if (getTargetContribution) {
      # if a dyad not on the diagonal is checked, return 0
      if (i != j) {
        return(0)
      }

      # get the diagonal value
      loops <- cache[[dep.var]]$valuedNetwork[i, j]

      # get the number of X'ers that move in to the node

      in.resources <-
        sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) -
        cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]

      return(loops * in.resources)
    }

    # get the change statistics
    # change statistics depends on either a new loop that is changed
    # or on the inflow proportion that the target node has
    # both need to know the proportion before

    inflow.res.before <-
      sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) -
      cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j]

    # first if a loop is formed
    if (i == j) {
      return(inflow.res.before * update)
    }

    # now if no loop is formed and the inflow number changes
    loops <- cache[[dep.var]]$valuedNetwork[j, j]

    inflow.res.after <-
      sum(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][, j]) -
      cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][j, j] +
      update * state[[resource.attribute.index]]$data[edge]

    change <- inflow.res.after - inflow.res.before
    return(change * loops)
  }


