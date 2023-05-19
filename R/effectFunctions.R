########## effectFunctions


# alter_covariate
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


# binary_transitivity
binary_transitivity <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           k = NULL,
           l = NULL,
           getTargetContribution = F) {
    if (getTargetContribution) {
      return(cache[[dep.var]]$twoPaths[i, j])
    }

    # TODO check for overlap between indexes
    # TODO also add the number of two-in starts and two-out starts to the change contribution
    #      to align it with the network statistic value
    change <- 0
    change <- change - cache[[dep.var]]$twoPaths[i, j]
    change <- change + cache[[dep.var]]$twoPaths[k, l]
    # subtract a potentially broken two-paths due to the removal of the first tie
    isPathBroken <-
      ((i == k && l %in% cache[[dep.var]]$outNeighbors[[j]]) ||
        j == l && (k %in% cache[[dep.var]]$outNeighbors[[j]]))
    change <- change + isPathBroken
    change
  }


# crowding_out_by_resource_inflow
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


# crowding_out_prop_covar_bin
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


# dyadic_covariate
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


# dyadic_covariate_resource_attribute
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


# dyadic_covariate_tie_weights_sigmoid
dyadic_covariate_tie_weights_sigmoid <-
  function(dep.var = 1,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           alpha) {
    if (i == j) {
      return(0)
    }
    if (alpha <= 0) {
      stop("Alpha parameter in sigmoid tie weights function must be positive")
    }

    if (getTargetContribution) {
      nResources <- cache[[dep.var]]$valuedNetwork[i, j]
      v <- 0
      for (zz in 0:nResources) {
        v <- v + zz / (zz + alpha)
      }
      return(v * state[[attribute.index]]$data[i, j])
    }

    value.old <- cache[[dep.var]]$valuedNetwork[i, j]
    value.new <- cache[[dep.var]]$valuedNetwork[i, j] + update
    if (update > 0) {
      v <- value.new / (value.new + alpha)
    }
    if (update < 0) {
      v <- -value.old / (value.old + alpha)
    }

    return(v * state[[attribute.index]]$data[i, j])
  }


# in_proportion_exponent_covar_bin
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


# in_ties_loops
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


# in_weights_exponent
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


# in_weights_exponent_covar
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


# loops
loops <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
    if (getTargetContribution) {
      return(cache[[dep.var]]$valuedNetwork[i, j] * (i == j))
    }

    if (i == j) {
      return(update)
    }
    return(0)
  }


# loops_GW
loops_GW <- function(dep.var = 1,
                     state,
                     cache,
                     i,
                     j,
                     edge,
                     update,
                     getTargetContribution = F,
                     alpha = 2) {
  if (alpha <= 0) {
    stop("Alpha parameter in GW loops weights function must be positive")
  }

  if (getTargetContribution) {
    if (i == j) {
      nResources <- cache[[dep.var]]$valuedNetwork[i, i]

      if (nResources == 0) {
        return(0)
      }

      v <- list()
      for (turn in 1:nResources) {
        ind_cont <- 0
        for (k in 1:turn) {
          ind_cont <- ind_cont + (1 / (alpha)^(k))
        }
        v[[turn]] <- ind_cont
      }
      return(sum(unlist(v)))
    }

    return(0)
  }

  if (i != j) {
    return(0)
  }

  if (i == j) {
    value.old <- cache[[dep.var]]$valuedNetwork[i, i]
    value.new <- cache[[dep.var]]$valuedNetwork[i, i] + update

    if (update > 0) {
      ind_cont <- 0
      for (k in 1:value.new) {
        ind_cont <- ind_cont + (1 / (alpha)^(k))
      }
      v <- ind_cont
    }

    if (update < 0) {
      ind_cont <- 0
      for (k in 1:value.old) {
        ind_cont <- ind_cont + (1 / (alpha)^(k))
      }
      v <- -ind_cont
    }

    return(v)
  }
}


# loops_node_covar
loops_node_covar <-
  function(dep.var = 1,
           attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
    if (i != j) {
      return(0)
    }

    if (getTargetContribution) {
      return(cache[[dep.var]]$valuedNetwork[i, j] * state[[attribute.index]]$data[i])
    }

    return(update * state[[attribute.index]]$data[i])
  }


# loops_resource_covar_node_covar
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
           getTargetContribution = F) {
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


# loops_resource_covar
loops_resource_covar <-
  function(dep.var = 1,
           resource.attribute.index,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F) {
    if (i != j) {
      return(0)
    }

    if (getTargetContribution) {
      return(cache[[dep.var]]$resourceNetworks[[resource.attribute.index]][i, j])
    }

    return(update * state[[resource.attribute.index]]$data[edge])
  }


# loops_weight_sigmoid
loops_weight_sigmoid <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           alpha) {
    if (alpha <= 0) {
      stop("Alpha parameter in sigmoid loops weights function must be positive")
    }

    if (getTargetContribution) {
      if (i == j) {
        nResources <- cache[[dep.var]]$valuedNetwork[i, i]
        v <- 0
        for (zz in 0:nResources) {
          v <- v + zz / (zz + alpha)
        }
        return(v)
      }

      return(0)
    }

    if (i != j) {
      return(0)
    }

    if (i == j) {
      value.old <- cache[[dep.var]]$valuedNetwork[i, i]
      value.new <- cache[[dep.var]]$valuedNetwork[i, i] + update

      if (update > 0) {
        v <- value.new / (value.new + alpha)
      }
      if (update < 0) {
        v <- -value.old / (value.old + alpha)
      }

      return(v)
    }
  }


# loopsNegExp
loopsNegExp <-
  loopsNegExp <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           half.cont = 0.1) {
    if (half.cont <= 0) {
      stop("half.cont parameter in GW loopsNegExp function must be positive")
    }

    getExpChange <- function(count, sum) {
      alpha <- -(1 / half.cont) * log(0.5)
      return(1 - exp(-alpha * count / sum))
    }

    if (getTargetContribution) {
      if (i == j) {
        nLoops <- cache[[dep.var]]$valuedNetwork[i, i]
        if (nLoops == 0) {
          return(0)
        }

        resInRow <- sum(cache[[dep.var]]$valuedNetwork[i, ])
        return(sum(getExpChange(1:nLoops, resInRow)))
      }
      return(0)
    }

    if (i != j) {
      return(0)
    }

    if (i == j) {
      nLoops <- cache[[dep.var]]$valuedNetwork[i, i]
      resInRow <- sum(cache[[dep.var]]$valuedNetwork[i, ])

      if (update == 1) {
        return(getExpChange((nLoops + 1), resInRow))
      }
      if (update == -1) {
        return(-getExpChange((nLoops), resInRow))
      }
    }
  }


# min_reciprocity
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


# min_reciprocity_resource_covar
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


# min_transitivity
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


# netflow_transitivity
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


# present_relations
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


# resource_covar_to_node_covar
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


# same_covariate
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


# sim_covariate
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


# staying_by_prop_bin_inflow
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


# staying_by_resource_inflow
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


# tie_weights_exponent
tie_weights_exponent <-
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

    if (getTargetContribution) {
      return((cache[[dep.var]]$valuedNetwork[i, j])^exponent)
    }

    v <- (cache[[dep.var]]$valuedNetwork[i, j] + update)^exponent -
      (cache[[dep.var]]$valuedNetwork[i, j])^exponent

    return(v)
  }


# tie_weights_sigmoid
tie_weights_sigmoid <-
  function(dep.var = 1,
           state,
           cache,
           i,
           j,
           edge,
           update,
           getTargetContribution = F,
           alpha) {
    if (i == j) {
      return(0)
    }
    if (alpha <= 0) {
      stop("Alpha parameter in sigmoid tie weights function must be positive")
    }

    if (getTargetContribution) {
      nResources <- cache[[dep.var]]$valuedNetwork[i, j]
      v <- 0
      for (zz in 0:nResources) {
        v <- v + zz / (zz + alpha)
      }
      return(v)
    }

    value.old <- cache[[dep.var]]$valuedNetwork[i, j]
    value.new <- cache[[dep.var]]$valuedNetwork[i, j] + update
    if (update > 0) {
      v <- value.new / (value.new + alpha)
    }
    if (update < 0) {
      v <- -value.old / (value.old + alpha)
    }

    return(v)
  }
