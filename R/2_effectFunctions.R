########## effectFunctions


#' crowding_out_prop_covar_bin
#' 
#' Is the tendency to stay in vs. move out of a location of individuals of type 
#' non-w dependent on the proportion of individuals of type w moving into the location? 
#' This is especially geared towards modelling how some locations become more or 
#' less attractive dependent on the change in composition for
#' particular groups. This models segregation dynamics.
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
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
           getTargetContribution = FALSE) {
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




#' staying_by_prop_bin_inflow
#' 
#' Is the tendency to stay in vs. move out of a location dependent on the proportion 
#' of individuals of type w that enter the location? This is especially geared 
#' towards modelling how some locations become more or less
#' attractive dependent on the change in composition.
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
#' @return Returns the change statistic or target statistic of the effect for 
#' internal use by the estimation algorithm.
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
           getTargetContribution = FALSE) {
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


