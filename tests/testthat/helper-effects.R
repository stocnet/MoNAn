# Helpers for testing that change statistics of effects match the difference
# in target statistics before and after a move. Based on target_change_match()
# in R/2x_effects_test.R, but with state, cache and resource covariates as
# arguments instead of the hard-coded lililiState / lililiCache.

# sum of target contributions over all origin-destination pairs
target_stat <- function(effect, args, state, cache) {
  f <- get(effect, mode = "function")
  nodes <- seq_along(state[[state[[state$dep.var]]$nodeSet[1]]]$ids)
  total <- 0
  for (i in nodes) {
    for (j in nodes) {
      total <- total + do.call(f, c(list(
        dep.var = 1, state = state, cache = cache,
        i = i, j = j, edge = 1, update = 1,
        getTargetContribution = TRUE
      ), args))
    }
  }
  total
}

# move resource `edge` from its current destination to `j_new`, the same way
# the sampler does in simulateOneStep(): remove with the old cache, update the
# cache, then add with the updated cache
check_move <- function(effect, args, edge, j_new,
                       state = lililiState, cache = lililiCache,
                       resourceCovariates = "sex") {
  f <- get(effect, mode = "function")
  dv <- state$dep.var
  i <- state[[dv]]$data[edge, 1]
  j <- state[[dv]]$data[edge, 2]

  change <- do.call(f, c(list(
    dep.var = 1, state = state, cache = cache,
    i = i, j = j, edge = edge, update = -1,
    getTargetContribution = FALSE
  ), args))

  updatedCache <- cache
  updatedCache[[1]] <- updateWeightedCache(
    cache[[1]], i, j,
    resourceID = edge, state = state, dep.var = 1, update = -1
  )
  # as in simulateOneStep(), the moving resource's row of the edgelist is NA
  # while the change of adding it is calculated
  removedState <- state
  removedState[[dv]]$data[edge, ] <- NA

  change <- change + do.call(f, c(list(
    dep.var = 1, state = removedState, cache = updatedCache,
    i = i, j = j_new, edge = edge, update = 1,
    getTargetContribution = FALSE
  ), args))

  newState <- state
  newState[[dv]]$data[edge, 2] <- j_new
  newCache <- createInternalCache(newState, resourceCovariates = resourceCovariates)

  c(
    i = i, j = j,
    change = change,
    before = target_stat(effect, args, state, cache),
    after = target_stat(effect, args, newState, newCache)
  )
}

# moves to test: random ones plus moves out of loops, into loops and no-ops.
# With MONAN_FULL_EFFECT_TESTS=true, every resource is moved once.
sample_moves <- function(n = 30, seed = 1, state = lililiState) {
  el <- state[[state$dep.var]]$data
  nNodes <- length(state[[state[[state$dep.var]]$nodeSet[1]]]$ids)
  set.seed(seed)

  if (identical(Sys.getenv("MONAN_FULL_EFFECT_TESTS"), "true")) {
    random <- seq_len(nrow(el))
  } else {
    random <- sample(nrow(el), n)
  }
  loopEdges <- which(el[, 1] == el[, 2])[1:3]
  nonLoopEdges <- which(el[, 1] != el[, 2])
  intoLoop <- nonLoopEdges[1:3]
  noOp <- nonLoopEdges[4:5]

  data.frame(
    edge = c(random, loopEdges, intoLoop, noOp),
    j_new = c(
      sample(nNodes, length(random), replace = TRUE),
      vapply(loopEdges, function(e) sample(setdiff(seq_len(nNodes), el[e, 1]), 1), numeric(1)),
      el[intoLoop, 1],
      el[noOp, 2]
    )
  )
}

expect_change_matches <- function(effect, args = list(), moves = sample_moves(),
                                  state = lililiState, cache = lililiCache,
                                  resourceCovariates = "sex") {
  for (k in seq_len(nrow(moves))) {
    r <- check_move(effect, args, moves$edge[k], moves$j_new[k],
                    state, cache, resourceCovariates)
    # tolerance relative to the size of the target statistic, as changes can
    # be tiny (e.g. AC effects with large weights) while targets are large
    diff <- r[["after"]] - r[["before"]]
    tol <- 1e-8 * max(1, abs(r[["before"]]))
    expect(
      is.finite(r[["change"]]) && abs(r[["change"]] - diff) <= tol,
      sprintf(
        "%s: moving edge %d (%d -> %d) to %d gives change statistic %.10g, but target statistics differ by %.10g",
        effect, moves$edge[k], r[["i"]], r[["j"]], moves$j_new[k], r[["change"]], diff
      )
    )
  }
}
