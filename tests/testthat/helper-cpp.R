# Helpers for testing the C++ sampler and effects against their R versions

# the C++ and R samplers give the same chain for the same seed: same edgelist,
# cache and random number state afterwards
compare_chains <- function(effects, parameters, n = 5000, allowLoops = TRUE,
                           multinomial = FALSE, state = myState, seed = 42) {
  dv <- state$dep.var
  oldOptions <- options(MoNAn2.useCpp = FALSE)
  on.exit(options(oldOptions), add = TRUE)
  set.seed(seed)
  r <- simulateNSteps(dv, state, state$cache, effects, parameters,
                      allowLoops = allowLoops, multinomialProposal = multinomial, n = n)
  rNext <- runif(1)
  options(MoNAn2.useCpp = TRUE)
  set.seed(seed)
  cpp <- simulateChain(dv, state, state$cache, effects, parameters,
                       burnIn = n, iterations = 0, thinning = 0,
                       allowLoops = allowLoops, multinomialProposal = multinomial)
  cppNext <- runif(1)

  expect_identical(cpp$state[[dv]]$data, toListState(r$state, dv)[[dv]]$data)
  expect_equal(cpp$cache, toListCache(r$cache), tolerance = 0)
  expect_identical(cppNext, rNext)
}

# every effect of `effects` gives exactly the same change statistic and target
# contribution in C++ as in R, for random moves including moves into and out
# of loops
expect_cpp_matches_r <- function(effects, state = myState, n = 100, seed = 3) {
  dv <- state$dep.var
  el <- state[[dv]]$data
  nNodes <- length(state[[state[[dv]]$nodeSet[1]]]$ids)
  set.seed(seed)
  E <- sample(nrow(el), n, TRUE)
  J <- sample(nNodes, n, TRUE)
  loopMoves <- seq_len(max(1, n %/% 10))
  J[loopMoves] <- el[E[loopMoves], 1]
  fromLoops <- which(el[, 1] == el[, 2])
  if (length(fromLoops) > 0) {
    E[loopMoves + length(loopMoves)] <- fromLoops[seq_along(loopMoves)]
  }
  U <- sample(c(-1, 1), n, TRUE)

  for (k in seq_along(effects$name)) {
    f <- effects$effectFormulas[[k]]
    for (m in 1:n) {
      expected <- c(
        change = f(dep.var = dv, state = state, cache = state$cache,
                   i = el[E[m], 1], j = J[m], edge = E[m], update = U[m]),
        target = f(dep.var = dv, state = state, cache = state$cache,
                   i = el[E[m], 1], j = J[m], edge = E[m], update = U[m],
                   getTargetContribution = TRUE)
      )
      expect_identical(
        cppEffectValues(effects, k, state, state$cache, el[E[m], 1], J[m], E[m], U[m]),
        expected,
        label = paste(effects$name[k], "edge", E[m], "to", J[m], "update", U[m])
      )
    }
  }
  # whole target statistics, which some C++ effects calculate differently
  old <- options(MoNAn2.useCpp = FALSE)
  rStatistics <- getNetworkStatistics(dv, state, state$cache, effects)
  options(old)
  expect_identical(cppNetworkStatistics(dv, state, state$cache, effects), rStatistics)
}
