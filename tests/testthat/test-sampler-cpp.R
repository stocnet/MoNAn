# The C++ sampler (src/sampler.cpp) must give exactly the same chain as the R
# sampler for the same seed: same edgelist, cache, statistics and random
# number state afterwards.

compare_chains <- function(effects, parameters, n = 5000, allowLoops = TRUE,
                           multinomial = FALSE, state = myState, seed = 42) {
  dv <- state$dep.var
  set.seed(seed)
  r <- simulateNSteps(dv, state, state$cache, effects, parameters,
                      allowLoops = allowLoops, multinomialProposal = multinomial, n = n)
  rNext <- runif(1)
  set.seed(seed)
  cpp <- simulateChain(dv, state, state$cache, effects, parameters,
                       burnIn = n, iterations = 0, thinning = 0,
                       allowLoops = allowLoops, multinomialProposal = multinomial)
  cppNext <- runif(1)

  expect_identical(cpp$state[[dv]]$data, toListState(r$state, dv)[[dv]]$data)
  expect_equal(cpp$cache, toListCache(r$cache), tolerance = 0)
  expect_identical(cppNext, rNext)
}

test_that("C++ and R chains are identical with the loops effect", {
  effects <- createEffects(myState) |> addEffect(loops)
  for (theta in c(0, 1.5, -1)) {
    compare_chains(effects, theta, allowLoops = TRUE)
    compare_chains(effects, theta, allowLoops = FALSE)
  }
})

test_that("C++ statistics equal the R statistics of each draw", {
  effects <- createEffects(myState) |> addEffect(loops)
  dv <- myState$dep.var
  set.seed(5)
  r <- simulateStatisticVectors(dv, myState, myState$cache, effects, 1,
                                burnIn = 200, iterations = 20, thinning = 200,
                                allowLoops = TRUE)
  set.seed(5)
  cpp <- simulateChain(dv, myState, myState$cache, effects, 1,
                       burnIn = 200, iterations = 20, thinning = 200, allowLoops = TRUE)
  expect_identical(unname(cpp$statistics), unname(r))
})

test_that("models with effects without a C++ version are recognised", {
  effects <- createEffects(myState) |> addEffect(loops) |> addEffect(transitivity_AC)
  expect_null(cppEffectSpec(effects, myState))
  expect_false(useCppSampler(effects, myState))
})
