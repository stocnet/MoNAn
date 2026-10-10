# The C++ sampler (src/sampler.cpp) must give exactly the same chain as the R
# sampler for the same seed: same edgelist, cache, statistics and random
# number state afterwards.

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
  old <- options(MoNAn2.useCpp = FALSE)
  on.exit(options(old), add = TRUE)
  set.seed(5)
  r <- simulateStatisticVectors(dv, myState, myState$cache, effects, 1,
                                burnIn = 200, iterations = 20, thinning = 200,
                                allowLoops = TRUE)
  options(MoNAn2.useCpp = TRUE)
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

test_that("C++ effects give the same values as the R effects", {
  effects <- createEffects(myState) |>
    addEffect(loops) |>
    addEffect(concentration_AC, alpha = 4) |>
    addEffect(reciprocity_AC, alpha = 2) |>
    addEffect(dyadic_covariate, node.attribute = "sameRegion") |>
    addEffect(alter_covariate, node.attribute = "size") |>
    addEffect(alter_covariate, node.attribute = "size", loop.contribution = FALSE) |>
    addEffect(resource_covar_to_node_covar, node.attribute = "region", edge.attribute = "sex") |>
    addEffect(resource_covar_to_node_covar, node.attribute = "size",
              edge.attribute = "second_or", loop.contribution = TRUE) |>
    addEffect(loops_resource_covar, edge.attribute = "second_or")
  dv <- myState$dep.var
  el <- myState[[dv]]$data
  set.seed(3)
  n <- 100
  E <- sample(nrow(el), n, TRUE)
  J <- sample(17, n, TRUE)
  J[1:10] <- el[E[1:10], 1]  # moves into loops
  U <- sample(c(-1, 1), n, TRUE)

  for (k in seq_along(effects$name)) {
    f <- effects$effectFormulas[[k]]
    for (m in 1:n) {
      expected <- c(
        change = f(dep.var = dv, state = myState, cache = myState$cache,
                   i = el[E[m], 1], j = J[m], edge = E[m], update = U[m]),
        target = f(dep.var = dv, state = myState, cache = myState$cache,
                   i = el[E[m], 1], j = J[m], edge = E[m], update = U[m],
                   getTargetContribution = TRUE)
      )
      expect_identical(
        cppEffectValues(effects, k, myState, myState$cache, el[E[m], 1], J[m], E[m], U[m]),
        expected,
        label = paste(effects$name[k], "edge", E[m], "to", J[m])
      )
    }
  }
})

test_that("C++ and R chains are identical for the example model", {
  compare_chains(myEffects, myResDN$estimates, n = 5000)
  compare_chains(myEffects, myResDN$estimates, n = 5000, allowLoops = FALSE)
  compare_chains(myEffects, rep(0, 7), n = 5000)
})

test_that("C++ and R chains are identical with the multinomial proposal", {
  compare_chains(myEffects, myResDN$estimates, n = 1000, multinomial = TRUE)
  compare_chains(myEffects, myResDN$estimates, n = 1000, multinomial = TRUE, allowLoops = FALSE)
  compare_chains(myEffects, rep(0, 7), n = 1000, multinomial = TRUE)
})

test_that("estimation gives the same results with the C++ and the R sampler", {
  skip_on_cran()
  alg <- createAlgorithm(myState, myEffects, burnInN1 = 300, thinningN1 = 100,
                         burnInN2 = 300, thinningN2 = 100, nsubN2 = 1,
                         initialIterationsN2 = 15, burnInN3 = 300, thinningN3 = 100,
                         iterationsN3 = 30)
  run <- function(useCpp) {
    old <- options(MoNAn2.useCpp = useCpp)
    on.exit(options(old))
    set.seed(8)
    estimateMobilityNetwork(myState, myEffects, alg)
  }
  expect_true(useCppSampler(myEffects, myState))
  resR <- run(FALSE)
  resCpp <- run(TRUE)
  expect_equal(resCpp$estimates, resR$estimates, tolerance = 1e-12)
  expect_equal(resCpp$covarianceMatrix, resR$covarianceMatrix, tolerance = 1e-10)
  expect_equal(resCpp$convergenceStatistics, resR$convergenceStatistics, tolerance = 1e-10)
})
