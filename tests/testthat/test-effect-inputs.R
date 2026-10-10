# Input checks of effects run once (when an effect is added and when
# estimation or simulation starts), not in every call of the effect function.

test_that("effects check their covariates when they are added", {
  checkState <- myState
  asym <- matrix(runif(17 * 17), 17, 17)
  checkState$asym <- createNetwork(asym, nodeSet = c("organisations", "organisations"))
  checkState$nonbinary <- createNetwork(asym + t(asym), nodeSet = c("organisations", "organisations"))
  checkState$sex_counts <- createNodeVariable(rpois(742, 2), nodeSet = "people",
                                              addSame = FALSE, addSim = FALSE)
  # addEffect() finds the state by its name, which does not see objects that
  # are local to this test
  assign("checkState", checkState, envir = globalenv())
  on.exit(rm("checkState", envir = globalenv()), add = TRUE)

  base <- createEffects(checkState)
  expect_error(addEffect(base, reciprocity_min_dyad_covar, node.attribute = "asym"),
               "must be symmetric")
  expect_error(addEffect(base, reciprocity_AC_dyad_covar, node.attribute = "asym"),
               "must be symmetric")
  expect_error(addEffect(base, reciprocity_GW_dyad_covar, node.attribute = "asym"),
               "must be symmetric")
  expect_error(addEffect(base, reciprocity_AC_dyad_covar_bin, node.attribute = "nonbinary"),
               "must be 0 or 1")
  expect_error(addEffect(base, reciprocity_GW_dyad_covar_bin, node.attribute = "nonbinary"),
               "must be 0 or 1")
  expect_error(addEffect(base, crowding_out_prop_covar_bin, edge.attribute = "sex_counts"),
               "binary covariates")
  expect_error(addEffect(base, staying_by_prop_bin_inflow, edge.attribute = "sex_counts"),
               "binary covariates")

  expect_no_error(addEffect(base, reciprocity_AC_dyad_covar_bin, node.attribute = "sameRegion"))
  expect_no_error(addEffect(base, crowding_out_prop_covar_bin, edge.attribute = "sex"))
})

test_that("estimation and simulation check effect inputs against their state", {
  # effects created for one state, used with a state whose covariate is invalid
  effects <- createEffects(myState) |>
    addEffect(reciprocity_min_dyad_covar, node.attribute = "sameRegion")
  badState <- myState
  badState$sameRegion$data[1, 2] <- 0
  badState$sameRegion$data[2, 1] <- 1

  alg <- createAlgorithm(badState, effects)
  expect_error(estimateMobilityNetwork(badState, effects, alg), "must be symmetric")
  expect_error(simulateMobilityNetworks(badState, effects, parameters = 0,
                                        allowLoops = TRUE, burnin = 10, thinning = 10,
                                        nSimulations = 1),
               "must be symmetric")
})

test_that("effects work when the dependent variable is not called 'transfers'", {
  renamed <- lililiState
  names(renamed)[names(renamed) == "transfers"] <- "moves"
  renamed$dep.var <- "moves"
  cache <- createInternalCache(renamed, resourceCovariates = c("sex", "second_or"))
  originalCache <- createInternalCache(lililiState, resourceCovariates = c("sex", "second_or"))
  el <- renamed$moves$data

  for (effect in c("avoiding_dissimilar_covar_cont",
                   "joining_similar_avoiding_dissimilar_covar_cont")) {
    f <- get(effect)
    for (edge in c(1, 100, 500)) {
      expect_identical(
        f("moves", "second_or", renamed, cache, el[edge, 1], el[edge, 2], edge, -1),
        f("transfers", "second_or", lililiState, originalCache, el[edge, 1], el[edge, 2], edge, -1)
      )
    }
    expect_identical(
      f("moves", "second_or", renamed, cache, 1, 2, 1, 1, getTargetContribution = TRUE),
      f("transfers", "second_or", lililiState, originalCache, 1, 2, 1, 1, getTargetContribution = TRUE)
    )
  }
})
