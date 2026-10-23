# Geyer-Thompson estimation (phase2Method = "GT") and its MPLE starting values

test_that("Geyer-Thompson step solves the moment equation in the weighted sample", {
  set.seed(1)
  # normal sample: exponential tilting by delta shifts the mean by Sigma %*% delta
  Sigma <- matrix(c(4, 1, 1, 2), 2)
  S <- matrix(rnorm(10000), ncol = 2) %*% chol(Sigma) +
    matrix(c(10, 20), 5000, 2, byrow = TRUE)
  obs <- c(10.6, 19.7)

  step <- geyerThompsonStep(S, obs, minESS = 0.1)
  w <- exp(S %*% step$delta)
  w <- as.vector(w / sum(w))

  expect_equal(step$method, "Geyer-Thompson")
  expect_false(step$restricted)
  expect_equal(colSums(w * S), obs, tolerance = 1e-6)
  expect_equal(step$delta, as.vector(solve(Sigma, obs - c(10, 20))), tolerance = 0.1)
})

test_that("Geyer-Thompson step is restricted when observed statistics are far away", {
  set.seed(2)
  S <- matrix(rnorm(2000), ncol = 2)

  far <- geyerThompsonStep(S, c(30, 0))
  expect_equal(far$method, "normal")
  expect_true(far$restricted)
  expect_true(all(is.finite(far$delta)))
  # the step moves towards the observed statistics, but not all the way
  expect_gt(far$delta[1], 0)
  expect_lt(far$delta[1], 30)

  near <- geyerThompsonStep(S, c(0.5, 0), minESS = 0.5)
  expect_equal(near$method, "Geyer-Thompson")
  expect_gte(near$ess, 0.5 * nrow(S) - 1e-8)
})

test_that("MPLE is finite and close to the reference estimates", {
  mple <- getMPLE(myState$dep.var, myState, myState$cache, myEffects, allowLoops = TRUE)

  expect_named(mple, myEffects$name)
  expect_true(all(is.finite(mple)))
  # the MPLE is a starting value; it should be within a few standard errors
  expect_true(all(abs(mple - myResDN$estimates) < 3 * myResDN$standardErrors))
})

test_that("MPLE does not change the state", {
  state <- myState
  frozen <- unserialize(serialize(state, NULL))
  getMPLE(state$dep.var, state, state$cache, myEffects, allowLoops = TRUE)
  expect_identical(state, frozen)
})

test_that("estimation with phase2Method = 'GT' runs and returns a result", {
  skip_on_cran()
  alg <- createAlgorithm(myState, myEffects, phase2Method = "GT",
                         burnInN2 = 500, thinningN2 = 200,
                         iterationsGT = 30, finalIterationsGT = 60, maxIterationsGT = 2,
                         burnInN3 = 500, thinningN3 = 200, iterationsN3 = 40)
  set.seed(3)
  invisible(capture.output(
    res <- estimateMobilityNetwork(myState, myEffects, alg)
  ))

  expect_s3_class(res, "result.monan")
  expect_named(res$estimates, myEffects$name)
  expect_true(all(is.finite(res$estimates)))
  expect_true(all(is.finite(res$standardErrors)))
})

test_that("algorithm objects without phase2Method use Robbins-Monro", {
  alg <- createAlgorithm(myState, myEffects)
  expect_equal(alg$phase2Method, "RM")
  expect_null(myAlg$phase2Method)
  expect_error(createAlgorithm(myState, myEffects, phase2Method = "XY"))
})

test_that("the C++ pseudo-likelihood equals the R pseudo-likelihood", {
  resourceEffects <- createEffects(cppTestState) |>
    addEffect(loops) |>
    addEffect(concentration_GW, alpha = 2) |>
    addEffect(reciprocity_min_resource_covar, edge.attribute = "second_or") |>
    addEffect(loops_resource_covar, edge.attribute = "sex") |>
    addEffect(transitivity_AC)
  for (effects in list(myEffects, resourceEffects)) {
    dv <- cppTestState$dep.var
    fR <- mpleObjectiveR(dv, cppTestState, cppTestState$cache, effects, TRUE)
    fC <- mpleObjectiveCpp(dv, cppTestState, cppTestState$cache, effects, TRUE)
    set.seed(1)
    for (theta in list(rep(0, length(effects$name)), runif(length(effects$name), -0.5, 0.5))) {
      r <- fR(theta)
      cc <- fC(theta)
      expect_equal(cc$logPL, r$logPL, tolerance = 1e-10)
      expect_equal(cc$gradient, r$gradient, tolerance = 1e-10)
      expect_equal(cc$hessian, r$hessian, tolerance = 1e-10)
      expect_equal(fC(theta, derivatives = FALSE)$logPL, r$logPL, tolerance = 1e-10)
    }
  }
})

test_that("MPLE is the same with the C++ and the R pseudo-likelihood", {
  dv <- myState$dep.var
  old <- options(MoNAn.useCpp = FALSE)
  mpleR <- getMPLE(dv, myState, myState$cache, myEffects, allowLoops = TRUE)
  options(old)
  mpleCpp <- getMPLE(dv, myState, myState$cache, myEffects, allowLoops = TRUE)
  expect_equal(mpleCpp, mpleR, tolerance = 1e-8)

  mpleNoLoops <- getMPLE(dv, myState, myState$cache,
                         createEffects(myState) |> addEffect(concentration_AC, alpha = 2),
                         allowLoops = FALSE)
  # the observed data contain loops, which are impossible without loops
  expect_null(mpleNoLoops)
})

test_that("the MPLE stores each distinct choice set once, weighted by its size", {
  dv <- myState$dep.var
  inputs <- cppChainInputs(dv, myState, myState$cache)
  sets <- mpleChoiceSetsCpp(inputs$edgelist, inputs$cache, inputs$resourceCovariates,
                            inputs$resourceSet, inputs$nodeSet,
                            cppEffectSpec(myEffects, inputs$state), "sex", TRUE, 1e9)
  el <- myState[[dv]]$data
  # myEffects use the resource covariate sex only
  expect_length(sets$weight, nrow(unique(cbind(el, myState$sex$data))))
  expect_equal(sum(sets$weight), nrow(el))
  expect_equal(diff(sets$start), rep(length(myState$organisations$ids), length(sets$weight)))
  # too little memory: nothing is computed
  small <- mpleChoiceSetsCpp(inputs$edgelist, inputs$cache, inputs$resourceCovariates,
                             inputs$resourceSet, inputs$nodeSet,
                             cppEffectSpec(myEffects, inputs$state), "sex", TRUE, 1000)
  expect_true(small$tooLarge)
  expect_null(small$X)
})

test_that("grouped choice sets give the same pseudo-likelihood as recomputing them", {
  # effects with resource covariates, including a continuous one, and effects
  # that look at the individual resources of a cell
  groupEffects <- createEffects(cppTestState) |>
    addEffect(loops_additional_origin, edge.attribute = "second_or") |>
    addEffect(joining_similar_avoiding_dissimilar_covar_cont, edge.attribute = "second_or") |>
    addEffect(concentration_AC_resource_covar_bin, edge.attribute = "sex", alpha = 2) |>
    addEffect(transitivity_AC, alpha = 4)
  dv <- cppTestState$dep.var
  for (effects in list(myEffects, groupEffects)) {
    for (allowLoops in c(TRUE, FALSE)) {
      state <- cppTestState
      if (!allowLoops) {
        # without loops, the observed destinations must not be loops
        el <- state[[dv]]$data
        el[el[, 1] == el[, 2], 2] <- el[el[, 1] == el[, 2], 2] %% 17 + 1
        state[[dv]] <- createEdgelist(el, nodeSet = state[[dv]]$nodeSet)
        state$cache <- createInternalCache(state, resourceCovariates = names(state$cache[[dv]]$resourceNetworks))
      }
      grouped <- mpleObjectiveCpp(dv, state, state$cache, effects, allowLoops)
      recomputed <- mpleObjectiveCpp(dv, state, state$cache, effects, allowLoops, memoryLimit = 0)
      set.seed(2)
      for (theta in list(rep(0, length(effects$name)), runif(length(effects$name), -0.5, 0.5))) {
        g <- grouped(theta)
        r <- recomputed(theta)
        expect_equal(g$logPL, r$logPL, tolerance = 1e-10)
        expect_equal(g$gradient, r$gradient, tolerance = 1e-10)
        expect_equal(g$hessian, r$hessian, tolerance = 1e-10)
        expect_equal(grouped(theta, derivatives = FALSE)$logPL, r$logPL, tolerance = 1e-10)
      }
    }
  }
})

test_that("the MPLE is the same with grouped and recomputed choice sets", {
  dv <- myState$dep.var
  grouped <- getMPLE(dv, myState, myState$cache, myEffects, allowLoops = TRUE)
  old <- options(MoNAn.mpleMemoryLimit = 0)
  on.exit(options(old), add = TRUE)
  recomputed <- getMPLE(dv, myState, myState$cache, myEffects, allowLoops = TRUE)
  expect_equal(grouped, recomputed, tolerance = 1e-10)
  expect_null(getMPLE(dv, myState, myState$cache,
                      createEffects(myState) |> addEffect(concentration_AC, alpha = 2),
                      allowLoops = FALSE))
})

test_that("MPLE uses the R pseudo-likelihood for effects written only in R", {
  assign("myLoopsMPLE", function(dep.var = 1, state, cache, i, j, edge, update,
                                 getTargetContribution = FALSE) {
    if (getTargetContribution) return(cache[[dep.var]]$valuedNetwork[i, j] * (i == j))
    if (i == j) return(update)
    0
  }, envir = globalenv())
  on.exit(rm("myLoopsMPLE", envir = globalenv()), add = TRUE)
  effects <- createEffects(myState) |> addEffect(myLoopsMPLE)
  expect_false(useCppSampler(effects, myState))
  mple <- getMPLE(myState$dep.var, myState, myState$cache, effects, allowLoops = TRUE)
  # same model as the loops effect alone
  expect_equal(unname(mple),
               unname(getMPLE(myState$dep.var, myState, myState$cache,
                              createEffects(myState) |> addEffect(loops), allowLoops = TRUE)),
               tolerance = 1e-8)
})
