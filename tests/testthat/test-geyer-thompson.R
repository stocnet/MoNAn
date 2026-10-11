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
