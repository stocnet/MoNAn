# parallel estimation: the workers load the installed MoNAn, so these tests
# only run when the installed package is tested (devtools::check, R CMD check
# with NOT_CRAN=true)
skipParallel <- function() {
  skip_on_cran()
  skip_if_not(testthat::is_checking(), "workers need the installed package")
}

dv <- myState$dep.var
smallAlg <- createAlgorithm(myState, myEffects, burnInN1 = 500, thinningN1 = 100,
                            burnInN2 = 500, thinningN2 = 100, nsubN2 = 1,
                            initialIterationsN2 = 20, burnInN3 = 500,
                            thinningN3 = 100, iterationsN3 = 101)

test_that("chains get distinct random number streams from the session's random numbers", {
  set.seed(3)
  s1 <- chainSeeds(4)
  kind <- RNGkind()
  set.seed(3)
  s2 <- chainSeeds(4)
  expect_identical(s1, s2)
  expect_identical(RNGkind(), kind)
  expect_length(unique(s1), 4)
  expect_identical(s1[[2]], parallel::nextRNGStream(s1[[1]]))
  expect_false(identical(chainSeeds(4)[[1]], s1[[1]]))
})

test_that("parallel estimation stores statistics in the order of the simulations", {
  skipParallel()
  set.seed(5)
  res <- estimateMobilityNetwork(myState, myEffects, smallAlg, parallel = TRUE,
                                 cpus = 2, returnDeps = TRUE)
  expect_length(res$deps, 101)
  expect_equal(dim(res$simulatedStatistics), c(101, length(myEffects$name)))
  expect_equal(unname(res$observedStatistics),
               unname(getNetworkStatistics(dv, res$state, res$cache, myEffects)))
  for (k in seq_along(res$deps)) {
    expect_equal(unname(res$simulatedStatistics[k, ]),
                 unname(getNetStatsFromDeps(dv, res$deps[[k]], res, myEffects)))
  }
  # the two chains (51 and 50 simulations) are not copies of each other
  expect_false(identical(res$deps[1:50], res$deps[52:101]))
  expect_no_error(scoreTest(res, createEffects(myState) |> addEffect(transitivity_AC, alpha = 8)))
  expect_no_error(extractTraces(res, myEffects))
})

test_that("parallel estimation is reproducible with set.seed and differs otherwise", {
  skipParallel()
  set.seed(5)
  a <- estimateMobilityNetwork(myState, myEffects, smallAlg, parallel = TRUE,
                               cpus = 2, returnDeps = TRUE)
  set.seed(5)
  b <- estimateMobilityNetwork(myState, myEffects, smallAlg, parallel = TRUE,
                               cpus = 2, returnDeps = TRUE)
  c <- estimateMobilityNetwork(myState, myEffects, smallAlg, parallel = TRUE,
                               cpus = 2, returnDeps = TRUE)
  expect_identical(a$estimates, b$estimates)
  expect_identical(a$deps, b$deps)
  expect_false(identical(a$estimates, c$estimates))
  expect_false(identical(a$deps, c$deps))
})

test_that("parallel Geyer-Thompson estimation is reproducible", {
  skipParallel()
  alg <- createAlgorithm(myState, myEffects, phase2Method = "GT", burnInN2 = 500,
                         thinningN2 = 100, iterationsGT = 100, burnInN3 = 500,
                         thinningN3 = 100, iterationsN3 = 100)
  set.seed(1)
  a <- estimateMobilityNetwork(myState, myEffects, alg, parallel = TRUE, cpus = 2)
  set.seed(1)
  b <- estimateMobilityNetwork(myState, myEffects, alg, parallel = TRUE, cpus = 2)
  expect_identical(a$estimates, b$estimates)
  expect_true(all(is.finite(a$standardErrors)))
})

test_that("parallel and serial simulations agree within Monte Carlo error", {
  skipParallel()
  alg <- createAlgorithm(myState, myEffects)
  cl <- startCluster(2)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  phase3 <- function(cl) {
    runPhase3(dv, myState, myState$cache, myEffects, myResDN$estimates,
              myResDN$observedStatistics, iterationsN3 = 400,
              burnInN3 = alg$burnInN3, thinningN3 = alg$thinningN3, cl = cl,
              allowLoops = alg$allowLoops, returnDeps = FALSE,
              multinomialProposal = alg$multinomialProposal, fish = FALSE)$statistics
  }
  set.seed(7)
  par <- phase3(cl)
  serial <- phase3(NULL)
  z <- (colMeans(par) - colMeans(serial)) /
    sqrt(apply(par, 2, var) / nrow(par) + apply(serial, 2, var) / nrow(serial))
  expect_true(all(abs(z) < 4.5))
})
