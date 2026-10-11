# scoreTest: score test conditional on the statistics of the model's effects

candidates <- createEffects(myState) |>
  addEffect(transitivity_AC, alpha = 8) |>
  addEffect(transitivity_min)

test_that("estimation with returnDeps stores the statistics of the model's effects", {
  dv <- myResDN$state$dep.var
  expect_equal(dim(myResDN$simulatedStatistics), c(length(myResDN$deps), length(myEffects$name)))
  expect_equal(unname(myResDN$observedStatistics),
               unname(getNetworkStatistics(dv, myResDN$state, myResDN$cache, myEffects)))
  # the stored statistics belong to the stored networks, in the same order
  for (k in c(1, 500, 1000)) {
    expect_equal(unname(myResDN$simulatedStatistics[k, ]),
                 unname(getNetStatsFromDeps(dv, myResDN$deps[[k]], myResDN, myEffects)))
  }
})

test_that("estimation without returnDeps does not store statistics", {
  alg <- createAlgorithm(myState, myEffects, burnInN1 = 200, thinningN1 = 50,
                         burnInN2 = 200, thinningN2 = 50, nsubN2 = 1,
                         initialIterationsN2 = 10, burnInN3 = 200, thinningN3 = 50,
                         iterationsN3 = 20)
  set.seed(1)
  res <- estimateMobilityNetwork(myState, myEffects, alg)
  expect_null(res$simulatedStatistics)
  expect_null(res$observedStatistics)
})

test_that("scoreTest gives the conditional score test", {
  st <- scoreTest(myResDN, candidates)
  s1 <- myResDN$simulatedStatistics
  d1 <- myResDN$observedStatistics - colMeans(s1)
  for (k in seq_along(candidates$name)) {
    x <- st$simulatedStatistics[, k]
    c21 <- as.numeric(cov(x, s1))
    w <- as.numeric(c21 %*% solve(cov(s1)))
    v <- var(x) - sum(w * c21)
    t <- ((st$observedStatistics[k] - mean(x)) - sum(w * d1)) / sqrt(v)
    expect_equal(unname(st$tStatistics[k]), unname(t))
    expect_equal(unname(st$pValues[k]), 2 * pnorm(-abs(unname(t))))
    expect_equal(unname(st$R2[k]), 1 - v / var(x))
  }
  expect_null(st$pValuesNonParametric)
  expect_output(print(st), "R2")
})

test_that("results without stored statistics need modelEffects", {
  old <- myResDN
  old$simulatedStatistics <- NULL
  old$observedStatistics <- NULL
  expect_error(scoreTest(old, candidates), "modelEffects")
  expect_error(scoreTest(old, candidates, modelEffects = createEffects(myState) |> addEffect(loops)),
               "effects of the estimated model")
  expect_equal(scoreTest(old, candidates, modelEffects = myEffects)$pValues,
               scoreTest(myResDN, candidates)$pValues)
})

test_that("the non-parametric value of earlier versions is available as an option", {
  st <- scoreTest(myResDN, candidates, nonParametric = TRUE)
  expected <- vapply(seq_along(candidates$name), function(k) {
    mean(st$simulatedStatistics[, k] <= st$observedStatistics[k])
  }, numeric(1))
  expect_equal(unname(st$pValuesNonParametric), expected)
  expect_output(print(st), "nonParametric")
})

test_that("scoreTest requires stored simulations", {
  noDeps <- myResDN
  noDeps$deps <- NULL
  expect_error(scoreTest(noDeps, candidates), "returnDeps = TRUE")
})

test_that("extractTraces uses the stored statistics and gives the same traces", {
  stored <- extractTraces(myResDN, myEffects)
  old <- myResDN
  old$simulatedStatistics <- NULL
  old$observedStatistics <- NULL
  recomputed <- extractTraces(old, myEffects)
  expect_equal(unname(stored$observedStats), unname(recomputed$observedStats))
  expect_equal(unname(stored$simulatedStats), unname(recomputed$simulatedStats))
  expect_identical(stored$effectNames, recomputed$effectNames)
  # other effects than those of the model are calculated from the simulations
  other <- extractTraces(myResDN, createEffects(myState) |> addEffect(loops))
  expect_length(other$observedStats, 1)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)
  expect_no_error(plot(stored))
})
