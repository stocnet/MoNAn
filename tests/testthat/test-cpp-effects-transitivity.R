# C++ versions of the transitivity and triad effects (src/effects.cpp)

transitivityEffects <- createEffects(myState) |>
  addEffect(transitivity_basic) |>
  addEffect(transitivity_min) |>
  addEffect(transitivity_netflow) |>
  addEffect(transitivity_AC, alpha = 1.1) |>
  addEffect(transitivity_AC, alpha = 2) |>
  addEffect(transitivity_GW, alpha = 1.5) |>
  addEffect(transitivity_GW, alpha = 1) |>
  addEffect(triad120C) |>
  addEffect(triad120D) |>
  addEffect(triad120U)

test_that("C++ transitivity and triad effects give the same values as the R effects", {
  expect_cpp_matches_r(transitivityEffects)
})

test_that("C++ and R chains are identical with the transitivity and triad effects", {
  parameters <- c(0.001, 0.01, 0.005, 0.1, 0.05, 0.05, 0.02, 0.0001, 0.0001, 0.0001)
  compare_chains(transitivityEffects, parameters, n = 3000)
  compare_chains(transitivityEffects, parameters, n = 300, multinomial = TRUE)
})
