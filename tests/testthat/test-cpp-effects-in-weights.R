# C++ versions of the in-weights AC/GW effects (src/effects.cpp)

inWeightsEffects <- createEffects(myState) |>
  addEffect(in_weights_AC, alpha = 2) |>
  addEffect(in_weights_AC, alpha = 5) |>
  addEffect(in_weights_GW, alpha = 2) |>
  addEffect(in_weights_GW, alpha = 1.5)

test_that("C++ in-weights effects give the same values as the R effects", {
  expect_cpp_matches_r(inWeightsEffects)
})

test_that("C++ and R chains are identical with the in-weights effects", {
  parameters <- c(0.2, 0.1, 0.05, 0.02)
  compare_chains(inWeightsEffects, parameters, n = 5000)
  compare_chains(inWeightsEffects, parameters, n = 500, multinomial = TRUE)
})
