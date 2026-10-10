# C++ versions of the concentration AC/GW effects (src/effects.cpp)

concentrationEffects <- createEffects(cppTestState) |>
  addEffect(concentration_GW, alpha = 2) |>
  addEffect(concentration_GW, alpha = 3) |>
  addEffect(concentration_AC_dyad_covar, node.attribute = "sameRegion", alpha = 2) |>
  addEffect(concentration_AC_dyad_covar, node.attribute = "size_similarity", alpha = 4) |>
  addEffect(concentration_GW_dyad_covar, node.attribute = "sameRegion", alpha = 2) |>
  addEffect(concentration_GW_dyad_covar, node.attribute = "size_similarity", alpha = 3) |>
  addEffect(concentration_AC_resource_covar_bin, edge.attribute = "sex", alpha = 2) |>
  addEffect(concentration_GW_resource_covar_bin, edge.attribute = "sex", alpha = 2)

test_that("C++ concentration effects give the same values as the R effects", {
  expect_cpp_matches_r(concentrationEffects, state = cppTestState)
})

test_that("C++ and R chains are identical with the concentration effects", {
  parameters <- c(0.3, 0.2, 0.2, 0.1, 0.1, 0.1, 0.2, -0.1)
  compare_chains(concentrationEffects, parameters, n = 5000, state = cppTestState)
  compare_chains(concentrationEffects, parameters, n = 500, state = cppTestState, multinomial = TRUE)
})
