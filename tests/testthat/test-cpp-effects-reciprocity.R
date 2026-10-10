# C++ versions of the reciprocity AC/GW effects (src/effects.cpp)

reciprocityEffects <- createEffects(cppTestState) |>
  addEffect(reciprocity_GW, alpha = 2) |>
  addEffect(reciprocity_GW, alpha = 3) |>
  addEffect(reciprocity_AC_dyad_covar, node.attribute = "size_similarity", alpha = 2) |>
  addEffect(reciprocity_AC_dyad_covar, node.attribute = "sameRegion", alpha = 4) |>
  addEffect(reciprocity_AC_dyad_covar_bin, node.attribute = "sameRegion", alpha = 2) |>
  addEffect(reciprocity_GW_dyad_covar, node.attribute = "size_similarity", alpha = 2) |>
  addEffect(reciprocity_GW_dyad_covar_bin, node.attribute = "sameRegion", alpha = 3)

test_that("C++ reciprocity effects give the same values as the R effects", {
  expect_cpp_matches_r(reciprocityEffects, state = cppTestState)
})

test_that("C++ and R chains are identical with the reciprocity effects", {
  parameters <- c(0.3, 0.1, 0.2, 0.1, 0.2, 0.1, -0.1)
  compare_chains(reciprocityEffects, parameters, n = 5000, state = cppTestState)
  compare_chains(reciprocityEffects, parameters, n = 500, state = cppTestState, multinomial = TRUE)
})
