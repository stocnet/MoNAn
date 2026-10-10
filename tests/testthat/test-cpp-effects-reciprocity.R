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

otherReciprocityEffects <- createEffects(cppTestState) |>
  addEffect(reciprocity_basic) |>
  addEffect(reciprocity_min) |>
  addEffect(reciprocity_min_dyad_covar, node.attribute = "sameRegion") |>
  addEffect(reciprocity_min_dyad_covar, node.attribute = "size_similarity") |>
  addEffect(reciprocity_min_resource_covar, edge.attribute = "sex") |>
  addEffect(reciprocity_min_resource_covar, edge.attribute = "second_or")

test_that("C++ basic and min reciprocity effects give the same values as R", {
  expect_cpp_matches_r(otherReciprocityEffects, state = cppTestState)
})

test_that("C++ and R chains are identical with the basic and min reciprocity effects", {
  parameters <- c(0.05, 0.2, 0.1, 0.1, 0.1, 0.005)
  compare_chains(otherReciprocityEffects, parameters, n = 5000, state = cppTestState)
  compare_chains(otherReciprocityEffects, parameters, n = 500, state = cppTestState,
                 multinomial = TRUE)
})
