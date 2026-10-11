# C++ versions of the exogenous, in-ties, present-relations and
# covariate-based effects (src/effects.cpp)

otherEffects <- createEffects(cppTestState) |>
  addEffect(same_covariate, node.attribute = "region") |>
  addEffect(sim_covariate, node.attribute = "size") |>
  addEffect(dyadic_covariate_resource_attribute, node.attribute = "sameRegion",
            edge.attribute = "sex") |>
  addEffect(dyadic_covariate_resource_attribute, node.attribute = "size_similarity",
            edge.attribute = "second_or") |>
  addEffect(in_ties_loops) |>
  addEffect(present_relations)

covariateEffects <- createEffects(cppTestState) |>
  addEffect(staying_by_prop_bin_inflow, edge.attribute = "sex") |>
  addEffect(crowding_out_prop_covar_bin, edge.attribute = "sex") |>
  addEffect(joining_similar_avoiding_dissimilar_covar_bin, edge.attribute = "sex") |>
  addEffect(avoiding_dissimilar_covar_bin, edge.attribute = "sex") |>
  addEffect(joining_similar_avoiding_dissimilar_covar_cont, edge.attribute = "second_or") |>
  addEffect(avoiding_dissimilar_covar_cont, edge.attribute = "second_or")

test_that("C++ exogenous, in-ties and present-relations effects give the same values as R", {
  expect_cpp_matches_r(otherEffects, state = cppTestState)
})

test_that("C++ covariate-based effects give the same values as R", {
  expect_cpp_matches_r(covariateEffects, state = cppTestState)
})

test_that("C++ and R chains are identical with these effects", {
  compare_chains(otherEffects, c(0.3, 0.2, 0.1, 0.01, 0.001, 0.2),
                 n = 5000, state = cppTestState)
  compare_chains(otherEffects, c(0.3, 0.2, 0.1, 0.01, 0.001, 0.2),
                 n = 500, state = cppTestState, multinomial = TRUE)
  compare_chains(covariateEffects, c(0.3, -0.2, 0.5, 0.3, 0.2, -0.05),
                 n = 5000, state = cppTestState)
  compare_chains(covariateEffects, c(0.3, -0.2, 0.5, 0.3, 0.2, -0.05),
                 n = 500, state = cppTestState, multinomial = TRUE)
})

test_that("all exported effects have a C++ version", {
  ns <- asNamespace("MoNAn")
  isEffect <- vapply(getNamespaceExports(ns), function(name) {
    f <- get(name, envir = ns)
    is.function(f) && all(c("state", "cache", "getTargetContribution") %in% names(formals(f)))
  }, logical(1))
  expect_setequal(names(isEffect)[isEffect], cppEffects)
})
