# C++ versions of the loops effects (src/effects.cpp)

loopsEffects <- createEffects(myState) |>
  addEffect(loops_AC, alpha = 2) |>
  addEffect(loops_AC, alpha = 1.5) |>
  addEffect(loops_GW, alpha = 2) |>
  addEffect(loops_GW, alpha = 3) |>
  addEffect(loops_additional_origin, edge.attribute = "second_or") |>
  addEffect(loops_x_loops_additional_origin, edge.attribute = "second_or") |>
  addEffect(loops_node_covar, node.attribute = "size") |>
  addEffect(loops_resource_covar_node_covar, edge.attribute = "sex", node.attribute = "size")

test_that("C++ loops effects give the same values as the R effects", {
  expect_cpp_matches_r(loopsEffects)
})

test_that("C++ and R chains are identical with the loops effects", {
  parameters <- c(0.5, 0.3, 0.4, -0.2, 0.6, 0.8, 0.05, -0.1)
  compare_chains(loopsEffects, parameters, n = 5000)
  compare_chains(loopsEffects, parameters, n = 5000, allowLoops = FALSE)
  compare_chains(loopsEffects, parameters, n = 500, multinomial = TRUE)
})
