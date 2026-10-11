# For every effect, the change statistic of moving a resource must equal the
# difference in target statistics before and after the move. A mismatch biases
# estimates without raising an error.
#
# Uses the internal example data lililiState / lililiCache (R/sysdata.rda):
#   node covariates:     size (continuous), region (binary)
#   dyadic covariates:   sameRegion (binary), size_similarity (continuous)
#   resource covariates: sex (binary), second_or (origin id, used as continuous)
#
# Set MONAN_FULL_EFFECT_TESTS=true to move every resource instead of a sample.

effect_specs <- list(
  # exogenous
  alter_covariate = list(attribute.index = "size"),
  alter_covariate_no_loops = list(.effect = "alter_covariate", attribute.index = "size",
                                  loop.contribution = FALSE),
  dyadic_covariate = list(attribute.index = "sameRegion"),
  dyadic_covariate_resource_attribute = list(attribute.index = "sameRegion",
                                             resource.attribute.index = "sex"),
  resource_covar_to_node_covar = list(resource.attribute.index = "sex", attribute.index = "region"),
  resource_covar_to_node_covar_loops = list(.effect = "resource_covar_to_node_covar",
                                            resource.attribute.index = "sex",
                                            attribute.index = "region",
                                            loop.contribution = TRUE),
  same_covariate = list(attribute.index = "region"),
  sim_covariate = list(attribute.index = "size"),

  # loops
  loops = list(),
  loops_AC = list(alpha = 2),
  loops_AC_alpha4 = list(.effect = "loops_AC", alpha = 4),
  loops_GW = list(alpha = 2),
  loops_GW_alpha4 = list(.effect = "loops_GW", alpha = 4),
  loops_additional_origin = list(resource.attribute.index = "second_or"),
  loops_node_covar = list(attribute.index = "size"),
  loops_resource_covar = list(resource.attribute.index = "sex"),
  loops_resource_covar_node_covar = list(resource.attribute.index = "sex", attribute.index = "size"),
  loops_x_loops_additional_origin = list(resource.attribute.index = "second_or"),

  # concentration
  concentration_AC = list(alpha = 2),
  concentration_AC_alpha4 = list(.effect = "concentration_AC", alpha = 4),
  concentration_AC_dyad_covar = list(attribute.index = "sameRegion", alpha = 2),
  concentration_AC_resource_covar_bin = list(resource.attribute.index = "sex", alpha = 2),
  concentration_GW = list(alpha = 2),
  concentration_GW_alpha4 = list(.effect = "concentration_GW", alpha = 4),
  concentration_GW_dyad_covar = list(attribute.index = "sameRegion", alpha = 2),
  concentration_GW_resource_covar_bin = list(resource.attribute.index = "sex", alpha = 2),
  concentration_basic = list(),
  concentration_basic_cube = list(),
  concentration_max = list(),
  concentration_basic_squared = list(),
  concentration_norm = list(),
  concentration_norm_squared = list(),
  concentration_prop_orig_cov = list(attribute.index = "size"),
  concentration_rankGW = list(lambda = 2),
  concentration_rankGW_lambda = list(.effect = "concentration_rankGW", lambda = 1.5),
  concentration_prop = list(),
  concentration_prop_AC = list(alpha = 2),
  concentration_prop_AC_alpha4 = list(.effect = "concentration_prop_AC", alpha = 4),

  # in-weights
  in_ties_loops = list(),
  in_weights_AC = list(alpha = 2),
  in_weights_AC_alpha4 = list(.effect = "in_weights_AC", alpha = 4),
  in_weights_GW = list(alpha = 2),
  in_weights_GW_alpha4 = list(.effect = "in_weights_GW", alpha = 4),
  present_relations = list(),

  # reciprocity
  reciprocity_AC = list(alpha = 2),
  reciprocity_AC_alpha4 = list(.effect = "reciprocity_AC", alpha = 4),
  reciprocity_AC_dyad_covar = list(attribute.index = "size_similarity", alpha = 2),
  reciprocity_AC_dyad_covar_bin = list(attribute.index = "sameRegion", alpha = 2),
  reciprocity_GW = list(alpha = 2),
  reciprocity_GW_alpha4 = list(.effect = "reciprocity_GW", alpha = 4),
  reciprocity_GW_dyad_covar = list(attribute.index = "size_similarity", alpha = 2),
  reciprocity_GW_dyad_covar_bin = list(attribute.index = "sameRegion", alpha = 2),
  reciprocity_basic = list(),
  reciprocity_min = list(),
  reciprocity_min_dyad_covar = list(attribute.index = "sameRegion"),
  reciprocity_min_resource_covar = list(resource.attribute.index = "sex"),

  # transitivity and triads
  transitivity_AC = list(alpha = 1.1),
  transitivity_AC_alpha4 = list(.effect = "transitivity_AC", alpha = 4),
  transitivity_GW = list(alpha = 1.1),
  transitivity_GW_alpha4 = list(.effect = "transitivity_GW", alpha = 4),
  transitivity_basic = list(),
  transitivity_min = list(),
  transitivity_netflow = list(),
  triad120C = list(),
  triad120D = list(),
  triad120U = list(),

  # endogenous, covariate-based
  avoiding_dissimilar_covar_bin = list(resource.attribute.index = "sex"),
  avoiding_dissimilar_covar_cont = list(resource.attribute.index = "second_or"),
  crowding_out_prop_covar_bin = list(resource.attribute.index = "sex"),
  joining_similar_avoiding_dissimilar_covar_bin = list(resource.attribute.index = "sex"),
  joining_similar_avoiding_dissimilar_covar_cont = list(resource.attribute.index = "second_or"),
  staying_by_prop_bin_inflow = list(resource.attribute.index = "sex")
)

# effects with a known mismatch between change and target statistics,
# as c(effect = "reason"); these are skipped until fixed
known_bugs <- c()

# same_covariate needs the categorical 'same' matrix, which lililiState lacks
same_state <- lililiState
same_state$region <- createNodeVariable(lililiState$region$data,
                                        nodeSet = "organisations", addSame = TRUE)

moves <- sample_moves()

for (spec_name in names(effect_specs)) {
  local({
    spec <- effect_specs[[spec_name]]
    effect <- if (is.null(spec$.effect)) spec_name else spec$.effect
    args <- spec[names(spec) != ".effect"]
    state <- if (effect == "same_covariate") same_state else lililiState

    test_that(paste(spec_name, "change statistic matches target difference"), {
      if (effect %in% names(known_bugs)) {
        skip(paste0("known bug in ", effect, ": ", known_bugs[[effect]]))
      }
      expect_change_matches(effect, args, moves, state = state)
    })
  })
}

# regression test: moves into a loop by resources whose additional origin is
# their origin used to fail with NA in loops_x_loops_additional_origin
test_that("additional-origin loop effects handle moves into a loop", {
  el <- lililiState$transfers$data
  candidates <- which(el[, 1] != el[, 2] & lililiState$second_or$data == el[, 1])
  expect_gt(length(candidates), 0)
  loop_moves <- data.frame(edge = candidates[1:3], j_new = el[candidates[1:3], 1])
  for (effect in c("loops_additional_origin", "loops_x_loops_additional_origin")) {
    expect_change_matches(effect, list(resource.attribute.index = "second_or"), loop_moves)
  }
})

test_that("every exported effect has a change/target test", {
  ns <- asNamespace("MoNAn")
  is_effect <- vapply(getNamespaceExports(ns), function(name) {
    f <- get(name, envir = ns)
    is.function(f) && all(c("state", "cache", "getTargetContribution") %in% names(formals(f)))
  }, logical(1))
  tested <- vapply(names(effect_specs), function(n) {
    if (is.null(effect_specs[[n]]$.effect)) n else effect_specs[[n]]$.effect
  }, character(1))

  expect_setequal(names(is_effect)[is_effect], unique(tested))
})
