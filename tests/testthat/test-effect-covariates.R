# The covariates of effects must be covariates of the state of the kind the
# effect reads. Effects index their data directly, so a nodeset, the
# dependent variable or a covariate of the wrong kind gave wrong statistics
# without a message, or crashed R in C++ (by reading outside the data).

# the example state with a 'same' matrix for size (for same_covariate)
covariateState <- myState
covariateState$sizeSame <- createNodeVariable(myState$size$data, nodeSet = "organisations",
                                              addSame = TRUE)

effectsWithCovariates <- local({
  ns <- asNamespace("MoNAn")
  Filter(function(name) {
    f <- get(name, envir = ns)
    is.function(f) && "getTargetContribution" %in% names(formals(f)) &&
      any(c("attribute.index", "resource.attribute.index") %in% names(formals(f)))
  }, getNamespaceExports(ns))
})

# a correct specification of each of these effects, and the inputs of the
# wrong kind for each covariate argument
correctArgs <- function(effect) {
  formalNames <- names(formals(get(effect, envir = asNamespace("MoNAn"))))
  args <- list()
  if ("attribute.index" %in% formalNames) {
    args$attribute.index <- if (attributeKinds[[effect]] == "dyadic") "sameRegion" else "size"
  }
  if ("resource.attribute.index" %in% formalNames) {
    args$resource.attribute.index <- if (grepl("_cont|additional_origin", effect)) "second_or" else "sex"
  }
  args
}
wrongInputs <- list(
  attribute.index = function(effect) {
    other <- if (attributeKinds[[effect]] == "dyadic") "size" else "sameRegion"
    c("organisations", "people", "sex", "transfers", "sizee", "dep.var", other)
  },
  resource.attribute.index = function(effect) {
    c("organisations", "people", "size", "sameRegion", "transfers", "sexx")
  }
)

test_that("every effect with a node attribute declares the kind of covariate it reads", {
  withAttribute <- Filter(function(name) {
    "attribute.index" %in% names(formals(get(name, envir = asNamespace("MoNAn"))))
  }, effectsWithCovariates)
  expect_setequal(withAttribute, names(attributeKinds))
})

test_that("effects accept covariates of the right kind", {
  for (effect in effectsWithCovariates) {
    spec <- c(list(effect), correctArgs(effect))
    if (effect == "same_covariate") spec$attribute.index <- "sizeSame"
    expect_no_error(createEffectsObject(list(spec), checkProcessState = covariateState))
  }
})

test_that("effects reject nodesets, the dependent variable and covariates of the wrong kind", {
  for (effect in effectsWithCovariates) {
    args <- correctArgs(effect)
    for (arg in names(args)) {
      for (value in wrongInputs[[arg]](effect)) {
        spec <- c(list(effect), args)
        spec[[arg]] <- value
        expect_error(createEffectsObject(list(spec), checkProcessState = covariateState),
                     label = paste(effect, arg, "=", value))
      }
    }
  }
})

test_that("a misspelled covariate is named in the error", {
  expect_error(createEffectsObject(list(list("alter_covariate", attribute.index = "sizee")),
                                   checkProcessState = covariateState),
               paste("Effect alter_covariate: node.attribute must be a covariate that applies",
                     "to the nodeset 'organisations'. Input that caused the error: 'sizee'"),
               fixed = TRUE)
})

test_that("covariates are checked again when estimation and simulation start", {
  # effects created for a state, used with a state where the covariate is a
  # network
  effects <- createEffects(myState) |> addEffect(alter_covariate, node.attribute = "size")
  badState <- myState
  badState$size <- createNetwork(myState$sameRegion$data, nodeSet = c("organisations", "organisations"))
  alg <- createAlgorithm(badState, effects)
  expect_error(estimateMobilityNetwork(badState, effects, alg),
               "must be a covariate that applies to the nodeset 'organisations'")
  expect_error(simulateMobilityNetworks(badState, effects, parameters = 0, allowLoops = TRUE,
                                        burnin = 10, thinning = 10, nSimulations = 1),
               "must be a covariate that applies to the nodeset 'organisations'")
})

test_that("the C++ effects check the size of their covariates", {
  # effects created without a state are not checked in R; the C++ effects
  # stop instead of reading outside the data
  for (spec in list(list("dyadic_covariate", attribute.index = "organisations"),
                    list("dyadic_covariate", attribute.index = "size"),
                    list("alter_covariate", attribute.index = "sameRegion"),
                    list("loops_node_covar", attribute.index = "people"))) {
    effects <- createEffectsObject(list(spec))
    expect_error(cppNetworkStatistics(myState$dep.var, myState, myState$cache, effects),
                 "attribute.index in", label = paste(spec, collapse = " "))
  }
})

test_that("the elements of a state have unique names", {
  expect_error(createProcessState(list(transfers = myState$transfers, people = myState$people,
                                       organisations = myState$organisations,
                                       organisations = myState$sameRegion),
                                  "transfers"),
               "must be unique; used more than once: organisations")
})
