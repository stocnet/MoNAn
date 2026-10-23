# MoNAn (development version)

## Performance

- The maximum pseudo-likelihood estimates (starting values) compute the change
  statistics once instead of in each Newton step, and store them only once for
  resources with the same origin, destination and values of the resource
  covariates used by the effects (about 5 times faster for the example data).
  If this would need more than `getOption("MoNAn.mpleMemoryLimit")` bytes
  (default 1 GB), the change statistics are recalculated in each step as
  before.

# MoNAn 2.1.1

## Bug fixes

- The tests comparing C++ and R effects allow rounding differences on R built
  without `long double` (CRAN's "noLD" check), where R sums in double
  precision; elsewhere they still require identical values.

# MoNAn 2.1.0

## Breaking changes / changes so that previous code produces different results

- `scoreTest` now performs a score test that is conditional on the statistics
  of the effects in the model: the deviation of a tested effect's statistic is
  adjusted for the remaining deviations of the model's statistics and divided
  by its standard deviation given the model's statistics. The previous
  parametric p-values ignored this correlation and were far too conservative; 
  the new p-values seem better calibrated in simple simulations.
  The result reports `tStatistics`, `pValues` and `R2` (the
  share of the tested statistic's variance explained by the model's
  statistics); the previous non-parametric value is available with
  `nonParametric = TRUE`.
- `estimateMobilityNetwork(..., returnDeps = TRUE)` stores the statistics of
  the model's effects in the observed and simulated networks
  (`observedStatistics`, `simulatedStatistics`), which `scoreTest` uses. For
  results estimated with earlier versions, specify the effects of the
  estimated model with `scoreTest(..., modelEffects = )`.
- The example results `myResDN` were re-estimated with these statistics.

## New features

- Parallel estimation (`parallel = TRUE`) uses the package 'parallel' (part of
  R) instead of 'snowfall', which is no longer a dependency. One cluster is
  started per estimation instead of one per phase. Each chain gets its own
  random number stream derived from the random numbers of the
  R session, so parallel estimations are now reproducible with `set.seed()`
  and the chains are independent.

- The concentration effects `concentration_basic_squared`,
  `concentration_norm`, `concentration_norm_squared`,
  `concentration_prop_orig_cov` and `concentration_rankGW`, removed in 2.0.0,
  are available again, unchanged from version 1.1.0, now also in C++.
  `concentration_norm` and `concentration_norm_squared` divide by the number
  of resources of the origin minus one, which is zero for origins with only one
  resource.

## Bug fixes

- `getIndegree` and `getTieWeights` failed or put counts in the wrong places
  when the levels in `gofMobilityNetwork(..., lvls = )` did not start at 1.
  Without `lvls`, `getIndegree` failed and `getTieWeights` returned
  distributions of different lengths for different networks, which `plot`
  misaligned; now all networks get the same levels, up to the largest value.
  `plot` of goodness of fit results now also works without `lvls`. Results for
  levels starting at 1 (e.g. `lvls = 1:100`) are unchanged.
- Plots of goodness of fit results label the x-axis with the values the
  counts refer to (e.g., indegree 0, 1, 2, ...) instead of positions:
  `getIndegree` and `getTieWeights` now name their results by these values.
  `plot` passes further arguments (e.g. `main`, `xlab`) to `boxplot`.

## Performance

- `extractTraces` uses the statistics stored by estimations with
  `returnDeps = TRUE` instead of recalculating them from the simulations.

- The maximum pseudo-likelihood estimates used as starting values are
  calculated in C++ if all effects have a C++ version (about 4-5 times faster
  for larger data), without storing the change statistics of all resources and
  destinations, so memory use does not grow with the data. The estimates are
  the same as with the R version, which is used otherwise.

# MoNAn 2.0.0

Version 2 runs the simulations in C++, which makes estimation typically 20-50
times faster, adds Geyer-Thompson MCMC maximum likelihood estimation, and
corrects several effects. Some changes are not backward compatible (see
"Breaking changes"), and results of models with corrected effects differ from
version 1.

## Breaking changes / changes so that previous code produces different results

- Removed effects:
  - `in_weights_exponent`: its change and target statistics did not match
    (the target statistic excludes loops and divides by the number of nodes,
    the change statistic does neither).
  - `concentration_basic_squared`, `concentration_norm`,
    `concentration_norm_squared`, `concentration_prop_orig_cov` and
    `concentration_rankGW`: older concentration effects that were only used
    for testing.
- `concentration_prop` has a new definition: each resource contributes the
  proportion of the resources from its origin that go to the same destination.
  The statistic of version 1 is no longer available.
- Corrections that change results; models that use them should be re-run:
  - `reciprocity_min_resource_covar`: the change statistic was +/-1 even when
    the moving resource did not have the attribute. It now uses the resource's
    value on the attribute, which also makes it correct for non-binary
    attributes.
  - `multinomialProposal = TRUE`: the statistics for the candidate destinations
    were calculated with the cache from before the resource was removed, which
    gave wrong proposal probabilities for all non-linear effects.
  - `simulateMobilityNetworks`: `allowLoops` was ignored after the burn-in
    (loops were never allowed).
  - Parallel estimation (`parallel = TRUE`): `allowLoops` was ignored in
    Phase 1. This only affected starting values and step sizes, not the final
    estimates.
- If neither `initialParameters` nor `prevAns` are given,
  `estimateMobilityNetwork` starts from maximum pseudo-likelihood estimates
  instead of zero (see below), so results for the same seed differ from
  version 1 (within Monte Carlo error).
- `simulateMobilityNetworks` stops if the number of parameters differs from
  the number of effects; before, parameters were silently recycled.
- The package contains compiled C++ code: installing it from source requires
  a C++ compiler.

## New features

- Phase 2 can use Geyer-Thompson MCMC maximum likelihood estimation instead of
  the Robbins-Monro algorithm: `createAlgorithm(..., phase2Method = "GT")`,
  with the new arguments `iterationsGT`, `maxIterationsGT` and
  `finalIterationsGT`. Robbins-Monro remains the default, and algorithm
  objects created with earlier versions use it.
  Note that Geyer-Thompson still uses a "normal" Phase 3 that calculates standard 
  errors and checks convergence.
- Maximum pseudo-likelihood estimates (the multinomial logit model for which
  `getMultinomialStatistics` provides the data) are used as starting values if
  neither `initialParameters` nor `prevAns` are given, for both Robbins-Monro
  and Geyer-Thompson estimation. With Robbins-Monro, Phase 1 then only
  determines the step sizes; with Geyer-Thompson, Phase 1 is not needed. If
  the estimates do not exist, estimation starts from zero as before.
- Fixed effects for locations: `monanDataCreate(..., fixedEffectDummies = TRUE)`
  adds a dummy variable for each location (except the first) to the process
  state, and `addDummyEffects()` / `addFixedEffects()` add them to the model.
- New effects: `concentration_max`, `concentration_basic_cube` and
  `concentration_prop_AC`.
- `options(MoNAn.useCpp = FALSE)` switches to the sampler written in R, which
  gives identical results for the same seed and is also used for models with
  effects written only in R.

## Performance

- The MCMC sampler and all effects run in C++. A sampler step is about 20-50
  times faster, depending on the effects; the estimation of the example model
  with the default algorithm takes 17 seconds instead of 13 minutes.
- The sampler updates the edgelist and the cache in place instead of copying
  them in every simulation step, so the time per step no longer grows with the
  number of locations and resources.
- The transitivity and triad effects no longer copy the network in every call
  and are vectorised, and the checks that covariates are symmetric or binary
  run once when an effect is added and when estimation or simulation starts,
  instead of in every call.

## Other bug fixes

- Estimation of models with only one effect.
- `multinomialProposal = TRUE` stopped with an error ("'x' must be an array of
  at least two dimensions") in models with only one effect.
- `avoiding_dissimilar_covar_cont` and
  `joining_similar_avoiding_dissimilar_covar_cont` only worked if the
  dependent variable was called "transfers".
- `loops_x_loops_additional_origin` stopped estimation and simulation with
  "missing value where TRUE/FALSE needed" when a resource whose additional
  origin is its origin was proposed to move into a loop.
- The example of `simulateMobilityNetworks` used 6 parameters for 7 effects.

## Documentation

- `monanDataCreate`: the `fixedEffectDummies` argument is documented (the
  documentation used the misspelled name `fixedEffectsDummies`).
- Removed outdated help pages of effects that are no longer part of the
  package.
- For developers, `ADDING_EFFECTS.md` in the GitHub repository explains how to
  add an effect in R and C++.

## Tests

- Tests check for every effect that the change statistic equals the difference
  in target statistics before and after a move, and that the C++ effects and
  sampler give exactly the same results as their R versions.

# MoNAn 1.1.0

Added effects that improve on the old GW effects (where AC stands for alternating cliques):

- concentration_AC, concentration_AC_dyad_covar,
  concentration_AC_resource_covar_bin,
  reciprocity_AC, reciprocity_AC_dyad_covar_bin, reciprocity_AC_dyad_covar,
  transitivity_AC,
  loops_AC,
  in_weights_AC

# MoNAn 1.0.1

# MoNAn 1.0.0

!! Backward breaking change !!

The "cache" is not visible to the user anymore, it is created automatically
and there is no more use for the function createWeightedCache.

The following functions are used differently:

- gofMobilityNetwork does not require the specification of simulations anymore!

New functions (aliases) for easier use

- monanDependent (createEdgelist), monanNodes (createNodeSet), 
  monanEdges (createNodeSet), monadicCovariate (createNodeAttribute),
  dyadicCovariate (createNetwork), monanEstimate (estimateMobilityNetwork),
  monanSimulate (simulateMobilityNetwork), monanGOF (gofMobilityNetwork)

Further changes:

- createEffect and addEffects; specification of node.attribute and
  edge.attribute instead of attribute.index and resource.attribute.index

- The results objects are much smaller and contain less redundant information.

- The State and Effects objects have print functions for easier inspection.

- The algorithm is now stored in the outcome object for later reference.

# MoNAn 0.3.0

New functions to specify the model: createEffects and addEffect
New function to generate the process state: monanDataCreate

Creation of alias: monanAlgorithmCreate for createAlgorithm
Creation of alias: monan07 for estimateMobilityNetwork

Implementation of new effects

# MoNAn 0.2.0

Implementation of parallelization in Phase 1 of the estimation.

# MoNAn 0.1.3

Various bugfixes and improvements.

Implementation of new effects.

User-friendly checks in creation of MoNAn objects.

# MoNAn 0.1.2

Original Submission
