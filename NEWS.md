# MoNAn2 (development version)

!! The package has been renamed from MoNAn to MoNAn2 !!

Load it with `library(MoNAn2)`. Saved objects that contain effect functions
(e.g. effects objects created with MoNAn) refer to the MoNAn namespace and
should be recreated with MoNAn2.

## Bug fixes that change results

Results from models that used the affected effects or options should be re-run.

- `concentration_prop`: the change statistic did not account for the change in
  the origin's row sum when a resource moves, so it did not match the target
  statistic. It is now calculated over the whole row of the origin.
  The target statistic is unchanged.
- `reciprocity_min_resource_covar`: the change statistic was +/-1 even when the
  moving resource did not have the attribute. It now uses the resource's value
  on the attribute, which also makes it correct for non-binary attributes.
  The target statistic is unchanged.
- `multinomialProposal = TRUE`: the statistics for the candidate destinations
  were calculated with the cache from before the resource was removed, which
  gave wrong proposal probabilities for all non-linear effects.
- `simulateMobilityNetworks`: `allowLoops` was ignored after the burn-in
  (loops were never allowed).
- Parallel estimation (`parallel = TRUE`): `allowLoops` was ignored in Phase 1.
  This only affected starting values and step sizes, not the final estimates.

## Other bug fixes

- `avoiding_dissimilar_covar_cont` and
  `joining_similar_avoiding_dissimilar_covar_cont` only worked if the
  dependent variable was called "transfers".

- `loops_x_loops_additional_origin` stopped estimation and simulation with
  "missing value where TRUE/FALSE needed" when a resource whose additional
  origin is its origin was proposed to move into a loop: it read the origin
  from the edgelist, where the sampler sets the moving resource to NA.

- `concentration_prop_AC` could not be used: `alpha` was not an argument of the
  function. `alpha` (default 2) has been added, and the change statistic is
  calculated over the whole row of the origin, as in `concentration_prop`.
- The startup message reads the version of the attached package instead of
  the hard-coded name "MoNAn".
- Parallel workers load MoNAn2 instead of MoNAn.
- The example data `myEffects` has been recreated with MoNAn2, as it referred
  to the MoNAn namespace.

## New features

- Phase 2 can use Geyer-Thompson MCMC maximum likelihood estimation instead of
  the Robbins-Monro algorithm: `createAlgorithm(..., phase2Method = "GT")`,
  with the new arguments `iterationsGT`, `maxIterationsGT` and
  `finalIterationsGT`. Starting values are maximum pseudo-likelihood
  estimates, so Phase 1 is not needed. Robbins-Monro remains the default, and
  algorithm objects created with earlier versions use it.
- If neither `initialParameters` nor `prevAns` are given,
  `estimateMobilityNetwork` now starts from maximum pseudo-likelihood
  estimates (the multinomial logit model for which `getMultinomialStatistics`
  provides the data) instead of zero, for both Robbins-Monro and
  Geyer-Thompson estimation. With Robbins-Monro, Phase 1 then only determines
  the step sizes. If the estimates do not exist, estimation starts from zero
  as before.

## Performance

- The MCMC sampler runs in C++ for models whose effects all have a C++
  version, currently `loops`, `loops_AC`, `loops_GW`, `loops_node_covar`,
  `loops_resource_covar`, `loops_resource_covar_node_covar`,
  `loops_additional_origin`, `loops_x_loops_additional_origin`,
  `concentration_AC`, `concentration_GW`, `concentration_AC_dyad_covar`,
  `concentration_GW_dyad_covar`, `concentration_AC_resource_covar_bin`,
  `concentration_GW_resource_covar_bin`, `reciprocity_AC`, `reciprocity_GW`,
  `reciprocity_AC_dyad_covar`, `reciprocity_AC_dyad_covar_bin`,
  `reciprocity_GW_dyad_covar`, `reciprocity_GW_dyad_covar_bin`,
  `in_weights_AC`, `in_weights_GW`, `transitivity_basic`, `transitivity_min`,
  `transitivity_netflow`, `transitivity_AC`, `transitivity_GW`, `triad120C`,
  `triad120D`, `triad120U`, `dyadic_covariate`, `alter_covariate` and
  `resource_covar_to_node_covar`. A sampler step is about 40-50 times faster; the
  estimation of the example model with the default algorithm takes 17 seconds.
  For the same seed, results are identical to those of the R sampler, which is
  used for all other models. `options(MoNAn2.useCpp = FALSE)` switches to the
  R sampler. The package now needs a C++ compiler when installed from source.

- The MCMC sampler updates the edgelist and the cache in place instead of
  copying them in every simulation step. The time per step no longer grows
  with the number of locations and resources (e.g. about 5 times faster with
  300 locations and 50,000 resources); for small data it is unchanged.
  Results are identical for the same seed. Objects returned to users (results,
  simulations, stored draws) are plain lists as before.
- Faster effects, with identical results: the transitivity and triad effects
  (`transitivity_basic`, `transitivity_min`, `transitivity_netflow`,
  `transitivity_AC`, `transitivity_GW`, `triad120C`, `triad120D`, `triad120U`)
  no longer copy the network in every call and are vectorised (e.g.
  `transitivity_GW` is about 25 times faster). The checks that covariates are
  symmetric or binary (reciprocity effects with dyadic covariates,
  `crowding_out_prop_covar_bin`, `staying_by_prop_bin_inflow`) now run once
  when the effect is added and when estimation or simulation starts, instead
  of in every call.

## Removed

- `in_weights_exponent` has been commented out and is no longer exported: its
  change and target statistics do not match (the target statistic excludes
  loops and divides by the number of nodes, the change statistic does neither).

## Documentation

- `monanDataCreate`: the `fixedEffectDummies` argument is now documented (the
  documentation used the misspelled name `fixedEffectsDummies`) and explains
  that the dummies are named "dummy2", "dummy3", etc. and can be added with
  `addDummyEffects()` / `addFixedEffects()`. Typos in the description and
  example fixed.
- Removed outdated help pages of effects that are no longer part of the package.

## Tests

- New tests check for every exported effect that the change statistic equals
  the difference in target statistics before and after a move (including moves
  into and out of loops), and that every exported effect is tested.
  Set the environment variable `MONAN_FULL_EFFECT_TESTS=true` to test every
  resource instead of a sample.

# MoNAn 1.1.1

Bug fixes if only one effect is included in an analysis.

Inclusion of option to add fixed effects for destination by default.

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

Implementation of paralellization in Phase 1 of the estimation.

# MoNAn 0.1.3

Various bugfixes and improvements.

Implementation of new effects.

User-friendly checks in creation of MoNAn objects.

# MoNAn 0.1.2

Original Submission
