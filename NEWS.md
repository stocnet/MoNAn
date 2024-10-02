# MoNAn 1.1.1

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
