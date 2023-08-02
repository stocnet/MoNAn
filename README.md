
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MoNAn

<!-- badges: start -->
<!-- badges: end -->

MoNAn is the software implementation of the statistical model outlined
in:

Block, P., Stadtfeld, C., & Robins, G. (2022). A statistical model for
the analysis of mobility tables as weighted networks with an application
to faculty hiring networks. Social Networks, 68, 264-278.

[Link to the study in Social
Networks](https://www.sciencedirect.com/science/article/pii/S0378873321000654).

[Pre-print with minor differences to the journal
version](https://osf.io/preprints/socarxiv/n86rx/).

The model allows for the analysis of emergent structures (known from
network analysis) in mobility tables, alongside exogenous predictors. It
combined features from classical log-linear models for mobility tables
and Exponential Random Graph Models (ERGMs). In the absence of emergent
structures, the models reduces to the former, when ignoring
characteristics of individuals, it is an ERGM for weighted networks.

Announcements about workshops etc. can be found
[here](https://www.suz.uzh.ch/de/institut/professuren/block/monan_software.html).

# Warning!!!

The package is still under development, especially the documentation.
When using the package you might encounter bugs or errors, or you might
not be able to do what you want. In that case please write the package
maintainer under his institutional email address (requires using
google).

# Installation

The package is under development. You can install the current
development version of MoNAn from GitHub using:

``` r
# install.packages("remotes")
remotes::install_github("perblock/MoNAn")
```

# Example

In this section we outline a simple example with synthetic data stored
in the MoNAn package.

``` r
library(MoNAn)

# packages for parallel computing
library(snowfall)
#> Loading required package: snow
```

## The Data

The example case uses synthetic data that represents mobility of 742
individuals between 17 organisations. The artificial data includes an
edgelist, i.e. a list of origins (column 1) and destinations (column 2)
of all individuals sorted by origin.

``` r
mobilityEdgelist[1:10,]
#>       [,1] [,2]
#>  [1,]    1    1
#>  [2,]    1    1
#>  [3,]    1    2
#>  [4,]    1    2
#>  [5,]    1    2
#>  [6,]    1    2
#>  [7,]    1    3
#>  [8,]    1    3
#>  [9,]    1    3
#> [10,]    1    3
```

The data further includes (artificial) characteristics of individuals,
which we will call sex. However, using continuous covariates is also
possible, although for some effects this would not be meaningful.

``` r
indSex[1:10]
#>  [1] 1 1 0 0 0 0 1 0 1 0
```

Characteristics of the organisations are the region in which they are
located (binary, northern/southern island) and the organisations’ size
that represents a composite measure of the number of workers, assets,
and revenue.

``` r
orgRegion
#>  [1] 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0
orgSize
#>  [1]  3.202903  5.178353  7.592492  7.368718  9.098950  7.977469  7.751859
#>  [8] 10.454544  7.811348  4.958048  7.740106  6.112728  6.640621  2.850886
#> [15]  8.955254  7.597493  9.313139
```

## Working with the package

### Preparing the data

First, we create MoNAn data objects from the introduced data files. A
necessary part is that we define and name the nodesets, i.e., we name
who (“people”) is mobile between what (“organisations”).

``` r
# extract number of individuals and organisations from the mobility data

N_ind <- nrow(mobilityEdgelist)
N_org <- length(unique(as.numeric(mobilityEdgelist)))

# create monan objects
people <- createNodeSet(N_ind)
organisations <- createNodeSet(N_org)
transfers <- createEdgelist(mobilityEdgelist, nodeSet = c("organisations", "organisations", "people"))
sameRegion <- outer(orgRegion, orgRegion, "==") * 1
sameRegion <- createNetwork(sameRegion, nodeSet = c("organisations", "organisations"))
region <- createNodeVariable(orgRegion, nodeSet = "organisations")
size <- createNodeVariable(orgSize, nodeSet = "organisations", addSim = TRUE)
sex <- createNodeVariable(indSex, nodeSet = "people")
```

We combine the data objects into the process state, i.e., a MoNAn object
that stores all information about the data that will be used in the
estimation later. Here we also need to specify what will be the outcome
variable to be modelled in the analysis.

``` r
myState <- createProcessState(
  list(
    transfers = transfers,
    
    people = people,
    organisations = organisations,
    
    sameRegion = sameRegion,
    region = region,
    size = size,
    sex = sex
  ),
  dependentVariable = "transfers"
)
```

We create a cache (a necessary object used in the estimation of the
model). In case variables of the individuals in the data are included in
the state (here: “sex”), they need to be explicitly mentioned in the
creation of the cache under “resourceCovariates”.

``` r
myCache <- createWeightedCache(myState, resourceCovariates = c("sex"))
```

### Specifying the model

The predictors in the model are called “Effects” and they are defined in
a list. Each effect itself is a list that contains the effect name and
the additional parameters it needs.

``` r
# create an effects object
myEffects <- createEffectsObject(
  list(
    list("loops"),
    list("reciprocity_min"),
    list("dyadic_covariate", attribute.index = "sameRegion"),
    list("alter_covariate", attribute.index = "size"),
    list("resource_covar_to_node_covar", 
         attribute.index = "region", resource.attribute.index = "sex"),
    list("loops_resource_covar", resource.attribute.index = "sex")
  ),
  checkProcessState = myState
)
```

### Optional: Pre-estimation

We can run a pseudo-likelihood estimation that gives a (biased) guess of
the model results. We use this to get improved initial estimates, which
increases the chances of model convergence in the first run of the
estimation considerably. To get pseudo-likelihood estimates, we need to
use functions from other libraries to estimate a multinomial logit model
(e.g., “dfidx” and “mlogit”)

``` r
# create multinomial statistics object pseudo-likelihood estimation
myStatisticsFrame <- getMultinomialStatistics(myState, myCache, myEffects)

### additional script to get pseudo-likelihood estimates, requires the dfidx and mlogit package
# library(dfidx)
# library(mlogit)
# my.mlogit.dataframe <- dfidx(myStatisticsFrame,
#                           shape = "long", 
#                           choice = "choice")
# 
# colnames(my.mlogit.dataframe) <- gsub(" ", "_", colnames(my.mlogit.dataframe))
# 
# IVs <- (colnames(my.mlogit.dataframe)[2:(ncol(myStatisticsFrame)-2)])
# 
# f <- as.formula(paste("choice ~ 1 + ", paste(IVs, collapse = " + "), "| 0"))
# 
# my.mlogit.results <- mlogit(formula = eval(f), data = my.mlogit.dataframe, heterosc = F)
# 
# summary(my.mlogit.results)
#
# initEst <- my.mlogit.results$coefficients[1:length(IVs)]
```

### Specifying the estimation algorithm

The estimation algorithm requires multiple parameters that guide the
simulation of the chain. Most values are generated by default from the
state and the effects object.

``` r
myAlg <- createAlgorithm(myState, myEffects, multinomialProposal = FALSE)
```

### Estimation

Now we can estimate the model. The first two lines indicate the
dependent variable, data (state), cache, and effects. The third line
specifies the initial estimates, where the previously obtained
pseudo-likelihood estimates can be used. The remaining lines define
whether parallel computing should be used, how much console output is
desired during the estimation, and whether the simulations from phase 3
should be stored.

Running the model takes a while (up to 10 minutes for this data with
parallel computing).

``` r
myResDN <- estimateMobilityNetwork(
  myState, myCache, myEffects, myAlg,
  initialParameters = NULL,
  parallel = T, cpus = 4,
  verbose = T,
  returnDeps = T,
  fish = F
)
#> Starting phase 1
#> Starting burn-in with 12614 steps
#> ........................
#> Starting phase 2
#> R Version:  R version 4.3.1 (2023-06-16)
#> snowfall 1.84-6.2 initialized (using snow 0.4-4): parallel execution on 4 CPUs.
#> Library MoNAn loaded.
#> Library MoNAn loaded in cluster.
#> Starting sub phase 1 
#> New parameters:
#> loops 
#>  1.96239947914881 
#>  reciprocity_min 
#>  0.5921525719097 
#>  dyadic_covariate sameRegion 
#>  1.40049574029235 
#>  alter_covariate size 
#>  0.033007136338196 
#>  resource_covar_to_node_covar region sex 
#>  -0.552542320949847 
#>  loops_resource_covar sex 
#>  0.0129469964664311 
#> Starting sub phase 2 
#> New parameters:
#> loops 
#>  2.39723742802551 
#>  reciprocity_min 
#>  0.735171317280412 
#>  dyadic_covariate sameRegion 
#>  1.59953402045278 
#>  alter_covariate size 
#>  0.0341741994300065 
#>  resource_covar_to_node_covar region sex 
#>  -0.615722236162676 
#>  loops_resource_covar sex 
#>  -0.258349863937289 
#> Starting sub phase 3 
#> New parameters:
#> loops 
#>  2.52006155407795 
#>  reciprocity_min 
#>  0.782692703309321 
#>  dyadic_covariate sameRegion 
#>  1.66117156519656 
#>  alter_covariate size 
#>  0.0360495088741273 
#>  resource_covar_to_node_covar region sex 
#>  -0.6380053370706 
#>  loops_resource_covar sex 
#>  -0.322942579704397 
#> Starting sub phase 4 
#> New parameters:
#> loops 
#>  2.55831223642456 
#>  reciprocity_min 
#>  0.801052602801524 
#>  dyadic_covariate sameRegion 
#>  1.67900023443131 
#>  alter_covariate size 
#>  0.0378054612771427 
#>  resource_covar_to_node_covar region sex 
#>  -0.643234381003142 
#>  loops_resource_covar sex 
#>  -0.348498091145485
#> 
#> Stopping cluster
#> Starting phase 3
#> snowfall 1.84-6.2 initialized (using snow 0.4-4): parallel execution on 4 CPUs.
#> Library MoNAn loaded.
#> Library MoNAn loaded in cluster.
#> 
#> Stopping cluster
```

In case a pseudo-likelihood estimates have been obtained previously,
this can be specified by

``` r
myResDN <- estimateMobilityNetwork(
  myState, myCache, myEffects, myAlg,
  initialParameters = initEst,
  parallel = T, cpus = 4,
  verbose = T,
  returnDeps = T,
  fish = F
)
```

Check convergence to see whether the results are reliable. In case the
maximum convergence ratio is above 0.1 (or 0.2 for less precise
estimates), another run is necessary.

``` r
max(abs(myResDN$convergenceStatistics))
#> [1] 0.06165639
```

If convergence is too high, update algorithm, re-run estimation with
previous results as starting values and check convergence:

``` r
# estimate mobility network model again based on previous results to improve convergence
# with an adjusted algorithm
myAlg <- createAlgorithm(myState, myEffects, multinomialProposal = TRUE, 
                         initialIterationsN2 = 500, nsubN2 = 1, initGain = 0.02, iterationsN3 = 1000)

myResDN <- estimateMobilityNetwork(
  myState, myCache, myEffects, myAlg,
  prevAns = myResDN,
  parallel = T, cpus = 4,
  verbose = T,
  returnDeps = T,
  fish = F
)
#> skipping phase 1 and taking values from prevAns
#> Starting phase 2
#> snowfall 1.84-6.2 initialized (using snow 0.4-4): parallel execution on 4 CPUs.
#> Library MoNAn loaded.
#> Library MoNAn loaded in cluster.
#> Starting sub phase 1 
#> New parameters:
#> loops 
#>  2.57401586997664 
#>  reciprocity_min 
#>  0.830449011051864 
#>  dyadic_covariate sameRegion 
#>  1.69209883761208 
#>  alter_covariate size 
#>  0.0369637391202157 
#>  resource_covar_to_node_covar region sex 
#>  -0.650002430600979 
#>  loops_resource_covar sex 
#>  -0.360896544627875
#> 
#> Stopping cluster
#> Starting phase 3
#> snowfall 1.84-6.2 initialized (using snow 0.4-4): parallel execution on 4 CPUs.
#> Library MoNAn loaded.
#> Library MoNAn loaded in cluster.
#> 
#> Stopping cluster
```

``` r
# check convergence
max(abs(myResDN$convergenceStatistics))
#> [1] 0.09751736
```

In case convergence is still poor, updating the algorithm might be
necessary. Otherwise, we can view the results. The first column is the
estimate, the second column the standard error, and the third column the
convergence ratio. All values in the final column should be below 0.1
(see above).

``` r
myResDN
#> Results
#>                                   Effects   Estimates StandardErrors
#> 1                                   loops  2.57401587     0.18651741
#> 2                         reciprocity_min  0.83044901     0.17844344
#> 3             dyadic_covariate sameRegion  1.69209884     0.10998693
#> 4                    alter_covariate size  0.03696374     0.02426595
#> 5 resource_covar_to_node_covar region sex -0.65000243     0.17541006
#> 6                loops_resource_covar sex -0.36089654     0.22183745
#>    Convergence
#> 1 -0.005483438
#> 2 -0.044616791
#> 3 -0.097517364
#> 4  0.003780352
#> 5 -0.021552036
#> 6 -0.004111089
```

## Diagnostics of the estimated model

The following two diagnostics indicate the extent to which the chain
mixes (i.e., whether the thinning was chosen appropriately). The
autoCorrelationTest indicates the degree to which the values of the
dependent variable of consecutive draws from the chain in phase 3 are
correlated. Here lower values are better. Values above 0.5 are very
problematic and indicate that a higher thinning is needed.

``` r
autoCorrelationTest(myResDN)
#> [1] 0.2153447
```

The output of extractTraces indicates the correlation of statistics
between subsequent draws from the chain in phase 3. The plot should show
data point randomly scattered around the target line, as shown below. If
patterns in the traces are discernible, a higher thinning is needed.

``` r
traces <- extractTraces(myResDN, myEffects)

par(mfrow = c(1,2))
plot(traces)
```

<img src="man/figures/README-unnamed-chunk-19-1.png" width="100%" /><img src="man/figures/README-unnamed-chunk-19-2.png" width="100%" /><img src="man/figures/README-unnamed-chunk-19-3.png" width="100%" />

## Score-tests to check model specification

Based on an estimated model, a score-type test is available that shows
whether statistics representing non-included effects are well
represented. If this is not the case, it is likely that including them
will result in significant estimates.

``` r
myEffects2 <- createEffectsObject(
  list(
    list("transitivity_min")
  )
)

test_ME.2 <- scoreTest(myResDN, myEffects2)
test_ME.2
#> Results
#>            Effects pValuesParametric pValuesNonParametric
#> 1 transitivity_min      4.727783e-09                    0
#> 
#>  Parametric p-values: small = more significant 
#>  Non-parametric p-values: further away from 0.5 = more significant
```

The interpretation is that there appears to be some transitive
clustering in the data that the model does not account for in its
current form.

## GOF testing

Akin to ERGMs, goodness of fit testing is available to see whether
auxiliary statistics are well captured by the model.

``` r
myGofIndegree <- gofDistributionNetwork(ans = myResDN, simulations = myResDN$deps, gofFunction = getIndegree, lvls = 1:70)
plot(myGofIndegree, lvls = 1:70)
```

<img src="man/figures/README-unnamed-chunk-21-1.png" width="100%" />

``` r

myGofTieWeight <- gofDistributionNetwork(ans = myResDN, simulations = myResDN$deps, gofFunction = getTieWeights, lvls = 1:20)
plot(myGofTieWeight, lvls = 1:20)
```

<img src="man/figures/README-unnamed-chunk-21-2.png" width="100%" />

## Simulation

The package also provides the possibility to exemplarily simulate
mobility networks based on the data and the specified effects and
parameters.

``` r
mySimDN <- simulateMobilityNetworks(
  myState,
  myCache,
  myEffects,
  parameters = c(2, 1, 1.5, 0.1, -1, -0.5),
  allowLoops = TRUE,
  burnin = 45000,
  thinning = 15000,
  nSimulations = 10
)

mySimDN[[1]]
```
