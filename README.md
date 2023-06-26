
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MoNAn

<!-- badges: start -->
<!-- badges: end -->

MoNAn is the software implementation of the statistical model outlined
in:

Block, P., Stadtfeld, C., & Robins, G. (2022). A statistical model for
the analysis of mobility tables as weighted networks with an application
to faculty hiring networks. Social Networks, 68, 264-278.

as found
[here](https://www.sciencedirect.com/science/article/pii/S0378873321000654).

## Installation

The package is still under development. You can install the development
version of MoNAn from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("perblock/MoNAn")
```

or using:

``` r
# install.packages("remotes")
remotes::install_github("perblock/MoNAn")
```

## Example

# This script runs a simple example with the data from the MoNAn package

library(MoNAn)

##### create data objects from internal data files, which are later combined to the process state

# create objects

transfers \<- createEdgelist(mobilityEdgelist, nodeSet =
c(“organisations”, “organisations”, “people”)) people \<-
createNodeSet(1:nrow(mobilityEdgelist)) organisations \<-
createNodeSet(1:length(orgRegion)) sameRegion \<- outer(orgRegion,
orgRegion, “==”) \* 1 sameRegion \<- createNetwork(sameRegion, nodeSet =
c(“organisations”, “organisations”)) region \<-
createNodeVariable(orgRegion, nodeSet = “organisations”) size \<-
createNodeVariable(orgSize, nodeSet = “organisations”, addSim = TRUE)
sex \<- createNodeVariable(indSex, nodeSet = “people”)

# combine created objects to the process state

myState \<- createProcessState(list( transfers = transfers,

people = people, organisations = organisations,

sameRegion = sameRegion, region = region, size = size, sex = sex ))

##### create cache

# define dependent variable and create cache object

myDependentVariable \<- “transfers” myCache \<-
createWeightedCache(myState, myDependentVariable, resourceCovariates =
c(“sex”))

##### create effects object

# create an effects object

myEffects \<- createEffectsObject( list( list(“loops”),
list(“min_reciprocity”), list(“dyadic_covariate”, attribute.index =
“sameRegion”), list(“alter_covariate”, attribute.index = “size”),
list(“resource_covar_to_node_covar”, attribute.index = “region”,
resource.attribute.index = “sex”), list(“loops_resource_covar”,
resource.attribute.index = “sex”) ) )

##### get multinomial statistics to estimate initial parameters using pseudo-likelihood estimation

# create statistics object, to be used, e.g., with the mlogit package

myStatisticsFrame \<- getMultinomialStatistics(myState, myCache,
myEffects, myDependentVariable)

##### estimate mobility network model

# estimate mobility network model

myResDN \<- estimateMobilityNetwork(myDependentVariable, myState,
myCache, myEffects, initialParameters = NULL, burnInN1 = 1500,
iterationsN1 = 50, thinningN1 = 750, gainN1 = 0.1, burnInN2 = 7500,
nsubN2 = 4, initGain = 0.2, thinningN2 = 1500, initialIterationsN2 = 25,
iterationsN3 = 500, burnInN3 = 7500, thinningN3 = 3750, parallel = T,
cpus = 4, allowLoops = T, verbose = T, returnDeps = T,
multinomialProposal = F, fish = F )

# check convergence

max(abs(myResDN\$convergenceStatistics))

# estimate mobility network model again based on previous results to improve convergence

myResDN \<- estimateMobilityNetwork(myDependentVariable, myState,
myCache, myEffects, prevAns = myResDN, burnInN1 = 1500, iterationsN1 =
50, thinningN1 = 750, gainN1 = 0.1, burnInN2 = 7500, nsubN2 = 4,
initGain = 0.2, thinningN2 = 3000, initialIterationsN2 = 40,
iterationsN3 = 500, burnInN3 = 15000, thinningN3 = 7500, parallel = T,
cpus = 4, allowLoops = T, verbose = T, returnDeps = T,
multinomialProposal = F, fish = F )

# chack convergence

max(abs(myResDN\$convergenceStatistics))

# view results

myResDN

##### regression diagnostics

autoCorrelationTest(myDependentVariable, myResDN)

traces \<- extractTraces(myDependentVariable, myResDN, myEffects)
plot(traces)

##### test whether other effects should be included

myEffects2 \<- createEffectsObject( list( list(“loops”),
list(“min_reciprocity”), list(“dyadic_covariate”, attribute.index =
“sameRegion”), list(“alter_covariate”, attribute.index = “size”),
list(“resource_covar_to_node_covar”, attribute.index = “region”,
resource.attribute.index = “sex”), list(“loops_resource_covar”,
resource.attribute.index = “sex”), list(“min_transitivity”) ) )

test_ME.2 \<- scoreTest(myDependentVariable, myResDN, myEffects2)
test_ME.2

##### goodness of fit

myGofIndegree \<- gofDistributionNetwork(ans = myResDN, simulations =
myResDN\$deps, gofFunction = getIndegree, lvls = 1:100)
plot(myGofIndegree)

myGofTieWeight \<- gofDistributionNetwork(ans = myResDN, simulations =
myResDN\$deps, gofFunction = getTieWeights, lvls = 1:30)
plot(myGofTieWeight)

This is a basic example which shows you how to solve a common problem:

``` r
library(MoNAn)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
