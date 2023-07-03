# This script runs a simple example with the data from the MoNAn package

library(MoNAn)

# packages for parallel computing
library(snowfall)

##### create data objects from internal data files, which are later combined to the process state #####

# create objects
people <- createNodeSet(nrow(mobilityEdgelist))
organisations <- createNodeSet(length(orgRegion))
transfers <- createEdgelist(mobilityEdgelist, nodeSet = c("organisations", "organisations", "people"))
sameRegion <- outer(orgRegion, orgRegion, "==") * 1
sameRegion <- createNetwork(sameRegion, nodeSet = c("organisations", "organisations"))
region <- createNodeVariable(orgRegion, nodeSet = "organisations")
size <- createNodeVariable(orgSize, nodeSet = "organisations", addSim = TRUE)
sex <- createNodeVariable(indSex, nodeSet = "people")

# combine created objects to the process state
myState <- createProcessState(list(
  transfers = transfers,
  
  people = people,
  organisations = organisations,
  
  sameRegion = sameRegion,
  region = region,
  size = size,
  sex = sex
))


##### create cache #####

# define dependent variable and create cache object
myDependentVariable <- "transfers"
myCache <- createWeightedCache(myState, myDependentVariable, resourceCovariates = c("sex"))


##### create effects object #####

# create an effects object
myEffects <- createEffectsObject(
  list(
    list("loops"),
    list("min_reciprocity"),
    list("dyadic_covariate", attribute.index = "sameRegion"),
    list("alter_covariate", attribute.index = "size"),
    list("resource_covar_to_node_covar", attribute.index = "region", resource.attribute.index = "sex"),
    list("loops_resource_covar", resource.attribute.index = "sex")
  )
)


##### get multinomial statistics to estimate initial parameters using pseudo-likelihood estimation #####

# create statistics object, to be used, e.g., with the mlogit package
myStatisticsFrame <- getMultinomialStatistics(myState, myCache, myEffects, myDependentVariable)

### additional script to get pseudo-likelihood estimates
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



##### estimate mobility network model #####

# estimate mobility network model
myResDN <- estimateMobilityNetwork(myDependentVariable,
  myState, myCache, myEffects,
  initialParameters = NULL,
  # in case a pseudo-likelihood estimation was run, replace with
  # initialParameters = initEst,
  burnInN1 = 1500, iterationsN1 = 50, thinningN1 = 750, gainN1 = 0.1,
  burnInN2 = 7500, nsubN2 = 4, initGain = 0.2, thinningN2 = 1500,
  initialIterationsN2 = 25,
  iterationsN3 = 500, burnInN3 = 7500, thinningN3 = 3750,
  parallel = T, cpus = 4,
  allowLoops = T,
  verbose = T,
  returnDeps = T,
  multinomialProposal = F,
  fish = F
)

# check convergence
max(abs(myResDN$convergenceStatistics))

# estimate mobility network model again based on previous results to improve convergence
myResDN <- estimateMobilityNetwork(myDependentVariable,
  myState, myCache, myEffects,
  prevAns = myResDN,
  burnInN1 = 1500, iterationsN1 = 50, thinningN1 = 750, gainN1 = 0.1,
  burnInN2 = 7500, nsubN2 = 4, initGain = 0.2, thinningN2 = 3000,
  initialIterationsN2 = 40,
  iterationsN3 = 500, burnInN3 = 15000, thinningN3 = 7500,
  parallel = T, cpus = 4,
  allowLoops = T,
  verbose = T,
  returnDeps = T,
  multinomialProposal = F,
  fish = F
)

# check convergence
max(abs(myResDN$convergenceStatistics))

# view results
myResDN


##### regression diagnostics #####

autoCorrelationTest(myDependentVariable, myResDN)

traces <- extractTraces(myDependentVariable, myResDN, myEffects)
plot(traces)


##### test whether other effects should be included #####

myEffects2 <- createEffectsObject(
  list(
    list("loops"),
    list("min_reciprocity"),
    list("dyadic_covariate", attribute.index = "sameRegion"),
    list("alter_covariate", attribute.index = "size"),
    list("resource_covar_to_node_covar", attribute.index = "region", resource.attribute.index = "sex"),
    list("loops_resource_covar", resource.attribute.index = "sex"),
    list("min_transitivity")
  )
)

test_ME.2 <- scoreTest(myDependentVariable, myResDN, myEffects2)
test_ME.2


##### goodness of fit #####

myGofIndegree <- gofDistributionNetwork(ans = myResDN, simulations = myResDN$deps, gofFunction = getIndegree, lvls = 1:100)
plot(myGofIndegree,  lvls = 20:70)

myGofTieWeight <- gofDistributionNetwork(ans = myResDN, simulations = myResDN$deps, gofFunction = getTieWeights, lvls = 1:30)
plot(myGofTieWeight, lvls = 1:15)
