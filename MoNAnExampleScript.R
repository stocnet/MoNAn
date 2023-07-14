# This script runs a simple example with the data from the MoNAn package

library(MoNAn)


##### create data objects from internal data files, which are later combined to the process state #####

# Create a process state out of the mobility data objects:
# create objects (which are later combined to the process state)
transfers <- createEdgelist(mobilityEdgelist,
                            nodeSet = c("organisations", "organisations", "people")
)
people <- createNodeSet(1:nrow(mobilityEdgelist))
organisations <- createNodeSet(1:length(orgRegion))
sameRegion <- outer(orgRegion, orgRegion, "==") * 1
sameRegion <- createNetwork(sameRegion,
                            nodeSet = c("organisations", "organisations")
)
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
    list("reciprocity_min"),
    list("dyadic_covariate", attribute.index = "sameRegion"),
    list("alter_covariate", attribute.index = "size"),
    list("resource_covar_to_node_covar",
         attribute.index = "region",
         resource.attribute.index = "sex"
    ),
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


##### create algorithm object #####

# define algorithm based on state and effects characteristics
myAlg <- createAlgorithm(myDependentVariable, myState, myEffects, multinomialProposal = FALSE)


##### estimate mobility network model #####

# estimate mobility network model
myResDN <- estimateMobilityNetwork(myDependentVariable,
  myState, myCache, myEffects, myAlg,
  initialParameters = NULL,
  # in case a pseudo-likelihood estimation was run, replace with
  # initialParameters = initEst,
  parallel = T, cpus = 4,
  verbose = T,
  returnDeps = T,
  fish = F
)

# check convergence
max(abs(myResDN$convergenceStatistics))

myResDN_old <- myResDN

# estimate mobility network model again based on previous results to improve convergence
# with an adjusted algorithm
myAlg <- createAlgorithm(myDependentVariable, myState, myEffects, multinomialProposal = TRUE, 
                         initialIterationsN2 = 500, nsubN2 = 1, initGain = 0.02, iterationsN3 = 1000)

myResDN <- estimateMobilityNetwork(myDependentVariable,
  myState, myCache, myEffects, myAlg,
  prevAns = myResDN,
  parallel = T, cpus = 4,
  verbose = T,
  returnDeps = T,
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
    list("reciprocity_min"),
    list("dyadic_covariate", attribute.index = "sameRegion"),
    list("alter_covariate", attribute.index = "size"),
    list("resource_covar_to_node_covar", attribute.index = "region", resource.attribute.index = "sex"),
    list("loops_resource_covar", resource.attribute.index = "sex"),
    list("transitivity_min")
  )
)

test_ME.2 <- scoreTest(myDependentVariable, myResDN, myEffects2)
test_ME.2


##### goodness of fit #####

myGofIndegree <- gofDistributionNetwork(ans = myResDN, simulations = myResDN$deps, gofFunction = getIndegree, lvls = 1:100)
plot(myGofIndegree,  lvls = 20:70)

myGofTieWeight <- gofDistributionNetwork(ans = myResDN, simulations = myResDN$deps, gofFunction = getTieWeights, lvls = 1:30)
plot(myGofTieWeight, lvls = 1:15)
