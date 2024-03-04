# This script runs a simple example with the data from the MoNAn package

library(MoNAn)


##### create data objects from internal data files, which are later combined to the process state #####

# extract number of individuals and organisations from the mobility data

N_ind <- nrow(mobilityEdgelist)
N_org <- length(unique(as.numeric(mobilityEdgelist)))

# Create a process state out of the mobility data objects:
# create objects (which are later combined to the process state)
transfers <- createEdgelist(mobilityEdgelist,
                            nodeSet = c("organisations", "organisations", "people")
)
people <- createNodeSet(N_ind)
organisations <- createNodeSet(N_org)
sameRegion <- outer(orgRegion, orgRegion, "==") * 1
sameRegion <- createNetwork(sameRegion,
                            nodeSet = c("organisations", "organisations")
)
region <- createNodeVariable(orgRegion, nodeSet = "organisations")
size <- createNodeVariable(orgSize, nodeSet = "organisations", addSim = TRUE)
sex <- createNodeVariable(indSex, nodeSet = "people")

# the following lines create an artificial second origin used for illustration
# in the examples in the manual
other_origin <- sample(1:17, 742, replace = T)
resample <- as.logical(sample(0:1, 742, replace = T, prob = c(0.88, 0.12)))
other_origin[resample] <- transfers$data[resample,2]

second_or <- createNodeVariable(other_origin, nodeSet = "people")

# combine created objects to the process state
myState <- monanDataCreate(transfers,
                           people,
                           organisations,
                           sameRegion,
                           region,
                           size,
                           sex,
                           second_or)


##### create cache #####

# cache object
myCache <- createWeightedCache(myState, resourceCovariates = c("sex"))


##### create effects object #####

# effects object
myEffects <- createEffects(myState) |>
  addEffect(loops) |>
  addEffect(reciprocity_min) |>
  addEffect(dyadic_covariate, attribute.index = "sameRegion") |>
  addEffect(alter_covariate, attribute.index = "size") |>
  addEffect(resource_covar_to_node_covar,
     attribute.index = "region",
     resource.attribute.index = "sex") |>
  addEffect(loops_resource_covar, resource.attribute.index = "sex")


##### get multinomial statistics to estimate initial parameters using pseudo-likelihood estimation #####

# create statistics object, to be used, e.g., with the mlogit package
myStatisticsFrame <- getMultinomialStatistics(myState, myCache, myEffects)

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
# form <- as.formula(paste("choice ~ 1 + ", paste(IVs, collapse = " + "), "| 0"))
# 
# my.mlogit.results <- mlogit(formula = eval(form), data = my.mlogit.dataframe, heterosc = FALSE)
# 
# summary(my.mlogit.results)
# 
# initEst <- my.mlogit.results$coefficients[1:length(IVs)]


##### create algorithm object #####

# define algorithm based on state and effects characteristics
myAlg <- monanAlgorithmCreate(myState, myEffects, nsubN2 = 3,
                              multinomialProposal = FALSE)


##### estimate mobility network model #####

# mobility network model 
myResDN <- estimateMobilityNetwork(
  myState, myCache, myEffects, myAlg,
  initialParameters = NULL,
  # in case a pseudo-likelihood estimation was run, replace with
  # initialParameters = initEst,
  parallel = TRUE, cpus = 4,
  verbose = TRUE,
  returnDeps = TRUE,
  fish = FALSE
)

# check convergence
max(abs(myResDN$convergenceStatistics))

myResDN_old <- myResDN

# estimate mobility network model again based on previous results to improve convergence
# with an adjusted algorithm
myAlg <- monanAlgorithmCreate(myState, myEffects, multinomialProposal = TRUE, 
                              initialIterationsN2 = 100, nsubN2 = 1, initGain = 0.05, iterationsN3 = 1000)

# monan07 is an alias for estimateMobilityNetwork
myResDN <- monan07(
  myState, myCache, myEffects, myAlg,
  prevAns = myResDN,
  parallel = TRUE, cpus = 4,
  verbose = TRUE,
  returnDeps = TRUE,
  fish = FALSE
)

# check convergence
max(abs(myResDN$convergenceStatistics))

# view results
myResDN


##### regression diagnostics #####

autoCorrelationTest(myResDN)

traces <- extractTraces(myResDN, myEffects)
plot(traces)


##### test whether other effects should be included #####

myEffects2 <- createEffects(myState) |>
  addEffect(transitivity_min)

test_ME.2 <- scoreTest(myResDN, myEffects2)
test_ME.2


##### goodness of fit #####

myGofIndegree <- gofMobilityNetwork(ans = myResDN, gofFunction = getIndegree, lvls = 1:100)
plot(myGofIndegree,  lvls = 20:70)

myGofTieWeight <- gofMobilityNetwork(ans = myResDN, gofFunction = getTieWeights, lvls = 1:30)
plot(myGofTieWeight, lvls = 1:15)


##### simulate mobility network #####

mySimDN <- simulateMobilityNetworks(myState,
                                    myCache,
                                    myEffects,
                                    parameters = c(2, 1, 1.5, 0.1, -1, -0.5),
                                    allowLoops = TRUE,
                                    burnin = 45000,
                                    thinning = 15000,
                                    nSimulations = 10
)

mySimDN[[1]]
