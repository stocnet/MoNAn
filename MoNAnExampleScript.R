# This script runs a simple example with the data from the MoNAn package

library(MoNAn)


##### create data objects from internal data files, 
# which are later combined to the process state #####

# extract number of individuals and organisations from the mobility data

N_ind <- nrow(mobilityEdgelist)
N_org <- length(unique(as.numeric(mobilityEdgelist)))

# Create a process state out of the mobility data objects:
# create objects (which are later combined to the process state)
transfers <- monanDependent(mobilityEdgelist,
                            nodes = "organisations",
                            edges = "people")

people <- monanEdges(N_ind)
organisations <- monanNodes(N_org)

sameRegion <- outer(orgRegion, orgRegion, "==") * 1
sameRegion <- dyadicCovar(sameRegion, nodes = "organisations")

region <- monadicCovar(orgRegion, nodes = "organisations")
size <- monadicCovar(orgSize, nodes = "organisations")
sex <- monadicCovar(indSex, edges = "people")


# the following lines create an artificial second origin used for illustration
# in the examples in the manual
other_origin <- sample(1:17, 742, replace = T)
resample <- as.logical(sample(0:1, 742, replace = T, prob = c(0.88, 0.12)))
other_origin[resample] <- transfers$data[resample,2]
second_or <- monadicCovar(other_origin, edges = "people")


# combine created objects to the process state
myState <- monanDataCreate(transfers,
                           people,
                           organisations,
                           sameRegion,
                           region,
                           size,
                           sex,
                           second_or,
                           fixedEffectDummies = TRUE)

# inspect the created object
myState

##### create effects object #####

# effects object
myEffects <- createEffects(myState) |>
  addEffect(loops) |>
  addEffect(concentration_AC, alpha = 4) |>
  addEffect(reciprocity_AC, alpha = 4) |>
  addEffect(dyadic_covariate, node.attribute = "sameRegion") |>
  addEffect(alter_covariate, node.attribute = "size") |>
  addEffect(resource_covar_to_node_covar,
            node.attribute = "region",
            edge.attribute = "sex") |>
  addEffect(loops_resource_covar, edge.attribute = "sex")

# inspect the created object
myEffects

# further effects object with fixed effects by location

myEffects_fe <- createEffects(myState) |>
  addEffect(loops) |>
  addEffect(concentration_AC, alpha = 4) |>
  addEffect(reciprocity_AC, alpha = 4) |>
  addEffect(dyadic_covariate, node.attribute = "sameRegion") |>
  addEffect(resource_covar_to_node_covar,
            node.attribute = "region",
            edge.attribute = "sex") |>
  addEffect(loops_resource_covar, edge.attribute = "sex") |>
  addFixedEffects(myState)

myEffects_fe

##### get multinomial statistics to estimate initial parameters using pseudo-likelihood estimation #####

# create statistics object, to be used, e.g., with the mlogit package
myStatisticsFrame <- getMultinomialStatistics(myState, myEffects)

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
myResDN <- monanEstimate(
  myState, myEffects, myAlg,
  initialParameters = NULL,
  # in case a pseudo-likelihood estimation was run, replace with
  # initialParameters = initEst,
  parallel = TRUE, cpus = 4,
  verbose = TRUE,
  returnDeps = TRUE,
  fish = FALSE
)

# myResDN_fe <- monanEstimate(
#   myState, myEffects_fe, myAlg,
#   initialParameters = NULL,
#   # in case a pseudo-likelihood estimation was run, replace with
#   # initialParameters = initEst,
#   parallel = TRUE, cpus = 4,
#   verbose = TRUE,
#   returnDeps = TRUE,
#   fish = FALSE
# )

# check convergence
max(abs(myResDN$convergenceStatistics))

myResDN_old <- myResDN

# estimate mobility network model again based on previous results to improve convergence
# with an adjusted algorithm
myAlg <- monanAlgorithmCreate(myState, myEffects, multinomialProposal = TRUE, 
                              initialIterationsN2 = 100, nsubN2 = 1, initGain = 0.05, iterationsN3 = 1000)

# monan07 is an alias for monanEstimate
myResDN <- monan07(
  myState, myEffects, myAlg,
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
  addEffect(transitivity_AC)

test_ME.2 <- scoreTest(myResDN, myEffects2)
test_ME.2


##### goodness of fit #####

myGofIndegree <- monanGOF(ans = myResDN, 
                          gofFunction = getIndegree, 
                          lvls = 1:100)
plot(myGofIndegree,  lvls = 20:70)

myGofTieWeight <- monanGOF(ans = myResDN, 
                           gofFunction = getTieWeights, 
                           lvls = 1:30)
plot(myGofTieWeight, lvls = 1:15)


##### simulate mobility network #####

mySimDN <- monanSimulate(myState,
                         myEffects,
                         parameters = c(2, 1, 1.5, 0.5, 0.1, -1, -0.5),
                         allowLoops = TRUE,
                         burnin = 45000,
                         thinning = 15000,
                         nSimulations = 10
)

mySimDN[[1]]
