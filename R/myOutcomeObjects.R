#' Exemplary Outcome Objects for the MoNAn Package
#'
#' @name myOutcomeObjects
#' @description
#' These are exemplary outcome objects for the MoNAn package and can be used in 
#' order not to run all precedent functions and thus save time. The following 
#' products are provided:
#' 
#' @rdname myOutcomeObjects
#' @format `myState`
#' An object of class `processState.monan` created by the function [createProcessState()].
"myState"


#' @rdname myOutcomeObjects
#' @format `myDependentVariable`
#' An object of class `character`.
"myDependentVariable"


#' @rdname myOutcomeObjects
#' @format `myCache`
#' An object of class `list` created by the function [createWeightedCache()].
"myCache"


#' @rdname myOutcomeObjects
#' @format `myEffects`
#' An object of class `effectsList.monan` created by the function [createEffectsObject()].
"myEffects"


#' @rdname myOutcomeObjects
#' @format `myAlg`
#' An object of class `algorithm.monan` created by the function [createAlgorithm()].
"myAlg"


#' @rdname myOutcomeObjects
#' @format `myResDN`
#' An object of class `result.monan` created by the function [estimateMobilityNetwork()].
"myResDN"


#' @rdname myOutcomeObjects
#' @format `mySimDN`
#' An object of class `sims.monan` created by the function [simulateMobilityNetworks()].
"mySimDN"
