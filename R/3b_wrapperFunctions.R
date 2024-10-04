####### wrapperFunctions


#' addDummyEffects
#' 
#' A function to add one [alter_covariate] effect for each location
#' dummy variable included in the state. Dummy variables in the state object
#' are named "dummyXX", where XX can take any value.
#' This is particularly useful when including fixed effects for each location
#' to explicitly model the indegree of each location.
#'
#' @param effectsObject The monan Effects object to which the dummy
#' location effects should be added.
#' @param state The monan state (data) object that the effects object applies to
#'
#' @return An object of type effectsList.monan
#' @export
#'
#' @examples
#' myEffects <- createEffects(myState) |>
#'   addEffect(loops) |>
#'   addDummyEffects(myState)
addDummyEffects <- function(effectsObject, state){
  stateName <- deparse(substitute(state))
  if(stateName != effectsObject$state){
    stop("in addDummyEffects: this is not the state specified in the effectsObject")
  }
  dummies_to_include <- names(state)[grep("dummy", names(state))]
  for (i in dummies_to_include){
    effectsObject <- addEffect(effectsObject, 
                               alter_covariate, 
                               node.attribute = i)
  }
  effectsObject
}

#' addFixedEffects
#'
#' @rdname addDummyEffects
addFixedEffects <- addDummyEffects


#' addEffect
#' 
#' A function to add additional effects to a moman effects object
#'
#' @param effectsObject The monan Effects object to which another
#' effect should be added.
#' @param effectName The name of the effect that should be added (e.g. loops).
#' @param ... Additional parameters of the effect, for example alpha,
#' attribute.index, or resource.attribute.index
#'
#' @return An object of type effectsList.monan
#' @export
#'
#' @examples
#' # Create effects object and add effects
#' myE1 <- createEffects(myState)
#' myE1 <- addEffect(myE1, loops)
#' myE1 <- addEffect(myE1, reciprocity_basic)
#' myE1 <- addEffect(myE1, effectName = same_covariate, attribute.index = "region")
#' 
#' # Or simpler
#' myE1 <- createEffects(myState) |>
#'   addEffect(loops) |>
#'   addEffect(reciprocity_basic) |>
#'   addEffect(same_covariate, attribute.index = "region")
addEffect <- function(effectsObject, effectName, ...){
  effectName <- deparse(substitute(effectName))
  eff.state <- effectsObject$state
  nEffects.prev <- length(effectsObject$effectFormulas)
  
  dots <- list(...)
  effect_arg <- c(effectName, dots)
  
  names(effect_arg)[names(effect_arg) == "node.attribute"] <- "attribute.index"
  names(effect_arg)[names(effect_arg) == "edge.attribute"] <- "resource.attribute.index"
  
  new.eff <- do.call(createEffectsObject, 
                     list( effectInit = list(effect_arg) , 
                           checkProcessState = eval(parse(text = eff.state)) ))
  
  effectsObject$effectFormulas[[(nEffects.prev + 1)]] <- new.eff$effectFormulas[[1]]
  if(new.eff$name[1] %in% effectsObject$name) stop(paste("Effect", new.eff$name[1], "already included in the effects object"))
  effectsObject$name[(nEffects.prev + 1)] <- new.eff$name[1]
  effectsObject
}



#' createEffects
#'
#' Generates an empty effects object to which new effects can be added 
#' consecutively
#'
#' @param state The state to which the model applies.
#'
#' @return An empty effects object of class effectsList.monan
#' @export
#'
#' @examples
#' #' myE1 <- createEffects(myState)
createEffects <- function(state){
  l <- list(effectFormulas = list(),
            name = c(),
            state = deparse(substitute(state)),
            nodeset1 = state[[state$dep.var]]$nodeSet[1],
            nodeset2 = state[[state$dep.var]]$nodeSet[3]
            )
  class(l) <- "effectsList.monan"
  l
}



#' monanDataCreate
#'
#' A function to create a moman process state, i.e., 
#' a MoNAn object that stores all information
#' about the data that will be used in the estimation. This includes the
#' outcome variable (edgelist), the nodesets, and all covariates.
#'
#' @param ... The monan objects to be included in the process State. 
#' This must include exactly one edgelist (dependent variable) and the
#' two nodesets associated with the edgelist.
#' Further allowed elements are (monadic or dyadic) covariates of locations and people
#' @param fixedEffectsDummies Add dummy variable for each location (minus the first)
#' in the process state (data object), so that fixed effects to control the
#' indegree of every location can be modelled directly.
#' 
#' @return An object of class "processState.monan".
#' @export
#'
#' @examples
#' #' # create objects (which are later combined to the process state)
#' transfers <- createEdgelist(mobilityEdgelist,
#'   nodeSet = c("organisations", "organisations", "people")
#' )
#' people <- createNodeSet(1:nrow(mobilityEdgelist))
#' organisations <- createNodeSet(1:length(orgRegion))
#' sameRegion <- outer(orgRegion, orgRegion, "==") * 1
#' sameRegion <- createNetwork(sameRegion,
#'   nodeSet = c("organisations", "organisations")
#' )
#' region <- createNodeVariable(orgRegion, nodeSet = "organisations")
#' size <- createNodeVariable(orgSize, nodeSet = "organisations", addSim = TRUE)
#' sex <- createNodeVariable(indSex, nodeSet = "people")
#' 
#' monanDataCreate(transfers, people, organisations, 
#'                 sameRegion, region, size, sex)
monanDataCreate <- function(..., fixedEffectDummies = FALSE){
  varnames <- unlist(lapply(substitute(list(...))[-1], deparse))
  elements <- list(...)
  names(elements) <- varnames
  dep.var <- NULL
  for (i in 1:length(elements)) {
    e <- elements[[i]]
    if (inherits(e, "edgelist.monan")) {
      if (is.null(dep.var)) {
        dep.var <- names(elements)[i]
        edgelist <- e
      }
      else {
        stop("More than one edgelist included; only one dependent variable possible")
      }
    }
  }
  if (is.null(dep.var)) {
    stop("No dependent variable (edgelist) included in input objects")
  }
  
  if (!fixedEffectDummies){
    state <- createProcessState(elements, dep.var)
  } else {
    
    createDummy <- function(length, position){
      vec <- rep(0, length)
      vec[position] <- 1
      vec
    }
    
    length <- length(elements[[edgelist$nodeSet[1] ]]$ids)
    dummies <- lapply(2:length, function(i) createDummy(length, i))
    
    fixedEffects <- lapply(dummies, monadicCovar, 
                           nodes = edgelist$nodeSet[1])
    
    dummy_name <- paste0("dummy", 2:length)
    names(fixedEffects) <- dummy_name
    
    elements_FE <- c(elements, fixedEffects)
    state <- createProcessState(elements_FE, dep.var)
    
  }
  
  state
}


