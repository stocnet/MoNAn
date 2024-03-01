####### wrapperFunctions


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
            state = deparse(substitute(state)))
  class(l) <- "effectsList.monan"
  l
}

#' addEffect
#' 
#' A function to add addtional effects to a moman effects object
#'
#' @param effectsObject The moman Effects object to which another
#' effect should be added.
#' @param effectName The name of the effect that should be added
#' as a character in quotes (e.g. "loops")
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
  
  new.eff <- createEffectsObject(effectInit = list(list(effectName, ...)),
                                 checkProcessState = eval(parse(text = eff.state)))
  effectsObject$effectFormulas[[(nEffects.prev + 1)]] <- new.eff$effectFormulas[[1]]
  if(new.eff$name[1] %in% effectsObject$name) stop(paste("Effect", new.eff$name[1], "already included in the effects object"))
  effectsObject$name[(nEffects.prev + 1)] <- new.eff$name[1]
  effectsObject
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
#'
#' @return An object of class "processState.monan".
#' @export
#'
#' @examples
#' monanDataCreate(transfers, people, organisations, 
#'                 sameRegion, region, size, sex)
monanDataCreate <- function(...){
  
  # create the named list to be transferred to createProcessState
  varnames <- unlist(lapply(substitute(list(...))[-1], deparse))
  elements <- list(...)
  names(elements) <- varnames
  
  # extract dependent variable
  dep.var <- NULL
  for (i in 1:length(elements)) {
    e <- elements[[i]]
    if (inherits(e, "edgelist.monan")){
      if(is.null(dep.var)) {
        dep.var <- names(elements)[i]
      } else {
        stop("More than one edgelist included; only one dependent variable possible")
      }
    }
  }
  if(is.null(dep.var)){
    stop("No dependent variable (edgelist) included in input objects")
  }
  
  # transfer both to createProcessState
  state <- createProcessState(elements, dep.var)
  state
}


