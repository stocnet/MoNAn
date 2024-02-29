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
#' myE1 <- createEffects(myState)
#' myE1 <- addEffect(myE1, "loops")
#' myE1 <- addEffect(myE1, "reciprocity_basic")
#' myE1 <- addEffect(myE1, effectName = "same_covariate", attribute.index = "region")
addEffect <- function(effectsObject, effectName, ...){
  eff.state <- effectsObject$state
  nEffects.prev <- length(effectsObject$effectFormulas)
  
  new.eff <- createEffectsObject(effectInit = list(list(effectName, ...)),
                                 checkProcessState = eval(parse(text = eff.state)))
  effectsObject$effectFormulas[[(nEffects.prev + 1)]] <- new.eff$effectFormulas[[1]]
  if(new.eff$name[1] %in% effectsObject$name) stop(paste("Effect", new.eff$name[1], "already included in the effects object"))
  effectsObject$name[(nEffects.prev + 1)] <- new.eff$name[1]
  effectsObject
}