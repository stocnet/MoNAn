########## effects test: a simple test for newly programmed effects


#' test_effect
#'  
#' Test for each person in the example data and one randomly selected
#' alternative destination whether change and target statistics match.
#'
#' @param effectName 
#' @param ... 
#'
#' @return A data frame containing for each person in the example data:
#' Origin; Destination; Suggested destination; Change statistic;
#' Target before change; Target after change; Whether this is an exact match;
#' Whether this is a match rounded to 8 digits.
#' In case the last column contains values that are not all "1",
#' the tested function does not work properly.
#'
#' @keywords internal
#' @examples
#' test_loop <- test_effect("loops")
test_effect <- function(effectName, ...){
  docu <- data.frame(orig = 0, dest = 0, new_dest = 0, 
                     change = 0, target_before = 0, target_after = 0, 
                     exact_match = 0, approx_match = 0)
  
  for(edge in 1:742){
    sug_dest <- sample(1:17, 1)
    x <- target_change_match(edge = edge, j_new = sug_dest, effectName, ...)
    if(x[4] !=1 ){
      paste("edge = ", edge, ", sug_dest = ", sug_dest, ": NO PRECISE MATCH")
    }
    docu[edge,] <- c(myState[[1]]$data[edge,1], myState[[1]]$data[edge,2], sug_dest, x)
  }
  if(all(docu[,7] == 1)){
    print("No mis-match found")
  }
  return(docu)
}

#' target_change_match
#' 
#' Function to test effects
#' 
#' Tests for specific resources and new proposed destination whether
#' the change statistic and the difference between two target statistics 
#' before and after change match
#'
#' @param edge 
#' @param j_new 
#' @param effectName 
#' @param ... 
#'
#' @return a vector containing (the change stat),(the target stat before change)
#' (the target stat after change),(whether this matches exactly),
#' (whether it matches rounded)
#'
#' @keywords internal
target_change_match <- function(edge, j_new, effectName, ...){
  i_test <- myState[[1]]$data[edge,1]
  j_test <- myState[[1]]$data[edge,2]
  
  count1 <- 0
  
  for(i in 1:17){
    for(j in 1:17){
      count1 <- count1 + match.fun(effectName)(dep.var = 1,
                                               state = myState, 
                                               cache = myCache, 
                                               i = i, j = j, 
                                               edge = 1, 
                                               update = 1, 
                                               getTargetContribution = T,
                                               ...)
    }
  }
  
  count1

  change <- match.fun(effectName)(dep.var = 1,
                                  state = myState, 
                                  cache = myCache, 
                                  i = i_test, j = j_test, 
                                  edge = edge, 
                                  update = -1,
                                  getTargetContribution = F,
                                  ...)
  myUpdtCache <- myCache
  myUpdtCache[[1]] <-
    updateWeightedCache(
      myCache[[1]],
      i_test,
      j_test,
      resourceID = edge,
      state = myState,
      dep.var = 1,
      update = -1
    )
  
  change <- change + match.fun(effectName)(dep.var = 1,
                                           state = myState, 
                                           cache = myUpdtCache, 
                                           i = i_test, j = j_new, 
                                           edge = edge, 
                                           update = 1,
                                           getTargetContribution = F,
                                           ...)
  
  myNewState <- myState
  
  myNewState[[1]]$data[edge,2] <- j_new
  
  # create cache
  myNewCache <- createWeightedCache(myNewState, resourceCovariates = c("sex"))
  
  count2 <- 0
  
  for(i in 1:17){
    for(j in 1:17){
      count2 <- count2 + match.fun(effectName)(dep.var = 1,
                                               state = myNewState, 
                                               cache = myNewCache, 
                                               i = i, j = j, 
                                               edge = 1, 
                                               update = 1, 
                                               getTargetContribution = T,
                                               ...)
    }
  }
  
  count3 <- (count1 + change)
  
  return(c(change, count1, count2, (count3 == count2), (round(count2, 8) == round(count3, 8))))
}


