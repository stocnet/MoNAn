#' Example Data for the MoNAn Package
#'
#' @name mobilityData
#' @description
#' These are example data for the MoNAn package and can be used to estimate a 
#' mobility network. The raw example data is synthetic (i.e., made up). 
#' This fictitious example contains 17 organisations representing 
#' a labour market that are located in two regions (north and south). 
#' 742 workers are employed in these organisations at two time-points. 
#' Some are mobile while others work in the same organisation at both time-points. 
#' The following objects are provided for this purpose:
#'
#' \describe{
#'   \item{`mobilityEdgelist`}{The data frame indicates the origin at time 1 
#'   (first column) and the destination at time 2 (second column) for 
#'   each of the 742 individuals between the 17 organisations. Note that some
#'   workers stay in their organisation, i.e. their origin equals their 
#'   destination.}
#'   \item{`orgRegion`}{Categorical characteristic describing whether the 
#'   organisation is located on the northern (1) or southern (0) region.}
#'   \item{`orgSize`}{Continuous measure representing the size of each 
#'   organisation based on assets and revenue.}
#'   \item{`indSex`}{Individual-level characteristics representing sex.}
#' }
#' 
#' @rdname mobilityData
#' @usage NULL
#' @format `mobilityEdgelist`
#' A data frame with 743 rows and 2 columns.
"mobilityEdgelist"


#' @rdname mobilityData
#' @usage NULL
#' @format `orgRegion`
#' An object with 17 values.
"orgRegion"


#' @rdname mobilityData
#' @usage NULL
#' @format `orgSize`
#' An object with 17 values.
"orgSize"


#' @rdname mobilityData
#' @usage NULL
#' @format `indSex`
#' An object with 742 values.
"indSex"
