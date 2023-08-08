#' Example Data for the MoNAn Package
#'
#' @name mobilityData
#' @description
#' These are example data for the MoNAn package and can be used to estimate a 
#' mobility network. The following objects are provided for this purpose:
#'
#' \describe{
#'   \item{`mobilityEdgelist`}{This data frame represents the origin at time 1 
#'   (first column) as well as the destination at time 2 (second column) for 
#'   each of the 742 individuals between the 17 organisations. Individuals could
#'   also stay in their organisation meaning that origin equals destination.}
#'   \item{`orgRegion`}{This characteristic describes whether the resp. 
#'   organisation is located on the northern (1) or southern (0) island.}
#'   \item{`orgSize`}{This composite measure represents the size of each 
#'   organisation including its number of workers, assets, and revenue.}
#'   \item{`indSex`}{This characteristic represents each individual's gender.}
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
