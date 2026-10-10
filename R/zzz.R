#' @importFrom utils packageVersion packageDescription
.onAttach <- function(libname, pkgname) {
  if (!interactive()) return()
  base::packageStartupMessage(
    pkgname, " (Mobility Network Analysis 2): version ", 
    utils::packageVersion(pkgname), "\n", "\n",
    "monan2 is the successor of MoNAn"
  )
}
