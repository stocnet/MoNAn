#' @importFrom utils packageVersion packageDescription
.onAttach <- function(libname, pkgname) {
  if (!interactive()) return()
  base::packageStartupMessage(
    pkgname, " (Mobility Network Analysis): version ", 
    utils::packageVersion("MoNAn"), "\n", "\n",
    "Please note syntax changes from version 1.0.0 (released March 2024)."
  )
}