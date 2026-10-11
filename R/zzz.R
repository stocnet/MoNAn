#' @importFrom utils packageVersion packageDescription
.onAttach <- function(libname, pkgname) {
  if (!interactive()) return()
  base::packageStartupMessage(
    pkgname, " (Mobility Network Analysis): version ", 
    utils::packageVersion(pkgname), "\n", "\n",
    "Version 2 runs the simulations in C++ and corrects several effects;\n",
    "results can differ from version 1. See NEWS.md for all changes."
  )
}
