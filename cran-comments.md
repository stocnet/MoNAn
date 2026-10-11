## Test environments

* local: macOS (aarch64), R 4.6.0 — `R CMD check --as-cran --run-donttest`
* win-builder: R-devel, R-release, R-oldrelease
* R-hub: linux, windows, macos, macos-arm64 (R-devel), clang-asan, clang-ubsan,
  gcc-asan, valgrind, gcc16 — all OK; no sanitizer or valgrind errors

## R CMD check results

0 errors | 0 warnings | 1 note

* Possibly misspelled words in DESCRIPTION: Geyer, Monro (see below).

## Comments

This is a major update of MoNAn (1.1.0 -> 2.0.0):

* The MCMC simulation now runs in compiled C++ code (using 'Rcpp'), so the
  package now contains compiled code.
* Several effects were corrected, which changes the results of models that use
  them, and a few effects were removed; all changes are listed in NEWS.md.
* There are no reverse dependencies on CRAN.
* Possibly misspelled words in DESCRIPTION: "Geyer" and "Monro" are names
  (Geyer-Thompson and Robbins-Monro algorithms); "analyse" and "modelling"
  are British spellings.
* The DOI link to the article in JRSSA in README.md
  (https://doi.org/10.1093/jrsssa/qnag057) may return 403 to automated checks
  because of the publisher's bot protection; it works in a browser.
