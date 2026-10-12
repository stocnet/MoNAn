## Test environments

* local: macOS (aarch64), R 4.6.0 — `R CMD check --as-cran --run-donttest`
* win-builder: R-devel, R-release, R-oldrelease
* R-hub: linux, windows, macos, macos-arm64 (R-devel), clang-asan, clang-ubsan,
  gcc-asan, valgrind, gcc16, nold

## R CMD check results

0 errors | 0 warnings | 1 note

* Days since last update: 1. This update fixes the test failures of 2.0.0 in
  the 'noLD' additional check, as requested in the email from Brian Ripley of
  12 October 2026 (see below), so it follows 2.0.0 quickly.

## Comments

This is an update of MoNAn (2.0.0 -> 2.1.1; 2.1.0 was not submitted).

* In response to the email from Brian Ripley of 12 October 2026, it fixes the
  test failures in the 'noLD' additional check for 2.0.0: the
  tests compare the C++ and R versions of the model's effects, which are
  identical when R sums in long double but differ in the last digits without
  long double. The tests now allow a relative difference of 1e-10 when
  `capabilities("long.double")` is FALSE. The R-hub 'nold' check passes.
* Parallel computation now uses the base package 'parallel' instead of
  'snowfall'. Examples use at most 2 cores.
* Further changes (a corrected score test, restored effects, bug fixes) are
  listed in NEWS.md.
* There are no reverse dependencies on CRAN.
* "Geyer" and "Monro" are names (Geyer-Thompson and Robbins-Monro algorithms);
  "analyse" and "modelling" are British spellings.
* The DOI link to the article in JRSSA in README.md
  (https://doi.org/10.1093/jrsssa/qnag057) may return 403 to automated checks
  because of the publisher's bot protection; it works in a browser.
