# nlmixr2est 7.2.1

This release fixes the issues flagged in the CRAN checks of 7.2.0:

* r-release-macos-x86_64, r-oldrel-macos-x86_64 and linux-arm64: the
  `nlmixr2GradFun()` example stopped with "one or more parameters have a zero
  gradient".  The Gill (1983) step search estimated the error in an objective
  near 0 as essentially 0, shrinking the finite-difference steps to roundoff;
  that error is now floored, and the example's derivatives are accurate on every
  platform.  A search ending in its "odd or nearly linear" case (which the
  example reaches on macOS) also reported the slope at its largest trial step;
  it now uses its starting step.
* noLD: a test compared a C++ BFGS update bitwise with an R reference; it now
  compares to a tolerance of 1e-10.

The first submission of 7.2.1 gave "Overall checktime 13 min > 10 min" on
r-devel-windows-x86_64.  The test files that compile models no longer run on
CRAN (they still run on our CI), and the events vignette's SAEM fit is shorter.

## Test environments

* local: Ubuntu 24.04, R 4.6.1 (x86_64-pc-linux-gnu), gcc/g++ 14.2.0, checked
  against the CRAN versions of all dependencies (including `rxode2` 5.1.8)
* GitHub Actions: Windows, macOS and Linux, R release and devel

## R CMD check results

`R CMD check --as-cran` on the submitted tarball: 0 ERRORs, 0 WARNINGs, 2 NOTEs:

* `checking CRAN incoming feasibility ... NOTE` (Days since last update: 1):
  this release fixes the check issues flagged on 7.2.0 above.
* `checking compilation flags used ... NOTE` (`-mno-omit-leaf-frame-pointer`):
  set by the local R installation's `Makeconf`, not by `src/Makevars`.

## Reverse dependencies

All 18 CRAN reverse dependencies were checked (CRAN versions, CRAN-mode
`R CMD check`) against this version: `admixr2`, `babelmixr2`, `ggPMX`,
`nlmixr2`, `nlmixr2auto`, `nlmixr2autoinit`, `nlmixr2extra`, `nlmixr2lib`,
`nlmixr2plot`, `nlmixr2rpt`, `nlmixr2save`, `nlmixr2scm`, `nlmixr2targets`,
`nlmixr2utils`, `nonmem2rx`, `shinyMixR`, `xpose.nlmixr2` and `xpose.xtras`.
All returned `Status: OK`.
