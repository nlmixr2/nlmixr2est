# nlmixr2est 7.1.0

## Test environments

* local: Ubuntu 24.04, R 4.6.1 (x86_64-pc-linux-gnu), gcc/g++ 14.2.0, checked
  against the CRAN versions of all dependencies (including `rxode2` 5.1.7)
* GitHub Actions: Windows, macOS and Linux, R release and devel

## R CMD check results

`R CMD check --as-cran` on the submitted tarball: 0 ERRORs, 0 WARNINGs,
3 NOTEs, none of which come from the package:

* `checking compilation flags used ... NOTE` (`-mno-omit-leaf-frame-pointer`):
  set by the local R installation's `Makeconf`, not by `src/Makevars`.
* `checking HTML version of manual ... NOTE`: HTML Tidy is not installed on
  the local machine.
* `checking examples ... NOTE` (examples > 5s): every example listed is inside
  `\donttest{}`; each fits a model and compiles it first.

## Reverse dependencies

All 16 CRAN reverse dependencies were checked (CRAN versions, CRAN-mode
`R CMD check`) against this version: `admixr2`, `babelmixr2`, `ggPMX`,
`nlmixr2`, `nlmixr2auto`, `nlmixr2autoinit`, `nlmixr2extra`, `nlmixr2lib`,
`nlmixr2plot`, `nlmixr2rpt`, `nlmixr2save`, `nlmixr2targets`, `nlmixr2utils`,
`shinyMixR`, `xpose.nlmixr2` and `xpose.xtras`.  All returned `Status: OK`.
